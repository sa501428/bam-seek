#include <bamseek/evidence.hpp>

#include <algorithm>
#include <cctype>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace bamseek {
namespace {

constexpr std::uint16_t flag_unmapped = 0x4;
constexpr std::uint16_t flag_reverse = 0x10;
constexpr std::uint16_t flag_secondary = 0x100;
constexpr std::uint16_t flag_duplicate = 0x400;
constexpr std::uint16_t flag_supplementary = 0x800;

bool allowed_resource_uri(const std::string& uri) {
    return uri.find("://") == std::string::npos || uri.starts_with("https://");
}

bool allowed_resource_uri(const std::optional<std::string>& uri) {
    return !uri || allowed_resource_uri(*uri);
}

std::unique_ptr<igv::AlignmentReader> open_allowed_alignments(const igv::Resource& resource) {
    if (!allowed_resource_uri(resource.uri) || !allowed_resource_uri(resource.index_uri) || !allowed_resource_uri(resource.reference_uri)) {
        throw std::runtime_error("BAM Seek accepts local paths and HTTPS resources only; HTTP and other URL schemes are not permitted.");
    }
    return igv::open_alignments(resource);
}

struct ReadLayout {
    std::map<std::int64_t, std::pair<char, int>> bases;
    std::map<std::int64_t, std::string> insertions_after;
    std::vector<igv::GenomicInterval> deletions;
};

ReadLayout layout_for(const igv::Alignment& alignment) {
    ReadLayout layout;
    auto reference_position = alignment.interval.start;
    std::size_t read_position = 0;
    for (const auto& operation : alignment.cigar) {
        const auto length = static_cast<std::size_t>(operation.length);
        switch (operation.operation) {
            case 'M': case '=': case 'X':
                for (std::size_t i = 0; i < length && read_position + i < alignment.sequence.size(); ++i) {
                    const auto quality = read_position + i < alignment.qualities.size()
                        ? static_cast<int>(static_cast<unsigned char>(alignment.qualities[read_position + i])) - 33
                        : 0;
                    layout.bases.emplace(reference_position + static_cast<std::int64_t>(i),
                                         std::pair{static_cast<char>(std::toupper(static_cast<unsigned char>(alignment.sequence[read_position + i]))), quality});
                }
                reference_position += static_cast<std::int64_t>(length);
                read_position += length;
                break;
            case 'I':
                if (read_position < alignment.sequence.size() && reference_position > alignment.interval.start) {
                    layout.insertions_after[reference_position - 1] = alignment.sequence.substr(read_position, length);
                    auto& insertion = layout.insertions_after[reference_position - 1];
                    std::transform(insertion.begin(), insertion.end(), insertion.begin(), [](unsigned char c) {
                        return static_cast<char>(std::toupper(c));
                    });
                }
                read_position += length;
                break;
            case 'D': case 'N':
                layout.deletions.push_back({alignment.interval.contig, reference_position,
                                            reference_position + static_cast<std::int64_t>(length)});
                reference_position += static_cast<std::int64_t>(length);
                break;
            case 'S':
                read_position += length;
                break;
            case 'H': case 'P':
                break;
            default:
                break;
        }
    }
    return layout;
}

bool deleted_at(const ReadLayout& layout, const std::int64_t position) {
    return std::any_of(layout.deletions.begin(), layout.deletions.end(), [position](const auto& deletion) {
        return deletion.start <= position && position < deletion.end;
    });
}

std::optional<int> sequence_matches(const ReadLayout& layout, const std::int64_t start, const std::string& sequence) {
    int minimum_quality = 10000;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const auto position = start + static_cast<std::int64_t>(i);
        const auto base = layout.bases.find(position);
        if (base == layout.bases.end() || deleted_at(layout, position) || base->second.first != sequence[i]) return std::nullopt;
        minimum_quality = std::min(minimum_quality, base->second.second);
    }
    return minimum_quality == 10000 ? std::optional<int>{} : minimum_quality;
}

struct AlleleCall { Allele allele; int minimum_base_quality; };

AlleleCall call_allele(const igv::Alignment& alignment, const VariantQuery& query) {
    const auto layout = layout_for(alignment);
    if (!layout.bases.contains(query.position) && !deleted_at(layout, query.position)) return {Allele::no_call, 0};
    const auto ref_match = sequence_matches(layout, query.position, query.reference);
    const auto alt_match = sequence_matches(layout, query.position, query.alternate);
    const auto anchor = query.position + static_cast<std::int64_t>(query.reference.size()) - 1;

    if (query.reference.size() == query.alternate.size()) {
        if (alt_match) return {Allele::alternate, *alt_match};
        if (ref_match) return {Allele::reference, *ref_match};
        return {Allele::other, 0};
    }
    if (query.alternate.size() > query.reference.size()) {
        const auto insertion = layout.insertions_after.find(anchor);
        const auto expected = query.alternate.substr(query.reference.size());
        if (ref_match && insertion != layout.insertions_after.end() && insertion->second == expected) return {Allele::alternate, *ref_match};
        if (ref_match && insertion == layout.insertions_after.end()) return {Allele::reference, *ref_match};
        return {Allele::other, ref_match.value_or(0)};
    }
    const auto deletion_start = query.position + static_cast<std::int64_t>(query.alternate.size());
    const auto deletion_end = query.position + static_cast<std::int64_t>(query.reference.size());
    const bool has_expected_deletion = std::any_of(layout.deletions.begin(), layout.deletions.end(), [&](const auto& deletion) {
        return deletion.start <= deletion_start && deletion.end >= deletion_end;
    });
    const auto prefix = sequence_matches(layout, query.position, query.alternate);
    if (prefix && has_expected_deletion) {
        return {Allele::alternate, *prefix};
    }
    if (ref_match) return {Allele::reference, *ref_match};
    return {Allele::other, ref_match.value_or(0)};
}

bool include_alignment(const igv::Alignment& alignment, const FilterSettings& filters) {
    if ((alignment.flags & flag_unmapped) != 0) return false;
    if (!filters.include_duplicates && (alignment.flags & flag_duplicate) != 0) return false;
    if (!filters.include_secondary && (alignment.flags & flag_secondary) != 0) return false;
    if (!filters.include_supplementary && (alignment.flags & flag_supplementary) != 0) return false;
    return alignment.mapping_quality >= filters.minimum_mapping_quality;
}

std::optional<std::string> tag_value(const igv::Alignment& alignment, const std::string& name) {
    const auto found = std::find_if(alignment.tags.begin(), alignment.tags.end(), [&](const auto& tag) { return tag.name == name; });
    if (found == alignment.tags.end() || found->value.empty()) return std::nullopt;
    return found->value;
}

std::string choose_tag(const std::vector<igv::Alignment>& alignments, const FilterSettings& filters) {
    if (filters.molecule_mode == MoleculeMode::raw_reads) return {};
    if (filters.molecule_mode == MoleculeMode::selected_tag) return filters.molecule_tag;
    for (const std::string candidate : {"MI", "RX", "UB"}) {
        if (std::any_of(alignments.begin(), alignments.end(), [&](const auto& alignment) { return tag_value(alignment, candidate).has_value(); })) return candidate;
    }
    return {};
}

void add_count(EvidenceCounts& counts, const Allele allele, const bool reverse) {
    switch (allele) {
        case Allele::reference: ++counts.reference_reads; break;
        case Allele::alternate:
            ++counts.alternate_reads;
            if (reverse) ++counts.alternate_reverse_reads; else ++counts.alternate_forward_reads;
            break;
        case Allele::other: ++counts.other_reads; break;
        case Allele::no_call: break;
    }
}

void add_molecule_count(EvidenceCounts& counts, const Allele allele) {
    switch (allele) {
        case Allele::reference: ++counts.reference_molecules; break;
        case Allele::alternate: ++counts.alternate_molecules; break;
        case Allele::other: case Allele::no_call: ++counts.other_molecules; break;
    }
}

VariantEvidence evaluate_variant(const igv::AlignmentReader& reader, const VariantQuery& variant, const FilterSettings& filters) {
    VariantEvidence evidence;
    evidence.query = variant;
    const auto alignments = reader.get(variant.query_window());
    evidence.molecule_tag_used = choose_tag(alignments, filters);
    std::unordered_map<std::string, std::vector<Allele>> molecule_calls;
    for (const auto& alignment : alignments) {
        if (!include_alignment(alignment, filters)) continue;
        const auto called = call_allele(alignment, variant);
        if (called.allele == Allele::no_call || called.minimum_base_quality < filters.minimum_base_quality) continue;
        const bool reverse = (alignment.flags & flag_reverse) != 0;
        std::string molecule_id;
        if (!evidence.molecule_tag_used.empty()) molecule_id = tag_value(alignment, evidence.molecule_tag_used).value_or("read:" + alignment.name);
        std::ostringstream summary;
        summary << allele_name(called.allele) << " mapQ=" << static_cast<int>(alignment.mapping_quality)
                << " baseQ=" << called.minimum_base_quality << " CIGAR=";
        for (const auto& operation : alignment.cigar) summary << operation.length << operation.operation;
        evidence.reads.push_back({alignment.name, called.allele, reverse, alignment.mapping_quality,
                                  called.minimum_base_quality, molecule_id, summary.str()});
        add_count(evidence.counts, called.allele, reverse);
        molecule_calls[molecule_id.empty() ? "read:" + alignment.name : molecule_id].push_back(called.allele);
    }
    for (const auto& [identifier, calls] : molecule_calls) {
        (void)identifier;
        const auto alt = static_cast<int>(std::count(calls.begin(), calls.end(), Allele::alternate));
        const auto ref = static_cast<int>(std::count(calls.begin(), calls.end(), Allele::reference));
        add_molecule_count(evidence.counts, alt > ref ? Allele::alternate : ref > alt ? Allele::reference : Allele::other);
    }
    evidence.passes_thresholds = evidence.counts.alternate_reads >= filters.minimum_alternate_reads
        && evidence.counts.allele_fraction() >= filters.minimum_variant_allele_fraction
        && (evidence.molecule_tag_used.empty() || evidence.counts.alternate_molecules >= filters.minimum_alternate_molecules);
    return evidence;
}

}  // namespace

double EvidenceCounts::allele_fraction() const noexcept {
    return depth() == 0 ? 0.0 : static_cast<double>(alternate_reads) / static_cast<double>(depth());
}

std::string allele_name(const Allele allele) {
    switch (allele) {
        case Allele::reference: return "REF";
        case Allele::alternate: return "ALT";
        case Allele::other: return "OTHER";
        case Allele::no_call: return "NO_CALL";
    }
    return "UNKNOWN";
}

EvidenceEngine::EvidenceEngine(igv::Resource resource) : reader_(open_allowed_alignments(resource)) {
    if (!reader_->indexed()) throw std::runtime_error("An indexed BAM, CRAM, or SAM resource is required (.bai, .csi, or .crai).");
    if (resource.reference_uri) reference_ = igv::open_fasta({.uri = *resource.reference_uri});
}

bool EvidenceEngine::indexed() const noexcept { return reader_->indexed(); }

PileupData EvidenceEngine::pileup(const VariantQuery& query, const FilterSettings& filters, const std::int64_t padding) const {
    if (padding < 0 || padding > 1000) throw std::invalid_argument("Pileup padding must be between 0 and 1,000 bases.");
    PileupData data;
    data.query = query;
    data.interval = query.query_window(padding);
    if (reference_) {
        data.reference_bases = reference_->get(data.interval);
        data.has_reference = true;
    }
    for (const auto& alignment : reader_->get(data.interval)) {
        if (include_alignment(alignment, filters)) data.alignments.push_back(alignment);
    }
    return data;
}

BatchEvidence EvidenceEngine::evaluate(const std::vector<Query>& queries, const FilterSettings& filters) const {
    BatchEvidence batch;
    for (const auto& query : queries) {
        if (const auto* region = std::get_if<RegionQuery>(&query)) {
            if (!reference_) {
                batch.results.emplace_back(RegionEvidence{*region, "Region scanning requires an indexed hg19 FASTA reference.", {}});
                continue;
            }
            if (region->interval.end - region->interval.start > 100000) {
                batch.results.emplace_back(RegionEvidence{*region, "Region is larger than 100 kb. Split it into smaller windows for targeted scanning.", {}});
                continue;
            }
            const auto reference_bases = reference_->get(region->interval);
            std::map<std::int64_t, std::map<char, int>> alternate_observations;
            const auto alignments = reader_->get(region->interval);
            for (const auto& alignment : alignments) {
                if (!include_alignment(alignment, filters)) continue;
                const auto layout = layout_for(alignment);
                for (const auto& [position, base] : layout.bases) {
                    if (position < region->interval.start || position >= region->interval.end || base.second < filters.minimum_base_quality) continue;
                    const auto offset = static_cast<std::size_t>(position - region->interval.start);
                    if (offset >= reference_bases.size()) continue;
                    const auto reference_base = static_cast<char>(std::toupper(static_cast<unsigned char>(reference_bases[offset])));
                    if (base.first != reference_base && reference_base != 'N' && base.first != 'N') ++alternate_observations[position][base.first];
                }
            }
            RegionEvidence region_evidence{*region, "No candidates met the selected thresholds.", {}};
            for (const auto& [position, observed] : alternate_observations) {
                const auto offset = static_cast<std::size_t>(position - region->interval.start);
                const auto reference_base = static_cast<char>(std::toupper(static_cast<unsigned char>(reference_bases[offset])));
                for (const auto& [alternate_base, preliminary_count] : observed) {
                    if (preliminary_count < filters.minimum_alternate_reads) continue;
                    auto candidate = evaluate_variant(*reader_, {region->source_text, region->interval.contig, position,
                                                                  std::string(1, reference_base), std::string(1, alternate_base)}, filters);
                    if (candidate.passes_thresholds) region_evidence.candidates.push_back(std::move(candidate));
                }
            }
            if (!region_evidence.candidates.empty()) region_evidence.note = std::to_string(region_evidence.candidates.size()) + " SNV candidate(s) met thresholds.";
            batch.results.emplace_back(std::move(region_evidence));
            continue;
        }
        const auto& variant = std::get<VariantQuery>(query);
        try {
            batch.results.emplace_back(evaluate_variant(*reader_, variant, filters));
        } catch (const std::exception& error) {
            batch.errors.push_back(variant.source_text + ": " + error.what());
        }
    }
    return batch;
}

}  // namespace bamseek
