#include <bamseek/evidence.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace bamseek {
namespace {

constexpr std::uint16_t flag_unmapped = 0x4;
constexpr std::uint16_t flag_reverse = 0x10;
constexpr std::uint16_t flag_paired = 0x1;
constexpr std::uint16_t flag_proper_pair = 0x2;
constexpr std::uint16_t flag_first_mate = 0x40;
constexpr std::uint16_t flag_second_mate = 0x80;
constexpr std::uint16_t flag_secondary = 0x100;
constexpr std::uint16_t flag_duplicate = 0x400;
constexpr std::uint16_t flag_supplementary = 0x800;

bool allowed_resource_uri(const std::string& uri) {
    return uri.find("://") == std::string::npos || uri.starts_with("https://");
}

bool allowed_resource_uri(const std::optional<std::string>& uri) {
    return !uri || allowed_resource_uri(*uri);
}

void validate_filters(const FilterSettings& filters) {
    if (filters.minimum_mapping_quality < 0 || filters.minimum_mapping_quality > 255
        || filters.minimum_base_quality < 0 || filters.minimum_base_quality > 255
        || filters.minimum_alternate_reads < 1 || filters.minimum_alternate_molecules < 1
        || filters.minimum_phasing_support < 1
        || !std::isfinite(filters.minimum_variant_allele_fraction)
        || filters.minimum_variant_allele_fraction < 0.0 || filters.minimum_variant_allele_fraction > 1.0
        || !std::isfinite(filters.maximum_phasing_conflict_fraction)
        || filters.maximum_phasing_conflict_fraction < 0.0 || filters.maximum_phasing_conflict_fraction > 1.0) {
        throw std::invalid_argument("Evidence filter settings are outside their valid ranges.");
    }
    if (filters.molecule_mode != MoleculeMode::raw_reads && filters.molecule_mode != MoleculeMode::auto_detect
        && filters.molecule_mode != MoleculeMode::selected_tag) {
        throw std::invalid_argument("Unknown molecule grouping mode.");
    }
    if (filters.molecule_mode == MoleculeMode::selected_tag
        && (filters.molecule_tag.size() != 2 || !std::isalpha(static_cast<unsigned char>(filters.molecule_tag[0]))
            || !std::isalnum(static_cast<unsigned char>(filters.molecule_tag[1])))) {
        throw std::invalid_argument("A selected BAM tag must contain exactly two valid tag characters.");
    }
}

void validate_variant_shape(const VariantQuery& query) {
    const auto valid_allele = [](const std::string& allele) {
        return !allele.empty() && allele.size() <= 10000 && std::all_of(allele.begin(), allele.end(), [](unsigned char base) {
            return base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'N';
        });
    };
    if (query.contig.empty() || query.position < 0 || !valid_allele(query.reference) || !valid_allele(query.alternate)
        || query.reference == query.alternate
        || (query.reference.size() != query.alternate.size() && !query.reference.starts_with(query.alternate)
            && !query.alternate.starts_with(query.reference))) {
        throw std::invalid_argument("Variant query is invalid or is not a supported left-anchored small variant.");
    }
}

std::string contig_key(std::string name) {
    std::transform(name.begin(), name.end(), name.begin(), [](const unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    if (name.starts_with("chr")) name.erase(0, 3);
    if (name == "m") name = "mt";
    return name;
}

std::optional<std::string> resolve_contig(
    const std::vector<igv::SequenceInfo>& sequences,
    const std::string& requested) {
    const auto exact = std::find_if(sequences.begin(), sequences.end(), [&](const auto& sequence) {
        return sequence.name == requested;
    });
    if (exact != sequences.end()) return exact->name;

    const auto key = contig_key(requested);
    std::optional<std::string> match;
    for (const auto& sequence : sequences) {
        if (contig_key(sequence.name) != key) continue;
        if (match) return std::nullopt;  // Do not guess if a header contains ambiguous aliases.
        match = sequence.name;
    }
    return match;
}

std::unique_ptr<igv::AlignmentReader> open_allowed_alignments(const igv::Resource& resource) {
    if (!allowed_resource_uri(resource.uri) || !allowed_resource_uri(resource.index_uri) || !allowed_resource_uri(resource.reference_uri)) {
        throw std::runtime_error("BAM Seek accepts local paths and HTTPS resources only; HTTP and other URL schemes are not permitted.");
    }
    return igv::open_alignments(resource);
}

struct ReadLayout {
    struct Insertion {
        std::string sequence;
        int minimum_quality{};
    };
    struct Deletion {
        igv::GenomicInterval interval;
        bool reference_skip = false;
    };
    std::map<std::int64_t, std::pair<char, int>> bases;
    std::map<std::int64_t, Insertion> insertions_after;
    std::vector<Deletion> deletions;
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
                    auto sequence = alignment.sequence.substr(read_position, length);
                    std::transform(sequence.begin(), sequence.end(), sequence.begin(), [](unsigned char c) {
                        return static_cast<char>(std::toupper(c));
                    });
                    int minimum_quality = std::numeric_limits<int>::max();
                    for (std::size_t i = 0; i < length; ++i) {
                        const int quality = read_position + i < alignment.qualities.size()
                            ? static_cast<int>(static_cast<unsigned char>(alignment.qualities[read_position + i])) - 33
                            : 0;
                        minimum_quality = std::min(minimum_quality, quality);
                    }
                    layout.insertions_after[reference_position - 1] = {std::move(sequence), minimum_quality};
                }
                read_position += length;
                break;
            case 'D': case 'N':
                layout.deletions.push_back({{alignment.interval.contig, reference_position,
                                             reference_position + static_cast<std::int64_t>(length)}, operation.operation == 'N'});
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
        return deletion.interval.start <= position && position < deletion.interval.end;
    });
}

std::optional<int> covered_quality(const ReadLayout& layout, const std::int64_t start, const std::size_t length) {
    int minimum_quality = std::numeric_limits<int>::max();
    for (std::size_t i = 0; i < length; ++i) {
        const auto position = start + static_cast<std::int64_t>(i);
        const auto base = layout.bases.find(position);
        if (base == layout.bases.end() || deleted_at(layout, position)) return std::nullopt;
        minimum_quality = std::min(minimum_quality, base->second.second);
    }
    return minimum_quality == std::numeric_limits<int>::max() ? std::optional<int>{} : minimum_quality;
}

std::optional<int> insertion_quality(const ReadLayout& layout, const std::int64_t first_anchor, const std::int64_t last_anchor) {
    int minimum_quality = std::numeric_limits<int>::max();
    for (auto insertion = layout.insertions_after.lower_bound(first_anchor);
         insertion != layout.insertions_after.end() && insertion->first <= last_anchor; ++insertion) {
        minimum_quality = std::min(minimum_quality, insertion->second.minimum_quality);
    }
    return minimum_quality == std::numeric_limits<int>::max() ? std::optional<int>{} : minimum_quality;
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
    const auto anchor_quality = covered_quality(layout, query.position, 1).value_or(0);

    if (query.reference.size() == query.alternate.size()) {
        const auto quality = covered_quality(layout, query.position, query.reference.size());
        if (!quality) return {Allele::other, anchor_quality};
        if (const auto inserted_quality = insertion_quality(layout, query.position, anchor)) {
            return {Allele::other, std::min(*quality, *inserted_quality)};
        }
        if (alt_match) return {Allele::alternate, *alt_match};
        if (ref_match) return {Allele::reference, *ref_match};
        return {Allele::other, *quality};
    }
    if (query.alternate.size() > query.reference.size()) {
        const auto insertion = layout.insertions_after.find(anchor);
        const auto expected = query.alternate.substr(query.reference.size());
        const bool deletion_at_boundary = std::any_of(layout.deletions.begin(), layout.deletions.end(), [anchor](const auto& deletion) {
            return deletion.interval.start <= anchor + 1 && anchor + 1 < deletion.interval.end;
        });
        if (ref_match && insertion != layout.insertions_after.end() && insertion->second.sequence == expected) {
            return {Allele::alternate, std::min(*ref_match, insertion->second.minimum_quality)};
        }
        if (ref_match && insertion == layout.insertions_after.end() && !deletion_at_boundary) return {Allele::reference, *ref_match};
        const auto other_quality = insertion == layout.insertions_after.end()
            ? ref_match.value_or(anchor_quality)
            : std::min(ref_match.value_or(anchor_quality), insertion->second.minimum_quality);
        return {Allele::other, other_quality};
    }
    const auto deletion_start = query.position + static_cast<std::int64_t>(query.alternate.size());
    const auto deletion_end = query.position + static_cast<std::int64_t>(query.reference.size());
    const bool has_expected_deletion = std::any_of(layout.deletions.begin(), layout.deletions.end(), [&](const auto& deletion) {
        return !deletion.reference_skip && deletion.interval.start == deletion_start && deletion.interval.end == deletion_end;
    });
    const auto prefix = sequence_matches(layout, query.position, query.alternate);
    const auto inserted_quality = insertion_quality(layout, query.position, deletion_end - 1);
    if (prefix && has_expected_deletion && !inserted_quality) {
        return {Allele::alternate, *prefix};
    }
    if (ref_match && !inserted_quality) return {Allele::reference, *ref_match};
    auto other_quality = ref_match.value_or(prefix.value_or(anchor_quality));
    if (inserted_quality) other_quality = std::min(other_quality, *inserted_quality);
    return {Allele::other, other_quality};
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

std::string choose_tag(const std::vector<const igv::Alignment*>& alignments, const FilterSettings& filters) {
    if (filters.molecule_mode == MoleculeMode::raw_reads) return {};
    if (filters.molecule_mode == MoleculeMode::selected_tag) return filters.molecule_tag;
    for (const std::string candidate : {"MI", "RX", "UB"}) {
        const auto tagged = std::count_if(alignments.begin(), alignments.end(), [&](const auto* alignment) {
            return tag_value(*alignment, candidate).has_value();
        });
        if (tagged > 0 && static_cast<double>(tagged) / static_cast<double>(alignments.size()) >= 0.9) return candidate;
    }
    return {};
}

std::string molecule_key(const igv::Alignment& alignment, const std::string& tag) {
    const auto value = tag_value(alignment, tag);
    if (!value) return {};
    if (tag == "MI") return *value;
    if (tag != "RX" && tag != "UB") return *value;
    auto fragment_start = alignment.interval.start;
    if (alignment.mate && alignment.mate->contig == alignment.interval.contig) {
        fragment_start = std::min(fragment_start, alignment.mate->start);
    }
    std::ostringstream key;
    if (tag == "UB") key << tag_value(alignment, "CB").value_or("no-CB") << '|';
    key << *value << '|' << alignment.interval.contig << ':' << fragment_start << '|' << std::abs(alignment.template_length);
    return key.str();
}

std::string fragment_key(const igv::Alignment& alignment) {
    if (!alignment.name.empty()) return "pair|" + alignment.name;
    std::ostringstream key;
    key << "alignment|" << alignment.interval.contig << ':' << alignment.interval.start << '-'
        << alignment.interval.end << '|' << alignment.flags;
    return key.str();
}

std::string alignment_record_key(const igv::Alignment& alignment) {
    std::ostringstream key;
    key << alignment.name << '|' << alignment.interval.contig << ':' << alignment.interval.start << '-'
        << alignment.interval.end << '|' << alignment.flags << '|';
    for (const auto& operation : alignment.cigar) key << operation.length << operation.operation;
    key << '|' << alignment.sequence;
    return key.str();
}

void add_count(EvidenceCounts& counts, const Allele allele, const bool reverse) {
    switch (allele) {
        case Allele::reference:
            ++counts.reference_reads;
            if (reverse) ++counts.reference_reverse_reads; else ++counts.reference_forward_reads;
            break;
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

Allele consensus_allele(const std::vector<Allele>& calls) {
    const auto alt = static_cast<int>(std::count(calls.begin(), calls.end(), Allele::alternate));
    const auto ref = static_cast<int>(std::count(calls.begin(), calls.end(), Allele::reference));
    const auto other = static_cast<int>(std::count(calls.begin(), calls.end(), Allele::other));
    return alt > ref && alt > other ? Allele::alternate
        : ref > alt && ref > other ? Allele::reference
        : Allele::other;
}

VariantEvidence evaluate_variant(const igv::AlignmentReader& reader, const VariantQuery& variant, const FilterSettings& filters) {
    VariantEvidence evidence;
    evidence.query = variant;
    const auto alignments = reader.get(variant.query_window());
    struct CalledRead {
        const igv::Alignment* alignment{};
        AlleleCall call{};
    };
    std::vector<CalledRead> called_reads;
    std::vector<const igv::Alignment*> callable_alignments;
    for (const auto& alignment : alignments) {
        ++evidence.overlapping_alignments;
        if (!include_alignment(alignment, filters)) {
            ++evidence.filtered_alignments;
            continue;
        }
        const auto called = call_allele(alignment, variant);
        if (called.allele == Allele::no_call) {
            ++evidence.uncallable_alignments;
            continue;
        }
        if (called.minimum_base_quality < filters.minimum_base_quality) {
            ++evidence.low_base_quality_alignments;
            continue;
        }
        called_reads.push_back({&alignment, called});
        callable_alignments.push_back(&alignment);
    }
    const auto selected_tag = choose_tag(callable_alignments, filters);
    std::unordered_map<std::string, std::vector<Allele>> molecule_calls;
    bool used_fragment_fallback = false;
    for (const auto& called_read : called_reads) {
        const auto& alignment = *called_read.alignment;
        const auto called = called_read.call;
        const bool reverse = (alignment.flags & flag_reverse) != 0;
        auto molecule_id = selected_tag.empty() ? std::string{} : molecule_key(alignment, selected_tag);
        if (molecule_id.empty()) {
            if (!selected_tag.empty()) ++evidence.reads_missing_molecule_tag;
            molecule_id = fragment_key(alignment);
            used_fragment_fallback = true;
        }
        std::ostringstream summary;
        summary << allele_name(called.allele) << " mapQ=" << static_cast<int>(alignment.mapping_quality)
                << " baseQ=" << called.minimum_base_quality << " CIGAR=";
        for (const auto& operation : alignment.cigar) summary << operation.length << operation.operation;
        evidence.reads.push_back({alignment.name, called.allele, reverse, alignment.mapping_quality,
                                  called.minimum_base_quality, molecule_id, summary.str()});
        add_count(evidence.counts, called.allele, reverse);
        molecule_calls[molecule_id].push_back(called.allele);
    }
    for (const auto& [identifier, calls] : molecule_calls) {
        (void)identifier;
        add_molecule_count(evidence.counts, consensus_allele(calls));
    }
    evidence.molecule_counts_available = !molecule_calls.empty();
    evidence.molecule_tag_used = selected_tag.empty() ? "paired fragments by read name" : selected_tag;
    if (!selected_tag.empty() && used_fragment_fallback) evidence.molecule_tag_used += " + pair fallback";
    const bool molecule_threshold_met = evidence.molecule_counts_available
        && evidence.counts.alternate_molecules >= filters.minimum_alternate_molecules;
    evidence.passes_thresholds = evidence.counts.alternate_reads >= filters.minimum_alternate_reads
        && evidence.counts.allele_fraction() >= filters.minimum_variant_allele_fraction
        && molecule_threshold_met;
    return evidence;
}

double log_choose(const int n, const int k) {
    if (k < 0 || k > n) return -std::numeric_limits<double>::infinity();
    return std::lgamma(static_cast<double>(n + 1)) - std::lgamma(static_cast<double>(k + 1))
        - std::lgamma(static_cast<double>(n - k + 1));
}

double hypergeometric_probability(const int x, const int row_one, const int column_one, const int total) {
    return std::exp(log_choose(column_one, x) + log_choose(total - column_one, row_one - x) - log_choose(total, row_one));
}

}  // namespace

std::string evidence_cache_key(const VariantQuery& query, const FilterSettings& filters) {
    std::ostringstream key;
    key << contig_key(query.contig) << ':' << query.position << ':'
        << query.reference.size() << ':' << query.reference << ':'
        << query.alternate.size() << ':' << query.alternate << '|'
        << filters.minimum_mapping_quality << ':' << filters.minimum_base_quality << ':'
        << filters.include_duplicates << ':' << filters.include_secondary << ':' << filters.include_supplementary << ':'
        << std::setprecision(std::numeric_limits<double>::max_digits10) << filters.minimum_variant_allele_fraction << ':'
        << filters.minimum_alternate_reads << ':' << filters.minimum_alternate_molecules << ':'
        << static_cast<int>(filters.molecule_mode) << ':' << filters.molecule_tag;
    return key.str();
}

double EvidenceCounts::allele_fraction() const noexcept {
    const auto denominator = informative_read_depth();
    return denominator == 0 ? 0.0 : static_cast<double>(alternate_reads) / static_cast<double>(denominator);
}

double EvidenceCounts::molecule_allele_fraction() const noexcept {
    const auto denominator = molecule_depth();
    return denominator == 0 ? 0.0 : static_cast<double>(alternate_molecules) / static_cast<double>(denominator);
}

std::optional<double> EvidenceCounts::strand_bias_p_value() const noexcept {
    if (alternate_reads == 0 || reference_reads == 0) return std::nullopt;
    const int forward = reference_forward_reads + alternate_forward_reads;
    const int reverse = reference_reverse_reads + alternate_reverse_reads;
    const int total = forward + reverse;
    if (forward == 0 || reverse == 0 || total == 0) return std::nullopt;
    const int minimum = std::max(0, alternate_reads - reverse);
    const int maximum = std::min(alternate_reads, forward);
    const auto observed = hypergeometric_probability(alternate_forward_reads, alternate_reads, forward, total);
    double p_value = 0.0;
    for (int x = minimum; x <= maximum; ++x) {
        const auto probability = hypergeometric_probability(x, alternate_reads, forward, total);
        if (probability <= observed + 1e-12) p_value += probability;
    }
    return std::min(1.0, p_value);
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

std::string phase_name(const PhaseClassification classification) {
    switch (classification) {
        case PhaseClassification::cis: return "cis";
        case PhaseClassification::trans: return "trans";
        case PhaseClassification::indeterminate: return "indeterminate";
        case PhaseClassification::too_far_apart: return "too far apart to phase";
    }
    return "indeterminate";
}

EvidenceEngine::EvidenceEngine(igv::Resource resource) : reader_(open_allowed_alignments(resource)) {
    if (!reader_->indexed()) throw std::runtime_error("An indexed BAM, CRAM, or SAM resource is required (.bai, .csi, or .crai).");
    if (resource.reference_uri) reference_ = igv::open_fasta({.uri = *resource.reference_uri});
}

bool EvidenceEngine::indexed() const noexcept { return reader_->indexed(); }

PhaseEvidence EvidenceEngine::phase_pair(
    const VariantQuery& first, const VariantQuery& second, const FilterSettings& filters) const {
    validate_filters(filters);
    validate_variant_shape(first);
    validate_variant_shape(second);

    PhaseEvidence evidence;
    evidence.first = first;
    evidence.second = second;
    evidence.gene = first.gene.empty() ? second.gene : first.gene;

    const auto first_contig = resolve_contig(reader_->sequences(), first.contig);
    const auto second_contig = resolve_contig(reader_->sequences(), second.contig);
    if (!first_contig) throw std::invalid_argument("Phasing contig is not present in the alignment header: " + first.contig);
    if (!second_contig) throw std::invalid_argument("Phasing contig is not present in the alignment header: " + second.contig);
    evidence.first.contig = *first_contig;
    evidence.second.contig = *second_contig;

    const auto validate_bounds_and_reference = [this](const VariantQuery& query, const std::string& requested_contig) {
        const auto sequence = std::find_if(reader_->sequences().begin(), reader_->sequences().end(), [&](const auto& item) {
            return item.name == query.contig;
        });
        if (sequence == reader_->sequences().end()
            || query.position + static_cast<std::int64_t>(query.reference.size()) > sequence->length) {
            throw std::invalid_argument("Phasing variant is outside the contig bounds.");
        }
        if (!reference_) return;
        const auto reference_contig = resolve_contig(reference_->sequences(), requested_contig);
        if (!reference_contig) {
            throw std::invalid_argument("Phasing contig is not present in the configured reference: " + requested_contig);
        }
        auto observed = reference_->get({*reference_contig, query.position,
            query.position + static_cast<std::int64_t>(query.reference.size())});
        std::transform(observed.begin(), observed.end(), observed.begin(), [](const unsigned char base) {
            return static_cast<char>(std::toupper(base));
        });
        if (observed != query.reference) {
            throw std::invalid_argument("Phasing query REF allele does not match the configured reference (observed "
                + observed + ").");
        }
    };
    validate_bounds_and_reference(evidence.first, first.contig);
    validate_bounds_and_reference(evidence.second, second.contig);

    if (*first_contig != *second_contig) {
        evidence.classification = PhaseClassification::too_far_apart;
        evidence.reason = "The variants are on different contigs, so no sequenced molecule can link them.";
        return evidence;
    }
    evidence.genomic_distance = evidence.first.position > evidence.second.position
        ? evidence.first.position - evidence.second.position
        : evidence.second.position - evidence.first.position;

    const auto first_alignments = reader_->get(evidence.first.query_window(0));
    const auto second_alignments = reader_->get(evidence.second.query_window(0));
    struct CalledAlignment {
        const igv::Alignment* alignment{};
        int locus{};
        Allele allele = Allele::no_call;
    };
    std::vector<CalledAlignment> called_alignments;
    std::vector<const igv::Alignment*> callable_alignments;
    const auto collect_calls = [&](const auto& alignments, const VariantQuery& variant, const int locus) {
        for (const auto& alignment : alignments) {
            if (!include_alignment(alignment, filters)) continue;
            const auto call = call_allele(alignment, variant);
            if (call.allele == Allele::no_call || call.minimum_base_quality < filters.minimum_base_quality) continue;
            called_alignments.push_back({&alignment, locus, call.allele});
            callable_alignments.push_back(&alignment);
        }
    };
    collect_calls(first_alignments, evidence.first, 0);
    collect_calls(second_alignments, evidence.second, 1);

    const auto selected_tag = choose_tag(callable_alignments, filters);
    struct PhaseCall {
        Allele allele = Allele::no_call;
        std::string alignment_record;
        std::uint16_t flags{};
    };
    struct MoleculeCalls {
        std::vector<PhaseCall> first;
        std::vector<PhaseCall> second;
        bool grouped_by_tag = false;
    };
    std::unordered_map<std::string, MoleculeCalls> molecule_calls;
    bool used_fragment_fallback = false;
    for (const auto& called : called_alignments) {
        auto molecule_id = selected_tag.empty() ? std::string{} : molecule_key(*called.alignment, selected_tag);
        const bool grouped_by_tag = !molecule_id.empty();
        if (molecule_id.empty()) {
            if (!selected_tag.empty()) ++evidence.reads_missing_molecule_tag;
            molecule_id = fragment_key(*called.alignment);
            used_fragment_fallback = true;
        }
        auto& calls = molecule_calls[molecule_id];
        calls.grouped_by_tag = calls.grouped_by_tag || grouped_by_tag;
        (called.locus == 0 ? calls.first : calls.second).push_back(
            {called.allele, alignment_record_key(*called.alignment), called.alignment->flags});
    }
    evidence.molecule_tag_used = selected_tag.empty() ? "paired fragments by read name" : selected_tag;
    if (!selected_tag.empty() && used_fragment_fallback) evidence.molecule_tag_used += " + pair fallback";

    const auto is_informative_allele = [](const Allele allele) {
        return allele == Allele::reference || allele == Allele::alternate;
    };
    for (const auto& [identifier, calls] : molecule_calls) {
        (void)identifier;
        const bool same_read_observes_both = std::any_of(calls.first.begin(), calls.first.end(), [&](const auto& first_call) {
            return std::any_of(calls.second.begin(), calls.second.end(), [&](const auto& second_call) {
                return first_call.alignment_record == second_call.alignment_record;
            });
        });
        const bool valid_read_pair_observes_both = std::any_of(calls.first.begin(), calls.first.end(), [&](const auto& first_call) {
            return std::any_of(calls.second.begin(), calls.second.end(), [&](const auto& second_call) {
                const bool valid_flags = (first_call.flags & (flag_paired | flag_proper_pair))
                        == (flag_paired | flag_proper_pair)
                    && (second_call.flags & (flag_paired | flag_proper_pair)) == (flag_paired | flag_proper_pair);
                const bool opposite_mates = ((first_call.flags & flag_first_mate) != 0
                        && (second_call.flags & flag_second_mate) != 0)
                    || ((first_call.flags & flag_second_mate) != 0
                        && (second_call.flags & flag_first_mate) != 0);
                return valid_flags && opposite_mates;
            });
        });
        if (!calls.grouped_by_tag && !same_read_observes_both && !valid_read_pair_observes_both) {
            ++evidence.counts.noninformative_molecules;
            continue;
        }
        const auto phase_consensus = [](const std::vector<PhaseCall>& observations) {
            std::vector<Allele> alleles;
            alleles.reserve(observations.size());
            for (const auto& observation : observations) alleles.push_back(observation.allele);
            return consensus_allele(alleles);
        };
        const auto first_call = phase_consensus(calls.first);
        const auto second_call = phase_consensus(calls.second);
        if (!is_informative_allele(first_call) || !is_informative_allele(second_call)) {
            ++evidence.counts.noninformative_molecules;
            continue;
        }
        if (first_call == Allele::alternate && second_call == Allele::alternate) {
            ++evidence.counts.alternate_alternate;
        } else if (first_call == Allele::alternate) {
            ++evidence.counts.alternate_reference;
        } else if (second_call == Allele::alternate) {
            ++evidence.counts.reference_alternate;
        } else {
            ++evidence.counts.reference_reference;
        }
    }

    const auto informative = evidence.counts.informative_molecules();
    evidence.direct_phasing_possible = informative > 0;
    if (!evidence.direct_phasing_possible) {
        evidence.classification = PhaseClassification::too_far_apart;
        evidence.reason = "No high-quality molecule directly observed both variant positions.";
        return evidence;
    }

    const int cis_support = evidence.counts.cis_supporting_molecules();
    const int trans_support = evidence.counts.trans_supporting_molecules();
    const int alternate_bearing = cis_support + trans_support;
    if (alternate_bearing == 0) {
        evidence.classification = PhaseClassification::indeterminate;
        evidence.reason = "Molecules span both positions, but none contains either alternate allele.";
        return evidence;
    }
    const auto cis_conflict_fraction = static_cast<double>(trans_support) / static_cast<double>(alternate_bearing);
    const auto trans_conflict_fraction = static_cast<double>(cis_support) / static_cast<double>(alternate_bearing);
    const bool cis_call = cis_support >= filters.minimum_phasing_support
        && cis_conflict_fraction <= filters.maximum_phasing_conflict_fraction;
    const bool trans_call = evidence.counts.alternate_reference >= filters.minimum_phasing_support
        && evidence.counts.reference_alternate >= filters.minimum_phasing_support
        && trans_conflict_fraction <= filters.maximum_phasing_conflict_fraction;

    if (cis_call && !trans_call) {
        evidence.classification = PhaseClassification::cis;
        evidence.discordant_molecules = trans_support;
        evidence.discordant_fraction = cis_conflict_fraction;
        evidence.reason = "ALT/ALT molecules meet the configured support and conflict thresholds.";
    } else if (trans_call && !cis_call) {
        evidence.classification = PhaseClassification::trans;
        evidence.discordant_molecules = cis_support;
        evidence.discordant_fraction = trans_conflict_fraction;
        evidence.reason = "Reciprocal ALT/REF and REF/ALT molecules meet the configured support and conflict thresholds.";
    } else {
        evidence.classification = PhaseClassification::indeterminate;
        evidence.discordant_molecules = std::min(cis_support, trans_support);
        evidence.discordant_fraction = static_cast<double>(evidence.discordant_molecules)
            / static_cast<double>(alternate_bearing);
        if (cis_support > 0 && trans_support > 0) {
            evidence.reason = "Conflicting cis- and trans-supporting molecules exceed the configured conflict threshold.";
        } else {
            evidence.reason = "Direct linkage was observed, but the configured molecule-support threshold was not met.";
        }
    }
    return evidence;
}

PileupData EvidenceEngine::pileup(const VariantQuery& query, const FilterSettings& filters, const std::int64_t padding) const {
    validate_filters(filters);
    validate_variant_shape(query);
    if (padding < 0 || padding > 1000) throw std::invalid_argument("Pileup padding must be between 0 and 1,000 bases.");
    PileupData data;
    data.query = query;
    const auto alignment_contig = resolve_contig(reader_->sequences(), query.contig);
    if (!alignment_contig) throw std::invalid_argument("Pileup contig is not present in the alignment header: " + query.contig);
    data.query.contig = *alignment_contig;
    data.interval = data.query.query_window(padding);
    data.minimum_base_quality = filters.minimum_base_quality;
    const auto sequence = std::find_if(reader_->sequences().begin(), reader_->sequences().end(), [&](const auto& item) {
        return item.name == *alignment_contig;
    });
    if (query.position + static_cast<std::int64_t>(query.reference.size()) > sequence->length) {
        throw std::invalid_argument("Pileup variant is outside the contig bounds.");
    }
    data.interval.end = std::min(data.interval.end, sequence->length);
    if (reference_) {
        const auto reference_contig = resolve_contig(reference_->sequences(), query.contig);
        if (!reference_contig) throw std::invalid_argument("Pileup contig is not present in the configured reference: " + query.contig);
        auto reference_interval = data.interval;
        reference_interval.contig = *reference_contig;
        data.reference_bases = reference_->get(reference_interval);
        data.has_reference = true;
    }
    constexpr std::size_t maximum_display_alignments = 5000;
    for (const auto& alignment : reader_->get(data.interval)) {
        if (!include_alignment(alignment, filters)) continue;
        ++data.total_alignments;
        if (data.alignments.size() < maximum_display_alignments) data.alignments.push_back(alignment);
    }
    data.truncated = data.total_alignments > data.alignments.size();
    return data;
}

BatchEvidence EvidenceEngine::evaluate(const std::vector<Query>& queries, const FilterSettings& filters) const {
    validate_filters(filters);
    BatchEvidence batch;
    for (const auto& query : queries) {
        try {
            auto resolved_query = query;
            auto* region = std::get_if<RegionQuery>(&resolved_query);
            if (!region) validate_variant_shape(std::get<VariantQuery>(resolved_query));
            else region->interval.validate();
            const auto requested_contig = region ? region->interval.contig : std::get<VariantQuery>(resolved_query).contig;
            const auto alignment_contig = resolve_contig(reader_->sequences(), requested_contig);
            if (!alignment_contig) throw std::runtime_error("contig is not present in the alignment header: " + requested_contig);
            if (region) region->interval.contig = *alignment_contig;
            else std::get<VariantQuery>(resolved_query).contig = *alignment_contig;
            const auto interval = region ? region->interval : std::get<VariantQuery>(resolved_query).query_window(0);
            const auto sequence = std::find_if(reader_->sequences().begin(), reader_->sequences().end(), [&](const auto& item) {
                return item.name == interval.contig;
            });
            if (interval.start < 0 || interval.end > sequence->length) throw std::runtime_error("query is outside the contig bounds");

            if (region) {
                if (!reference_) {
                    batch.results.emplace_back(RegionEvidence{*region, "Region scanning requires an indexed hg19 FASTA reference.", {}});
                    continue;
                }
                if (region->interval.end - region->interval.start > 100000) {
                    batch.results.emplace_back(RegionEvidence{*region, "Region is larger than 100 kb. Split it into smaller windows for targeted scanning.", {}});
                    continue;
                }
                const auto reference_contig = resolve_contig(reference_->sequences(), requested_contig);
                if (!reference_contig) throw std::runtime_error("contig is not present in the configured reference: " + requested_contig);
                auto reference_interval = region->interval;
                reference_interval.contig = *reference_contig;
                const auto reference_bases = reference_->get(reference_interval);
                if (reference_bases.size() != static_cast<std::size_t>(region->interval.end - region->interval.start)) {
                    throw std::runtime_error("reference did not return the complete requested region");
                }
                std::map<std::int64_t, std::map<char, int>> alternate_observations;
                const auto alignments = reader_->get(region->interval);
                for (const auto& alignment : alignments) {
                    if (!include_alignment(alignment, filters)) continue;
                    const auto layout = layout_for(alignment);
                    for (const auto& [position, base] : layout.bases) {
                        if (position < region->interval.start || position >= region->interval.end || base.second < filters.minimum_base_quality) continue;
                        const auto offset = static_cast<std::size_t>(position - region->interval.start);
                        const auto reference_base = static_cast<char>(std::toupper(static_cast<unsigned char>(reference_bases[offset])));
                        if (base.first != reference_base && reference_base != 'N' && base.first != 'N') ++alternate_observations[position][base.first];
                    }
                }
                RegionEvidence region_evidence{*region, "No candidates met the selected thresholds.", {}};
                constexpr std::size_t maximum_region_candidates = 500;
                bool truncated = false;
                for (const auto& [position, observed] : alternate_observations) {
                    const auto offset = static_cast<std::size_t>(position - region->interval.start);
                    const auto reference_base = static_cast<char>(std::toupper(static_cast<unsigned char>(reference_bases[offset])));
                    for (const auto& [alternate_base, preliminary_count] : observed) {
                        if (preliminary_count < filters.minimum_alternate_reads) continue;
                        auto candidate = evaluate_variant(*reader_, {region->source_text, region->interval.contig, position,
                                                                      std::string(1, reference_base), std::string(1, alternate_base), {}, {}, {}, {}, {}}, filters);
                        if (!candidate.passes_thresholds) continue;
                        if (region_evidence.candidates.size() == maximum_region_candidates) {
                            truncated = true;
                            break;
                        }
                        region_evidence.candidates.push_back(std::move(candidate));
                    }
                    if (truncated) break;
                }
                if (!region_evidence.candidates.empty()) {
                    region_evidence.note = std::to_string(region_evidence.candidates.size()) + " SNV candidate(s) met thresholds";
                    region_evidence.note += truncated ? " (result limit reached)." : ".";
                }
                batch.results.emplace_back(std::move(region_evidence));
                continue;
            }

            const auto& variant = std::get<VariantQuery>(resolved_query);
            if (reference_) {
                const auto reference_contig = resolve_contig(reference_->sequences(), requested_contig);
                if (!reference_contig) throw std::runtime_error("contig is not present in the configured reference: " + requested_contig);
                auto observed_reference = reference_->get({*reference_contig, variant.position,
                    variant.position + static_cast<std::int64_t>(variant.reference.size())});
                std::transform(observed_reference.begin(), observed_reference.end(), observed_reference.begin(), [](unsigned char base) {
                    return static_cast<char>(std::toupper(base));
                });
                if (observed_reference != variant.reference) {
                    throw std::runtime_error("query REF allele does not match the configured reference (observed " + observed_reference + ")");
                }
            }
            batch.results.emplace_back(evaluate_variant(*reader_, variant, filters));
        } catch (const std::exception& error) {
            const auto source = std::visit([](const auto& item) { return item.source_text; }, query);
            batch.errors.push_back(source + ": " + error.what());
        }
    }
    return batch;
}

}  // namespace bamseek
