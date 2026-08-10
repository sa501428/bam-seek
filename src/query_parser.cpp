#include <bamseek/query.hpp>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <string_view>

namespace bamseek {
namespace {

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r");
    return value.substr(first, last - first + 1);
}

std::string uppercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return value;
}

bool parse_positive(std::string_view text, std::int64_t& result) {
    std::string stripped;
    stripped.reserve(text.size());
    for (const char character : text) if (character != ',') stripped += character;
    const auto [end, error] = std::from_chars(stripped.data(), stripped.data() + stripped.size(), result);
    return error == std::errc{} && end == stripped.data() + stripped.size() && result > 0;
}

bool valid_allele(const std::string& allele) {
    return !allele.empty() && std::all_of(allele.begin(), allele.end(), [](unsigned char c) {
        return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N';
    });
}

bool supported_alleles(const std::string& reference, const std::string& alternate) {
    if (!valid_allele(reference) || !valid_allele(alternate)) return false;
    if (reference == alternate || reference.size() > 10000 || alternate.size() > 10000) return false;
    return reference.size() == alternate.size() || reference.starts_with(alternate) || alternate.starts_with(reference);
}

bool safe_position(const std::int64_t one_based_position) {
    return one_based_position <= std::numeric_limits<std::int64_t>::max() - 20000;
}

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        const auto tab = line.find('\t', start);
        fields.push_back(trim(line.substr(start, tab == std::string::npos ? tab : tab - start)));
        if (tab == std::string::npos) break;
        start = tab + 1;
    }
    return fields;
}

std::string join_columns(const std::vector<std::string>& columns, const std::size_t start) {
    std::string joined;
    for (std::size_t i = start; i < columns.size(); ++i) {
        if (!joined.empty()) joined += ' ';
        joined += columns[i];
    }
    return joined;
}

struct ClinicalInput {
    std::string gene;
    std::string transcript;
    std::string coding_change;
    std::string protein_change;
    std::size_t consumed{};
};

std::optional<ClinicalInput> parse_clinical_prefix(const std::vector<std::string>& columns) {
    if (columns.size() < 2) return std::nullopt;
    ClinicalInput input;
    input.gene = columns[0];
    std::size_t index = 1;
    const auto coding_separator = columns[index].find(":c.");
    if (coding_separator != std::string::npos) {
        input.transcript = columns[index].substr(0, coding_separator);
        input.coding_change = columns[index].substr(coding_separator + 1);
        ++index;
    } else if (columns[index].starts_with("c.")) {
        input.coding_change = columns[index++];
    } else if (index + 1 < columns.size() && columns[index + 1].starts_with("c.")) {
        input.transcript = columns[index++];
        input.coding_change = columns[index++];
    } else {
        return std::nullopt;
    }
    if (index < columns.size() && columns[index].starts_with("p.")) input.protein_change = columns[index++];
    input.consumed = index;
    return input;
}

bool clinical_mapping_matches(const ClinicalVariantMapping& mapping, const ClinicalInput& input) {
    if (uppercase(mapping.gene) != uppercase(input.gene) || mapping.coding_change != input.coding_change) return false;
    if (!input.transcript.empty() && mapping.transcript != input.transcript) return false;
    return input.protein_change.empty() || mapping.protein_change.empty() || mapping.protein_change == input.protein_change;
}

}  // namespace

igv::GenomicInterval VariantQuery::query_window(const std::int64_t padding) const {
    return {contig, std::max<std::int64_t>(0, position - padding), position + static_cast<std::int64_t>(reference.size()) + padding};
}

ParsedQueries parse_queries(const std::string& text, const std::vector<ClinicalVariantMapping>& clinical_mappings) {
    ParsedQueries result;
    std::istringstream lines(text);
    std::string raw_line;
    int line_number = 0;
    while (std::getline(lines, raw_line)) {
        ++line_number;
        if (result.queries.size() + result.errors.size() >= 10000) {
            result.errors.push_back("Input is limited to 10,000 query lines per batch");
            break;
        }
        const auto line = trim(raw_line);
        if (line.empty() || line.starts_with('#')) continue;
        std::istringstream fields(line);
        std::vector<std::string> columns;
        for (std::string column; fields >> column;) columns.push_back(column);
        const auto genomic_hgvs = line.find(":g.");
        if (columns.size() == 1 && genomic_hgvs != std::string::npos) {
            const auto contig = line.substr(0, genomic_hgvs);
            const auto change = line.substr(genomic_hgvs + 3);
            const auto allele_start = change.find_first_of("ACGTNacgtn");
            const auto arrow = change.find('>', allele_start);
            std::int64_t one_based_position{};
            if (allele_start == std::string::npos || arrow == std::string::npos
                || !parse_positive(change.substr(0, allele_start), one_based_position) || !safe_position(one_based_position)) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": invalid genomic HGVS substitution");
                continue;
            }
            const auto reference = uppercase(change.substr(allele_start, arrow - allele_start));
            const auto alternate = uppercase(change.substr(arrow + 1));
            if (!supported_alleles(reference, alternate)) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": unsupported genomic HGVS allele");
                continue;
            }
            result.queries.emplace_back(VariantQuery{line, contig, one_based_position - 1, reference, alternate, {}, {}, {}, {}});
            continue;
        }
        if (const auto clinical = parse_clinical_prefix(columns)) {
            if (clinical->consumed < columns.size()) {
                const auto genomic_text = join_columns(columns, clinical->consumed);
                auto genomic = parse_queries(genomic_text);
                if (genomic.errors.empty() && genomic.queries.size() == 1 && std::holds_alternative<VariantQuery>(genomic.queries.front())) {
                    auto variant = std::get<VariantQuery>(std::move(genomic.queries.front()));
                    variant.source_text = line;
                    variant.gene = clinical->gene;
                    variant.transcript = clinical->transcript;
                    variant.coding_change = clinical->coding_change;
                    variant.protein_change = clinical->protein_change;
                    result.queries.emplace_back(std::move(variant));
                } else {
                    result.errors.push_back("Line " + std::to_string(line_number) + ": clinical notation must be followed by one genomic variant");
                }
                continue;
            }
            std::vector<const ClinicalVariantMapping*> matches;
            for (const auto& mapping : clinical_mappings) {
                if (clinical_mapping_matches(mapping, *clinical)) matches.push_back(&mapping);
            }
            if (matches.size() == 1) {
                const auto& mapping = *matches.front();
                result.queries.emplace_back(VariantQuery{line, mapping.contig, mapping.position, mapping.reference, mapping.alternate,
                    clinical->gene, clinical->transcript.empty() ? mapping.transcript : clinical->transcript,
                    clinical->coding_change, clinical->protein_change.empty() ? mapping.protein_change : clinical->protein_change});
            } else if (matches.empty()) {
                result.errors.push_back("Line " + std::to_string(line_number)
                    + ": no genomic mapping found for clinical notation; provide an inline genomic allele or a local mapping TSV");
            } else {
                result.errors.push_back("Line " + std::to_string(line_number)
                    + ": clinical notation maps to multiple loci/transcripts; include a transcript accession or genomic allele");
            }
            continue;
        }
        if (columns.size() == 4 || columns.size() >= 5) {
            const bool vcf_like = columns.size() >= 5;
            const auto& contig = columns[0];
            const auto& position_text = columns[1];
            const auto& reference_text = columns[vcf_like ? 3 : 2];
            const auto& alternate_text = columns[vcf_like ? 4 : 3];
            std::int64_t one_based_position{};
            const auto reference = uppercase(reference_text);
            const auto alternate = uppercase(alternate_text);
            if (!parse_positive(position_text, one_based_position) || !safe_position(one_based_position) || !supported_alleles(reference, alternate)) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": invalid variant (small indels must be left-anchored)");
            } else {
                result.queries.emplace_back(VariantQuery{line, contig, one_based_position - 1, reference, alternate, {}, {}, {}, {}});
            }
            continue;
        }
        const auto colon = line.find(':');
        if (colon == std::string::npos || colon == 0) {
            result.errors.push_back("Line " + std::to_string(line_number) + ": expected contig:position or contig:start-end");
            continue;
        }
        const auto contig = line.substr(0, colon);
        const auto payload = trim(line.substr(colon + 1));
        const auto arrow = payload.find('>');
        if (arrow != std::string::npos) {
            const auto space = payload.find_first_of(" \t");
            const auto allele_start = space == std::string::npos ? 0 : space + 1;
            const auto position_text = trim(payload.substr(0, allele_start));
            const auto allele_text = trim(payload.substr(allele_start));
            const auto allele_arrow = allele_text.find('>');
            std::int64_t one_based_position{};
            if (allele_arrow == std::string::npos || !parse_positive(position_text, one_based_position) || !safe_position(one_based_position)) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": expected chr:position REF>ALT");
                continue;
            }
            const auto reference = uppercase(trim(allele_text.substr(0, allele_arrow)));
            const auto alternate = uppercase(trim(allele_text.substr(allele_arrow + 1)));
            if (!supported_alleles(reference, alternate)) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": invalid variant (small indels must be left-anchored)");
                continue;
            }
            result.queries.emplace_back(VariantQuery{line, contig, one_based_position - 1, reference, alternate, {}, {}, {}, {}});
            continue;
        }
        const auto dash = payload.find('-');
        if (dash != std::string::npos) {
            std::int64_t start{};
            std::int64_t end{};
            if (!parse_positive(trim(payload.substr(0, dash)), start) || !parse_positive(trim(payload.substr(dash + 1)), end) || end < start) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": invalid one-based region");
                continue;
            }
            result.queries.emplace_back(RegionQuery{line, {contig, start - 1, end}});
            continue;
        }
        result.errors.push_back("Line " + std::to_string(line_number) + ": expected chr:position REF>ALT or chr:start-end");
    }
    return result;
}

LoadedClinicalMappings load_clinical_mappings(const std::string& local_path) {
    LoadedClinicalMappings loaded;
    if (local_path.empty()) return loaded;
    if (local_path.find("://") != std::string::npos) {
        loaded.errors.push_back("Clinical mapping tables must be local files.");
        return loaded;
    }
    std::ifstream input(local_path);
    if (!input) {
        loaded.errors.push_back("Could not open clinical mapping TSV: " + local_path);
        return loaded;
    }
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto cleaned = trim(line);
        if (cleaned.empty() || cleaned.starts_with('#')) continue;
        const auto fields = split_tabs(cleaned);
        if (fields.size() >= 2 && uppercase(fields[0]) == "GENE" && uppercase(fields[1]) == "TRANSCRIPT") continue;
        if (fields.size() != 8) {
            loaded.errors.push_back("Clinical mapping TSV line " + std::to_string(line_number) + ": expected 8 tab-separated columns");
            continue;
        }
        std::int64_t one_based_position{};
        const auto reference = uppercase(fields[6]);
        const auto alternate = uppercase(fields[7]);
        if (fields[0].empty() || fields[2].empty() || fields[4].empty()
            || !parse_positive(fields[5], one_based_position) || !safe_position(one_based_position)
            || !supported_alleles(reference, alternate)) {
            loaded.errors.push_back("Clinical mapping TSV line " + std::to_string(line_number) + ": invalid mapping values");
            continue;
        }
        if (loaded.mappings.size() == 100000) {
            loaded.errors.push_back("Clinical mapping TSV is limited to 100,000 entries");
            break;
        }
        loaded.mappings.push_back({fields[0], fields[1], fields[2], fields[3], fields[4], one_based_position - 1, reference, alternate});
    }
    return loaded;
}

}  // namespace bamseek
