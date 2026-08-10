#pragma once

#include <cstdint>
#include <string>
#include <variant>
#include <vector>

#include <igv/interval.hpp>

namespace bamseek {

struct VariantQuery {
    std::string source_text;
    std::string contig;
    std::int64_t position{};  // Zero-based, left-most VCF base.
    std::string reference;
    std::string alternate;
    std::string gene;
    std::string transcript;
    std::string coding_change;
    std::string protein_change;

    [[nodiscard]] igv::GenomicInterval query_window(std::int64_t padding = 8) const;
};

struct ClinicalVariantMapping {
    std::string gene;
    std::string transcript;
    std::string coding_change;
    std::string protein_change;
    std::string contig;
    std::int64_t position{};  // Zero-based genomic coordinate.
    std::string reference;
    std::string alternate;
};

struct LoadedClinicalMappings {
    std::vector<ClinicalVariantMapping> mappings;
    std::vector<std::string> errors;
};

struct RegionQuery {
    std::string source_text;
    igv::GenomicInterval interval;
};

using Query = std::variant<VariantQuery, RegionQuery>;

struct ParsedQueries {
    std::vector<Query> queries;
    std::vector<std::string> errors;
};

[[nodiscard]] ParsedQueries parse_queries(
    const std::string& text,
    const std::vector<ClinicalVariantMapping>& clinical_mappings = {});
[[nodiscard]] LoadedClinicalMappings load_clinical_mappings(const std::string& local_path);

}  // namespace bamseek
