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

    [[nodiscard]] igv::GenomicInterval query_window(std::int64_t padding = 8) const;
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

[[nodiscard]] ParsedQueries parse_queries(const std::string& text);

}  // namespace bamseek
