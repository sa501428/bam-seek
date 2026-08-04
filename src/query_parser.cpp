#include <bamseek/query.hpp>

#include <algorithm>
#include <cctype>
#include <charconv>
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
    return reference.size() == alternate.size() || reference.starts_with(alternate) || alternate.starts_with(reference);
}

}  // namespace

igv::GenomicInterval VariantQuery::query_window(const std::int64_t padding) const {
    return {contig, std::max<std::int64_t>(0, position - padding), position + static_cast<std::int64_t>(reference.size()) + padding};
}

ParsedQueries parse_queries(const std::string& text) {
    ParsedQueries result;
    std::istringstream lines(text);
    std::string raw_line;
    int line_number = 0;
    while (std::getline(lines, raw_line)) {
        ++line_number;
        const auto line = trim(raw_line);
        if (line.empty() || line.starts_with('#')) continue;
        std::istringstream fields(line);
        std::vector<std::string> columns;
        for (std::string column; fields >> column;) columns.push_back(column);
        if (columns.size() == 4 || columns.size() >= 5) {
            const bool vcf_like = columns.size() >= 5;
            const auto& contig = columns[0];
            const auto& position_text = columns[1];
            const auto& reference_text = columns[vcf_like ? 3 : 2];
            const auto& alternate_text = columns[vcf_like ? 4 : 3];
            std::int64_t one_based_position{};
            const auto reference = uppercase(reference_text);
            const auto alternate = uppercase(alternate_text);
            if (!parse_positive(position_text, one_based_position) || !supported_alleles(reference, alternate)) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": invalid variant (small indels must be left-anchored)");
            } else {
                result.queries.emplace_back(VariantQuery{line, contig, one_based_position - 1, reference, alternate});
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
            if (allele_arrow == std::string::npos || !parse_positive(position_text, one_based_position)) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": expected chr:position REF>ALT");
                continue;
            }
            const auto reference = uppercase(trim(allele_text.substr(0, allele_arrow)));
            const auto alternate = uppercase(trim(allele_text.substr(allele_arrow + 1)));
            if (!supported_alleles(reference, alternate)) {
                result.errors.push_back("Line " + std::to_string(line_number) + ": invalid variant (small indels must be left-anchored)");
                continue;
            }
            result.queries.emplace_back(VariantQuery{line, contig, one_based_position - 1, reference, alternate});
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

}  // namespace bamseek
