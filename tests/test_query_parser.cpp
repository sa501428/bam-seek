#include <bamseek/query.hpp>

#include <cstdlib>
#include <iostream>

int main() {
    const auto parsed = bamseek::parse_queries("chr7:140453136 A>T\nchr1:1,000-1,100\nchr2:bad A>T\nchr3\t20\t.\tC\tG\n");
    if (parsed.queries.size() != 3 || parsed.errors.size() != 1) {
        std::cerr << "Unexpected parser result\n";
        return EXIT_FAILURE;
    }
    const auto& variant = std::get<bamseek::VariantQuery>(parsed.queries.front());
    if (variant.position != 140453135 || variant.reference != "A" || variant.alternate != "T") return EXIT_FAILURE;
    const auto& region = std::get<bamseek::RegionQuery>(parsed.queries[1]);
    if (region.interval.start != 999 || region.interval.end != 1100) return EXIT_FAILURE;
    const auto invalid = bamseek::parse_queries("chr1:10 A>A\nchr1:10 A><script>\n");
    if (!invalid.queries.empty() || invalid.errors.size() != 2) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
