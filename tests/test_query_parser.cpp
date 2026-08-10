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

    const auto inline_clinical = bamseek::parse_queries("BRAF c.1799T>A p.V600E chr7:140453136 A>T\n");
    if (!inline_clinical.errors.empty() || inline_clinical.queries.size() != 1) return EXIT_FAILURE;
    const auto& inline_variant = std::get<bamseek::VariantQuery>(inline_clinical.queries.front());
    if (inline_variant.gene != "BRAF" || inline_variant.coding_change != "c.1799T>A"
        || inline_variant.protein_change != "p.V600E" || inline_variant.contig != "chr7") return EXIT_FAILURE;

    const auto mappings = bamseek::load_clinical_mappings(BAM_SEEK_CLINICAL_MAPPING);
    if (!mappings.errors.empty() || mappings.mappings.size() != 4) return EXIT_FAILURE;
    const auto mapped = bamseek::parse_queries("BRAF NM_004333.6:c.1799T>A p.V600E\n", mappings.mappings);
    if (!mapped.errors.empty() || mapped.queries.size() != 1) return EXIT_FAILURE;
    const auto& mapped_variant = std::get<bamseek::VariantQuery>(mapped.queries.front());
    if (mapped_variant.position != 140453135 || mapped_variant.transcript != "NM_004333.6") return EXIT_FAILURE;
    const auto jak2 = bamseek::parse_queries("JAK2 NM_004972.4:c.1849G>T p.V617F\n", mappings.mappings);
    if (!jak2.errors.empty() || jak2.queries.size() != 1) return EXIT_FAILURE;
    const auto& jak2_variant = std::get<bamseek::VariantQuery>(jak2.queries.front());
    if (jak2_variant.position != 5073769 || jak2_variant.reference != "G" || jak2_variant.alternate != "T") return EXIT_FAILURE;

    const auto ambiguous = bamseek::parse_queries("TP53 c.743G>A p.R248Q\n", mappings.mappings);
    if (!ambiguous.queries.empty() || ambiguous.errors.size() != 1) return EXIT_FAILURE;
    const auto genomic_hgvs = bamseek::parse_queries("chr7:g.140453136A>T\n");
    if (!genomic_hgvs.errors.empty() || std::get<bamseek::VariantQuery>(genomic_hgvs.queries.front()).position != 140453135) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
