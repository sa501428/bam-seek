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

    const auto igv_report = bamseek::parse_queries(
        "IGV\t9\t139390861\tNOTCH1\tstop_gained\tc.7330C>T\tp.Q2444*\tG\tA\t3.1\n"
        "IGV\t9\t139399347\tNOTCH1\tinframe_insertion\tc.4793_4795dupGCG\tp.R1598_V1599insG\tA\tACGC\t1.6\n"
        "10\t89717711\tPTEN\tframeshift_variant\tc.736_737insCACATTGTCTTATATAAAACGTCGGCTGTAGTCGT\tp.L247Hfs*5\tC\tCCACATTGTCTTATATAAAACGTCGGCTGTAGTCGT\t30.7\n"
        "IGV\n"
        "8\t117866581\tRAD21\tmissense_variant\tc.1064C>T\tp.P355L\tG\tA\t1.4\n"
        "IGV\t1\t100\tGENE1\tinframe_deletion\tc.10_11del\tp.K4del\tATC\tA\t12.3\n");
    if (!igv_report.errors.empty() || igv_report.queries.size() != 5) {
        std::cerr << "IGV report rows were not parsed\n";
        return EXIT_FAILURE;
    }
    const auto& notch_snv = std::get<bamseek::VariantQuery>(igv_report.queries[0]);
    const auto& notch_insertion = std::get<bamseek::VariantQuery>(igv_report.queries[1]);
    const auto& pten_insertion = std::get<bamseek::VariantQuery>(igv_report.queries[2]);
    const auto& deletion = std::get<bamseek::VariantQuery>(igv_report.queries[4]);
    if (notch_snv.contig != "9" || notch_snv.position != 139390860 || notch_snv.gene != "NOTCH1"
        || notch_snv.variant_type != "stop_gained" || notch_snv.coding_change != "c.7330C>T"
        || notch_snv.protein_change != "p.Q2444*" || notch_snv.reference != "G" || notch_snv.alternate != "A"
        || notch_insertion.reference != "A" || notch_insertion.alternate != "ACGC"
        || pten_insertion.reference != "C" || !pten_insertion.alternate.starts_with("C")
        || deletion.reference != "ATC" || deletion.alternate != "A" || deletion.variant_type != "inframe_deletion") {
        std::cerr << "IGV report metadata or left-anchored alleles were parsed incorrectly\n";
        return EXIT_FAILURE;
    }

    const auto transcript_report = bamseek::parse_queries(
        "IGV\n"
        "5\t176943930\tDDX41\tENST00000507955.1\tMISSENSE\tc.17C>T\tp.P6L\tG\tA\tGTTCGXGTTCC\t4.4\t19\t45\n"
        "IGV\n"
        "4\t106196838\tTET2\tENST00000380013.4\tFRAMESHIFT\tc.5171_5172insT\tp.P1725Sfs*4\tA\tAT\tTCCTTXTCCCA\t3.3\t600\n");
    if (!transcript_report.errors.empty() || transcript_report.queries.size() != 2) return EXIT_FAILURE;
    const auto& ddx41 = std::get<bamseek::VariantQuery>(transcript_report.queries[0]);
    const auto& tet2 = std::get<bamseek::VariantQuery>(transcript_report.queries[1]);
    if (ddx41.transcript != "ENST00000507955.1" || ddx41.variant_type != "MISSENSE"
        || ddx41.reference != "G" || ddx41.alternate != "A"
        || tet2.transcript != "ENST00000380013.4" || tet2.variant_type != "FRAMESHIFT"
        || tet2.reference != "A" || tet2.alternate != "AT") return EXIT_FAILURE;

    const auto unanchored_report = bamseek::parse_queries(
        "IGV\t1\t100\tGENE\tFRAMESHIFT\tc.1_2insT\tp.X1fs\t-\tT\t50\n");
    if (!unanchored_report.queries.empty() || unanchored_report.errors.size() != 1) {
        std::cerr << "Unanchored IGV indels must be rejected\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
