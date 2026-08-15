#include <bamseek/comparison.hpp>

#include <cstdlib>
#include <iostream>

namespace {

bamseek::VariantQuery variant(std::string gene, std::string coding, std::string protein,
                              std::string contig, const std::int64_t position) {
    return {gene + ' ' + coding, std::move(contig), position, "A", "T", std::move(gene), {},
        std::move(coding), std::move(protein), "missense_variant"};
}

bamseek::ComparativeEvidence evidence(const bamseek::VariantQuery& query, const std::string& bam,
                                      const bool current_bam, const bamseek::VariantOrigin origin,
                                      const int alt_reads, const int ref_reads, const int alt_molecules) {
    bamseek::EvidenceCounts counts;
    counts.alternate_reads = alt_reads;
    counts.reference_reads = ref_reads;
    counts.alternate_molecules = alt_molecules;
    counts.reference_molecules = ref_reads;
    bamseek::VariantEvidence observed;
    observed.query = query;
    observed.counts = counts;
    observed.molecule_counts_available = true;
    return {bam, current_bam, origin, std::move(observed)};
}

}  // namespace

int main() {
    const auto current_a = variant("JAK2", "c.1849G>T", "p.V617F", "chr9", 5073769);
    const auto current_a_duplicate = variant("jak2", "c.1849G>T", "p.V617F", "9", 5073769);
    const auto shared_current = variant("TP53", "c.743G>A", "p.R248Q", "chr17", 7674220);
    const auto shared_historical = variant("tp53", "c.999G>A", "p.R248Q", "17", 7674220);
    const auto historical_only = variant("TET2", "c.5162T>G", "p.L1721W", "chr4", 105235420);
    const auto classified = bamseek::classify_variants(
        {current_a, current_a_duplicate, shared_current}, {shared_historical, historical_only, historical_only});
    if (classified.size() != 3 || classified[0].origin != bamseek::VariantOrigin::current
        || classified[1].origin != bamseek::VariantOrigin::both
        || classified[2].origin != bamseek::VariantOrigin::historical) {
        std::cerr << "Variant origin deduplication failed\n";
        return EXIT_FAILURE;
    }

    const auto current_b = variant("NOTCH1", "c.7330C>T", "p.Q2444*", "chr9", 139390860);
    const std::vector<bamseek::ComparativeEvidence> results{
        evidence(current_a, "old-A.bam", false, bamseek::VariantOrigin::current, 2, 98, 1),
        evidence(current_b, "old-A.bam", false, bamseek::VariantOrigin::current, 0, 75, 0),
        evidence(shared_current, "old-A.bam", false, bamseek::VariantOrigin::both, 10, 90, 5),
        evidence(current_a, "old-B.bam", false, bamseek::VariantOrigin::current, 0, 80, 0),
        evidence(historical_only, "current.bam", true, bamseek::VariantOrigin::historical, 0, 120, 0),
    };
    const auto narrative = bamseek::comparison_narrative(results);
    if (narrative.find("previous sample old-A.bam") == std::string::npos
        || narrative.find("JAK2 c.1849G>T (p.V617F, 2/100 consensus reads, 1 molecule)") == std::string::npos
        || narrative.find("NOTCH1 c.7330C>T (p.Q2444*, 0/75 consensus reads)") == std::string::npos
        || narrative.find("previous sample old-B.bam") == std::string::npos
        || narrative.find("previously reported TET2 c.5162T>G (p.L1721W, 0/120 consensus reads, 0 molecules) mutation") == std::string::npos
        || narrative.find("TP53") != std::string::npos) {
        std::cerr << narrative << '\n';
        return EXIT_FAILURE;
    }

    const auto historical_second = variant("ASXL1", "c.1934dupG", "p.G646Wfs*12", "chr20", 32434637);
    const auto plural_narrative = bamseek::comparison_narrative({
        evidence(current_a, "old-C.bam", false, bamseek::VariantOrigin::current, 2, 98, 1),
        evidence(current_b, "old-C.bam", false, bamseek::VariantOrigin::current, 3, 97, 2),
        evidence(historical_only, "current.bam", true, bamseek::VariantOrigin::historical, 0, 120, 0),
        evidence(historical_second, "current.bam", true, bamseek::VariantOrigin::historical, 1, 119, 1),
    });
    if (plural_narrative.find("variants reported in this sample were present") == std::string::npos
        || plural_narrative.find("Accordingly, their presence now represents") == std::string::npos
        || plural_narrative.find("mutations above the level of noise") == std::string::npos) {
        std::cerr << "Plural comparison grammar failed\n" << plural_narrative << '\n';
        return EXIT_FAILURE;
    }

    bamseek::PhaseEvidence cis_phase;
    cis_phase.first = variant("JAK2", "c.1849G>T", "p.V617F", "chr9", 5073769);
    cis_phase.second = variant("JAK2", "c.1860C>A", "p.D620E", "chr9", 5073780);
    cis_phase.gene = "JAK2";
    cis_phase.genomic_distance = 11;
    cis_phase.direct_phasing_possible = true;
    cis_phase.counts.alternate_alternate = 3;
    cis_phase.counts.reference_reference = 7;
    cis_phase.classification = bamseek::PhaseClassification::cis;
    cis_phase.reason = "ALT/ALT molecules meet the configured support and conflict thresholds.";
    auto too_far_phase = cis_phase;
    too_far_phase.first = variant("TET2", "c.100A>T", {}, "chr4", 100);
    too_far_phase.second = variant("TET2", "c.900A>T", {}, "chr4", 900);
    too_far_phase.gene = "TET2";
    too_far_phase.genomic_distance = 800;
    too_far_phase.direct_phasing_possible = false;
    too_far_phase.counts = {};
    too_far_phase.classification = bamseek::PhaseClassification::too_far_apart;
    too_far_phase.reason = "No high-quality molecule directly observed both variant positions.";
    const auto phase_text = bamseek::phasing_narrative({
        {"current.bam", true, cis_phase}, {"old-A.bam", false, too_far_phase}});
    if (phase_text.find("JAK2 c.1849G>T and c.1860C>A are in cis") == std::string::npos
        || phase_text.find("AB=3, Ab=0, aB=0, ab=7") == std::string::npos
        || phase_text.find("TET2 c.100A>T and c.900A>T are too far apart to phase") == std::string::npos) {
        std::cerr << "Phasing narrative failed\n" << phase_text << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
