#include <bamseek/evidence.hpp>

#include <cstdlib>
#include <cmath>
#include <iostream>

int main() {
    bamseek::EvidenceEngine engine({.uri = BAM_SEEK_TEST_BAM, .reference_uri = BAM_SEEK_TEST_REFERENCE});
    bamseek::FilterSettings filters;
    filters.minimum_mapping_quality = 0;
    filters.minimum_base_quality = 0;
    filters.molecule_mode = bamseek::MoleculeMode::raw_reads;

    bamseek::VariantQuery variant{"chr1:3 G>A", "chr1", 2, "G", "A"};
    const auto batch = engine.evaluate({variant}, filters);
    if (!batch.errors.empty() || batch.results.size() != 1) return EXIT_FAILURE;
    const auto& evidence = std::get<bamseek::VariantEvidence>(batch.results.front());
    if (evidence.counts.reference_reads != 1 || evidence.counts.alternate_reads != 0 || evidence.counts.depth() != 1) {
        std::cerr << "Unexpected targeted evidence counts: ref=" << evidence.counts.reference_reads
                  << " alt=" << evidence.counts.alternate_reads << " other=" << evidence.counts.other_reads
                  << " depth=" << evidence.counts.depth() << " reads=" << evidence.reads.size() << '\n';
        return EXIT_FAILURE;
    }

    const auto pileup = engine.pileup(variant, filters, 2);
    if (!pileup.has_reference || pileup.alignments.size() != 1 || pileup.reference_bases != "ACGTA") {
        std::cerr << "Unexpected pileup data\n";
        return EXIT_FAILURE;
    }

    bamseek::RegionQuery region{"chr1:1-12", {"chr1", 0, 12}};
    const auto region_batch = engine.evaluate({region}, filters);
    if (!region_batch.errors.empty() || region_batch.results.size() != 1) return EXIT_FAILURE;
    const auto& region_evidence = std::get<bamseek::RegionEvidence>(region_batch.results.front());
    if (!region_evidence.candidates.empty()) {
        std::cerr << "Reference-matching test reads should not yield candidates\n";
        return EXIT_FAILURE;
    }
    bool rejected_insecure_resource = false;
    try {
        bamseek::EvidenceEngine insecure({.uri = "http://example.com/sample.bam"});
        (void)insecure;
    } catch (const std::runtime_error&) {
        rejected_insecure_resource = true;
    }
    if (!rejected_insecure_resource) {
        std::cerr << "Insecure remote resource should have been rejected before access\n";
        return EXIT_FAILURE;
    }

    bamseek::EvidenceEngine evidence_engine({.uri = BAM_SEEK_EVIDENCE_BAM, .reference_uri = BAM_SEEK_EVIDENCE_REFERENCE});
    bamseek::FilterSettings molecule_filters;
    molecule_filters.minimum_mapping_quality = 20;
    molecule_filters.minimum_base_quality = 20;
    molecule_filters.molecule_mode = bamseek::MoleculeMode::selected_tag;
    molecule_filters.molecule_tag = "MI";
    const auto snv_batch = evidence_engine.evaluate({bamseek::VariantQuery{"chr1:15 A>T", "chr1", 14, "A", "T"}}, molecule_filters);
    if (!snv_batch.errors.empty() || snv_batch.results.size() != 1) return EXIT_FAILURE;
    const auto& snv = std::get<bamseek::VariantEvidence>(snv_batch.results.front());
    if (snv.counts.reference_reads != 2 || snv.counts.alternate_reads != 3 || snv.counts.other_reads != 2
        || snv.counts.reference_forward_reads != 1 || snv.counts.reference_reverse_reads != 1
        || snv.counts.alternate_forward_reads != 2 || snv.counts.alternate_reverse_reads != 1
        || std::abs(snv.counts.allele_fraction() - 0.6) > 1e-12
        || !snv.molecule_counts_available || snv.counts.alternate_molecules != 3 || snv.counts.reference_molecules != 2
        || snv.counts.other_molecules != 2
        || std::abs(snv.counts.molecule_allele_fraction() - 0.6) > 1e-12
        || snv.reads_missing_molecule_tag != 1 || !snv.counts.strand_bias_p_value()) {
        std::cerr << "SNV, strand, or molecule evidence regression\n";
        return EXIT_FAILURE;
    }

    auto auto_filters = molecule_filters;
    auto_filters.molecule_mode = bamseek::MoleculeMode::auto_detect;
    const auto auto_batch = evidence_engine.evaluate(
        {bamseek::VariantQuery{"chr1:15 A>T", "chr1", 14, "A", "T"}}, auto_filters);
    const auto& auto_snv = std::get<bamseek::VariantEvidence>(auto_batch.results.front());
    if (!auto_snv.molecule_counts_available || auto_snv.counts.alternate_molecules != 3
        || auto_snv.molecule_tag_used != "read pairs/fragments") {
        std::cerr << "Auto grouping must fall back to read pairs/fragments\n";
        return EXIT_FAILURE;
    }

    const auto paired_batch = evidence_engine.evaluate(
        {bamseek::VariantQuery{"chr1:75 A>T", "chr1", 74, "A", "T"}}, auto_filters);
    if (!paired_batch.errors.empty() || paired_batch.results.size() != 1) return EXIT_FAILURE;
    const auto& paired = std::get<bamseek::VariantEvidence>(paired_batch.results.front());
    if (paired.counts.reference_reads != 2 || paired.counts.alternate_reads != 3
        || paired.counts.informative_read_depth() != 5 || std::abs(paired.counts.allele_fraction() - 0.6) > 1e-12
        || paired.counts.reference_molecules != 1 || paired.counts.alternate_molecules != 2
        || paired.counts.molecule_depth() != 3
        || std::abs(paired.counts.molecule_allele_fraction() - (2.0 / 3.0)) > 1e-12) {
        std::cerr << "Read-pair molecule consensus regression\n";
        return EXIT_FAILURE;
    }

    bamseek::EvidenceCounts biallelic_denominator;
    biallelic_denominator.reference_reads = 6;
    biallelic_denominator.alternate_reads = 4;
    biallelic_denominator.other_reads = 90;
    biallelic_denominator.reference_molecules = 3;
    biallelic_denominator.alternate_molecules = 2;
    biallelic_denominator.other_molecules = 45;
    if (std::abs(biallelic_denominator.allele_fraction() - 0.4) > 1e-12
        || std::abs(biallelic_denominator.molecule_allele_fraction() - 0.4) > 1e-12) {
        std::cerr << "OTHER/N observations must be excluded from both VAF denominators\n";
        return EXIT_FAILURE;
    }

    auto raw_filters = molecule_filters;
    raw_filters.molecule_mode = bamseek::MoleculeMode::raw_reads;
    const auto insertion_batch = evidence_engine.evaluate(
        {bamseek::VariantQuery{"chr1:35 A>AT", "chr1", 34, "A", "AT"}}, raw_filters);
    const auto& insertion = std::get<bamseek::VariantEvidence>(insertion_batch.results.front());
    if (!insertion_batch.errors.empty() || insertion.counts.alternate_reads != 1 || insertion.counts.reference_reads != 1
        || insertion.counts.other_reads != 1 || insertion.counts.depth() != 3) {
        std::cerr << "Insertion quality or adjacent deletion classification regression\n";
        return EXIT_FAILURE;
    }

    const auto deletion_batch = evidence_engine.evaluate(
        {bamseek::VariantQuery{"chr1:55 AAA>A", "chr1", 54, "AAA", "A"}}, raw_filters);
    const auto& deletion = std::get<bamseek::VariantEvidence>(deletion_batch.results.front());
    if (!deletion_batch.errors.empty() || deletion.counts.alternate_reads != 1 || deletion.counts.reference_reads != 1
        || deletion.counts.other_reads != 1) {
        std::cerr << "Exact deletion classification regression\n";
        return EXIT_FAILURE;
    }

    const auto mismatch = evidence_engine.evaluate(
        {bamseek::VariantQuery{"chr1:15 C>T", "chr1", 14, "C", "T"}}, raw_filters);
    if (!mismatch.results.empty() || mismatch.errors.size() != 1) {
        std::cerr << "Reference mismatch should be reported as an error\n";
        return EXIT_FAILURE;
    }
    auto invalid_filters = raw_filters;
    invalid_filters.minimum_variant_allele_fraction = 1.1;
    bool rejected_invalid_filters = false;
    try {
        (void)evidence_engine.evaluate({bamseek::VariantQuery{"chr1:15 A>T", "chr1", 14, "A", "T"}}, invalid_filters);
    } catch (const std::invalid_argument&) {
        rejected_invalid_filters = true;
    }
    if (!rejected_invalid_filters) {
        std::cerr << "Invalid filter ranges must be rejected\n";
        return EXIT_FAILURE;
    }

    struct Jak2Fixture {
        const char* path;
        int depth;
        int alternate_reads;
        double vaf;
        const char* header_contig;
    };
    const Jak2Fixture jak2_fixtures[]{
        {BAM_SEEK_JAK2_NEGATIVE_BAM, 120, 0, 0.0, "chr9"},
        {BAM_SEEK_JAK2_25PCT_BAM, 150, 38, 38.0 / 150.0, "chr9"},
        {BAM_SEEK_JAK2_28PCT_BAM, 180, 50, 50.0 / 180.0, "9"},
    };
    bamseek::FilterSettings jak2_filters;
    jak2_filters.minimum_mapping_quality = 20;
    jak2_filters.minimum_base_quality = 20;
    jak2_filters.minimum_variant_allele_fraction = 0.05;
    jak2_filters.molecule_mode = bamseek::MoleculeMode::auto_detect;
    const bamseek::VariantQuery jak2_v617f{
        "JAK2 c.1849G>T p.V617F chr9:5073770 G>T", "chr9", 5073769, "G", "T",
        "JAK2", "NM_004972.4", "c.1849G>T", "p.V617F"};
    for (const auto& fixture : jak2_fixtures) {
        bamseek::EvidenceEngine jak2_engine({.uri = fixture.path});
        const auto result = jak2_engine.evaluate({jak2_v617f}, jak2_filters);
        if (!result.errors.empty() || result.results.size() != 1) {
            std::cerr << "JAK2 fixture could not be evaluated: " << fixture.path << '\n';
            return EXIT_FAILURE;
        }
        const auto& observed = std::get<bamseek::VariantEvidence>(result.results.front());
        if (observed.counts.depth() != fixture.depth || observed.counts.alternate_reads != fixture.alternate_reads
            || std::abs(observed.counts.allele_fraction() - fixture.vaf) > 1e-12
            || !observed.molecule_counts_available || observed.counts.alternate_molecules != fixture.alternate_reads
            || observed.passes_thresholds != (fixture.alternate_reads > 0)) {
            std::cerr << "Incorrect JAK2 evidence for " << fixture.path << ": depth=" << observed.counts.depth()
                      << " alt=" << observed.counts.alternate_reads << " VAF=" << observed.counts.allele_fraction()
                      << " alt molecules=" << observed.counts.alternate_molecules << '\n';
            return EXIT_FAILURE;
        }
        const auto jak2_pileup = jak2_engine.pileup(jak2_v617f, jak2_filters, 0);
        if (jak2_pileup.total_alignments != static_cast<std::size_t>(fixture.depth)
            || jak2_pileup.query.contig != fixture.header_contig) {
            std::cerr << "Incorrect aliased JAK2 pileup for " << fixture.path << '\n';
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
