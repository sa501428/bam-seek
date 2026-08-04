#include <bamseek/evidence.hpp>

#include <cstdlib>
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
    if (snv.counts.reference_reads != 2 || snv.counts.alternate_reads != 3 || snv.counts.other_reads != 1
        || snv.counts.reference_forward_reads != 1 || snv.counts.reference_reverse_reads != 1
        || snv.counts.alternate_forward_reads != 2 || snv.counts.alternate_reverse_reads != 1
        || !snv.molecule_counts_available || snv.counts.alternate_molecules != 2 || snv.counts.reference_molecules != 2
        || snv.reads_missing_molecule_tag != 1 || !snv.counts.strand_bias_p_value()) {
        std::cerr << "SNV, strand, or molecule evidence regression\n";
        return EXIT_FAILURE;
    }

    auto auto_filters = molecule_filters;
    auto_filters.molecule_mode = bamseek::MoleculeMode::auto_detect;
    const auto auto_batch = evidence_engine.evaluate(
        {bamseek::VariantQuery{"chr1:15 A>T", "chr1", 14, "A", "T"}}, auto_filters);
    const auto& auto_snv = std::get<bamseek::VariantEvidence>(auto_batch.results.front());
    if (auto_snv.molecule_counts_available || auto_snv.counts.alternate_molecules != 0) {
        std::cerr << "Incomplete auto-detected tags must not produce molecule counts\n";
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
    return EXIT_SUCCESS;
}
