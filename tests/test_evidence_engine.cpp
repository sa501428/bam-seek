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

    bamseek::RegionQuery region{"chr1:1-12", {"chr1", 0, 12}};
    const auto region_batch = engine.evaluate({region}, filters);
    if (!region_batch.errors.empty() || region_batch.results.size() != 1) return EXIT_FAILURE;
    const auto& region_evidence = std::get<bamseek::RegionEvidence>(region_batch.results.front());
    if (!region_evidence.candidates.empty()) {
        std::cerr << "Reference-matching test reads should not yield candidates\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
