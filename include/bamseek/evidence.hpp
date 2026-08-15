#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <bamseek/query.hpp>
#include <igv/readers.hpp>
#include <igv/resource.hpp>

namespace bamseek {

// raw_reads groups alignments into paired fragments by read name; it does not use a UMI.
enum class MoleculeMode { raw_reads, auto_detect, selected_tag };
enum class Allele { reference, alternate, other, no_call };
enum class PhaseClassification { cis, trans, indeterminate, too_far_apart };

struct FilterSettings {
    int minimum_mapping_quality = 20;
    int minimum_base_quality = 20;
    bool include_duplicates = false;
    bool include_secondary = false;
    bool include_supplementary = false;
    double minimum_variant_allele_fraction = 0.0;
    int minimum_alternate_reads = 1;
    int minimum_alternate_molecules = 1;
    int minimum_phasing_support = 2;
    double maximum_phasing_conflict_fraction = 0.1;
    MoleculeMode molecule_mode = MoleculeMode::raw_reads;
    std::string molecule_tag;
};

[[nodiscard]] std::string evidence_cache_key(const VariantQuery& query, const FilterSettings& filters);

struct ReadEvidence {
    std::string read_name;
    Allele allele = Allele::no_call;
    bool reverse_strand = false;
    int mapping_quality{};
    int minimum_base_quality{};
    std::string molecule_id;
    std::string summary;
};

struct EvidenceCounts {
    int reference_reads{};
    int alternate_reads{};
    int other_reads{};
    int reference_forward_reads{};
    int reference_reverse_reads{};
    int alternate_forward_reads{};
    int alternate_reverse_reads{};
    int reference_molecules{};
    int alternate_molecules{};
    int other_molecules{};

    [[nodiscard]] int depth() const noexcept { return reference_reads + alternate_reads + other_reads; }
    [[nodiscard]] int informative_read_depth() const noexcept { return reference_reads + alternate_reads; }
    [[nodiscard]] int molecule_depth() const noexcept { return reference_molecules + alternate_molecules; }
    [[nodiscard]] int total_molecule_depth() const noexcept {
        return reference_molecules + alternate_molecules + other_molecules;
    }
    [[nodiscard]] double allele_fraction() const noexcept;
    [[nodiscard]] double molecule_allele_fraction() const noexcept;
    [[nodiscard]] std::optional<double> strand_bias_p_value() const noexcept;
};

struct VariantEvidence {
    VariantQuery query;
    EvidenceCounts counts;
    int overlapping_alignments{};
    int filtered_alignments{};
    int uncallable_alignments{};
    int low_base_quality_alignments{};
    std::string molecule_tag_used;
    bool molecule_counts_available = false;
    int reads_missing_molecule_tag{};
    bool passes_thresholds = false;
    std::vector<ReadEvidence> reads;
};

struct RegionEvidence {
    RegionQuery query;
    std::string note;
    std::vector<VariantEvidence> candidates;
};

struct PhaseCounts {
    int alternate_alternate{};  // AB
    int alternate_reference{};  // Ab
    int reference_alternate{};  // aB
    int reference_reference{};  // ab
    int noninformative_molecules{};

    [[nodiscard]] int informative_molecules() const noexcept {
        return alternate_alternate + alternate_reference + reference_alternate + reference_reference;
    }
    [[nodiscard]] int cis_supporting_molecules() const noexcept { return alternate_alternate; }
    [[nodiscard]] int trans_supporting_molecules() const noexcept {
        return alternate_reference + reference_alternate;
    }
};

struct PhaseEvidence {
    VariantQuery first;
    VariantQuery second;
    std::string gene;
    std::int64_t genomic_distance = -1;
    bool direct_phasing_possible = false;
    PhaseCounts counts;
    PhaseClassification classification = PhaseClassification::indeterminate;
    int discordant_molecules{};
    double discordant_fraction{};
    std::string molecule_tag_used;
    int reads_missing_molecule_tag{};
    std::string reason;
};

using QueryEvidence = std::variant<VariantEvidence, RegionEvidence>;

struct BatchEvidence {
    std::vector<QueryEvidence> results;
    std::vector<std::string> errors;
};

struct PileupData {
    VariantQuery query;
    igv::GenomicInterval interval;
    std::string reference_bases;
    bool has_reference = false;
    int minimum_base_quality{};
    std::vector<igv::Alignment> alignments;
    std::size_t total_alignments{};
    bool truncated = false;
};

class EvidenceEngine {
public:
    explicit EvidenceEngine(igv::Resource resource);
    [[nodiscard]] bool indexed() const noexcept;
    [[nodiscard]] BatchEvidence evaluate(const std::vector<Query>& queries, const FilterSettings& filters) const;
    [[nodiscard]] PhaseEvidence phase_pair(
        const VariantQuery& first, const VariantQuery& second, const FilterSettings& filters) const;
    [[nodiscard]] PileupData pileup(const VariantQuery& query, const FilterSettings& filters, std::int64_t padding = 80) const;

private:
    std::unique_ptr<igv::AlignmentReader> reader_;
    std::unique_ptr<igv::ReferenceReader> reference_;
};

[[nodiscard]] std::string allele_name(Allele allele);
[[nodiscard]] std::string phase_name(PhaseClassification classification);

}  // namespace bamseek
