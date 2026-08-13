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

struct FilterSettings {
    int minimum_mapping_quality = 20;
    int minimum_base_quality = 20;
    bool include_duplicates = false;
    bool include_secondary = false;
    bool include_supplementary = false;
    double minimum_variant_allele_fraction = 0.0;
    int minimum_alternate_reads = 1;
    int minimum_alternate_molecules = 1;
    MoleculeMode molecule_mode = MoleculeMode::raw_reads;
    std::string molecule_tag;
};

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
    [[nodiscard]] PileupData pileup(const VariantQuery& query, const FilterSettings& filters, std::int64_t padding = 80) const;

private:
    std::unique_ptr<igv::AlignmentReader> reader_;
    std::unique_ptr<igv::ReferenceReader> reference_;
};

[[nodiscard]] std::string allele_name(Allele allele);

}  // namespace bamseek
