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

enum class MoleculeMode { raw_reads, auto_detect, selected_tag };
enum class Allele { reference, alternate, other, no_call };

struct FilterSettings {
    int minimum_mapping_quality = 20;
    int minimum_base_quality = 20;
    bool include_duplicates = false;
    bool include_secondary = false;
    bool include_supplementary = false;
    double minimum_variant_allele_fraction = 0.0005;
    int minimum_alternate_reads = 1;
    int minimum_alternate_molecules = 1;
    MoleculeMode molecule_mode = MoleculeMode::auto_detect;
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
    int alternate_forward_reads{};
    int alternate_reverse_reads{};
    int reference_molecules{};
    int alternate_molecules{};
    int other_molecules{};

    [[nodiscard]] int depth() const noexcept { return reference_reads + alternate_reads + other_reads; }
    [[nodiscard]] double allele_fraction() const noexcept;
};

struct VariantEvidence {
    VariantQuery query;
    EvidenceCounts counts;
    std::string molecule_tag_used;
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

class EvidenceEngine {
public:
    explicit EvidenceEngine(igv::Resource resource);
    [[nodiscard]] bool indexed() const noexcept;
    [[nodiscard]] BatchEvidence evaluate(const std::vector<Query>& queries, const FilterSettings& filters) const;

private:
    std::unique_ptr<igv::AlignmentReader> reader_;
    std::unique_ptr<igv::ReferenceReader> reference_;
};

[[nodiscard]] std::string allele_name(Allele allele);

}  // namespace bamseek
