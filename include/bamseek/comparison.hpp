#pragma once

#include <bamseek/evidence.hpp>

#include <string>
#include <vector>

namespace bamseek {

enum class VariantOrigin { current, historical, both };

struct ClassifiedVariant {
    VariantQuery query;
    VariantOrigin origin = VariantOrigin::current;
};

struct ComparativeEvidence {
    std::string bam_name;
    bool bam_is_current = false;
    VariantOrigin variant_origin = VariantOrigin::current;
    VariantEvidence evidence;
};

struct ComparativePhaseEvidence {
    std::string bam_name;
    bool bam_is_current = false;
    PhaseEvidence evidence;
};

[[nodiscard]] bool same_annotated_variant(const VariantQuery& left, const VariantQuery& right);
[[nodiscard]] std::vector<ClassifiedVariant> classify_variants(
    const std::vector<VariantQuery>& current, const std::vector<VariantQuery>& historical);
[[nodiscard]] std::string comparison_narrative(const std::vector<ComparativeEvidence>& evidence);
[[nodiscard]] std::string phasing_narrative(const std::vector<ComparativePhaseEvidence>& evidence);

}  // namespace bamseek
