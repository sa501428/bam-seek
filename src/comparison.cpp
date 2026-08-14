#include <bamseek/comparison.hpp>

#include <algorithm>
#include <cctype>
#include <map>
#include <sstream>
#include <tuple>

namespace bamseek {
namespace {

std::string uppercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char character) {
        return static_cast<char>(std::toupper(character));
    });
    return value;
}

bool same_nonempty(const std::string& left, const std::string& right) {
    return !left.empty() && !right.empty() && uppercase(left) == uppercase(right);
}

void merge_metadata(VariantQuery& target, const VariantQuery& source) {
    if (target.gene.empty()) target.gene = source.gene;
    if (target.transcript.empty()) target.transcript = source.transcript;
    if (target.coding_change.empty()) target.coding_change = source.coding_change;
    if (target.protein_change.empty()) target.protein_change = source.protein_change;
    if (target.variant_type.empty()) target.variant_type = source.variant_type;
}

void add_classified(std::vector<ClassifiedVariant>& classified, const VariantQuery& query, const VariantOrigin origin) {
    const auto existing = std::find_if(classified.begin(), classified.end(), [&](const auto& candidate) {
        return same_annotated_variant(candidate.query, query);
    });
    if (existing == classified.end()) {
        classified.push_back({query, origin});
        return;
    }
    merge_metadata(existing->query, query);
    if (existing->origin != origin) existing->origin = VariantOrigin::both;
}

std::string natural_list(const std::vector<std::string>& values) {
    if (values.empty()) return {};
    if (values.size() == 1) return values.front();
    if (values.size() == 2) return values[0] + " and " + values[1];
    std::string joined;
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index > 0) joined += index + 1 == values.size() ? ", and " : ", ";
        joined += values[index];
    }
    return joined;
}

std::string annotated_name(const VariantQuery& query) {
    std::string label = query.gene;
    if (!query.coding_change.empty()) {
        if (!label.empty()) label += ' ';
        label += query.coding_change;
    } else if (!query.protein_change.empty()) {
        if (!label.empty()) label += ' ';
        label += query.protein_change;
    }
    if (!label.empty()) return label;
    return query.contig + ':' + std::to_string(query.position + 1) + ' '
        + query.reference + '>' + query.alternate;
}

std::string evidence_item(const VariantEvidence& evidence, const bool force_zero_reads, const bool include_molecules) {
    const auto& counts = evidence.counts;
    std::string label = annotated_name(evidence.query);
    std::vector<std::string> details;
    if (!evidence.query.coding_change.empty() && !evidence.query.protein_change.empty()) {
        details.push_back(evidence.query.protein_change);
    }
    details.push_back(std::to_string(force_zero_reads ? 0 : counts.alternate_reads) + '/'
        + std::to_string(counts.alternate_reads + counts.reference_reads) + " consensus reads");
    if (include_molecules) {
        const auto molecules = evidence.molecule_counts_available
            ? std::to_string(counts.alternate_molecules) : std::string("N/A");
        details.push_back(molecules + (molecules == "1" ? " molecule" : " molecules"));
    }
    std::string detail_text;
    for (const auto& detail : details) {
        if (!detail_text.empty()) detail_text += ", ";
        detail_text += detail;
    }
    return label + " (" + detail_text + ')';
}

std::string detected_paragraph(const std::string& bam, const std::vector<const VariantEvidence*>& variants) {
    std::vector<std::string> labels;
    labels.reserve(variants.size());
    for (const auto* variant : variants) labels.push_back(evidence_item(*variant, false, true));
    const bool plural = variants.size() > 1;
    return "Targeted manual review of sequence reads from previous sample " + bam
        + " from this patient shows that the " + natural_list(labels) + (plural
            ? " variants reported in this sample were present at that time, but at VAFs below the validated threshold for calling variants and, therefore, were not reported at that time. Accordingly, their presence now represents an increase in relative abundance from the previous sample, but not de novo new mutations for this patient."
            : " variant reported in this sample was present at that time, but at a VAF below the validated threshold for calling a variant and, therefore, was not reported at that time. Accordingly, its presence now represents an increase in relative abundance from the previous sample, but not a de novo new mutation for this patient.");
}

std::string absent_paragraph(const std::string& bam, const std::vector<const VariantEvidence*>& variants) {
    std::vector<std::string> labels;
    labels.reserve(variants.size());
    for (const auto* variant : variants) labels.push_back(evidence_item(*variant, true, false));
    const bool plural = variants.size() > 1;
    return "Targeted manual review of sequence reads from previous sample " + bam
        + " from this patient shows that the " + natural_list(labels) + (plural
            ? " variants reported in this sample were not present and, therefore, their presence in this sample represents new mutations for this patient."
            : " variant reported in this sample was not present and, therefore, its presence in this sample represents a new mutation for this patient.");
}

std::string lost_paragraph(const std::vector<const VariantEvidence*>& variants) {
    std::vector<std::string> labels;
    labels.reserve(variants.size());
    for (const auto* variant : variants) labels.push_back(evidence_item(*variant, false, true));
    return "A detailed investigation in IGV does NOT identify the previously reported " + natural_list(labels)
        + (variants.size() > 1 ? " mutations" : " mutation") + " above the level of noise of the assay.";
}

}  // namespace

bool same_annotated_variant(const VariantQuery& left, const VariantQuery& right) {
    const bool same_gene = same_nonempty(left.gene, right.gene);
    if (same_gene && (same_nonempty(left.coding_change, right.coding_change)
                      || same_nonempty(left.protein_change, right.protein_change))) return true;
    if (!left.gene.empty() || !right.gene.empty()) return false;
    return uppercase(left.contig) == uppercase(right.contig) && left.position == right.position
        && uppercase(left.reference) == uppercase(right.reference)
        && uppercase(left.alternate) == uppercase(right.alternate);
}

std::vector<ClassifiedVariant> classify_variants(
    const std::vector<VariantQuery>& current, const std::vector<VariantQuery>& historical) {
    std::vector<ClassifiedVariant> classified;
    classified.reserve(current.size() + historical.size());
    for (const auto& query : current) add_classified(classified, query, VariantOrigin::current);
    for (const auto& query : historical) add_classified(classified, query, VariantOrigin::historical);
    return classified;
}

std::string comparison_narrative(const std::vector<ComparativeEvidence>& evidence) {
    struct HistoricalBamEvidence {
        std::vector<const VariantEvidence*> detected;
        std::vector<const VariantEvidence*> absent;
    };
    std::map<std::string, HistoricalBamEvidence> historical_bams;
    std::vector<const VariantEvidence*> lost_from_current;
    for (const auto& item : evidence) {
        if (item.variant_origin == VariantOrigin::both) continue;
        if (!item.bam_is_current && item.variant_origin == VariantOrigin::current) {
            auto& group = historical_bams[item.bam_name];
            if (item.evidence.molecule_counts_available && item.evidence.counts.alternate_molecules > 0) {
                group.detected.push_back(&item.evidence);
            } else {
                group.absent.push_back(&item.evidence);
            }
        } else if (item.bam_is_current && item.variant_origin == VariantOrigin::historical) {
            lost_from_current.push_back(&item.evidence);
        }
    }

    std::vector<std::string> paragraphs;
    for (const auto& [bam, group] : historical_bams) {
        if (!group.detected.empty()) paragraphs.push_back(detected_paragraph(bam, group.detected));
        if (!group.absent.empty()) paragraphs.push_back(absent_paragraph(bam, group.absent));
    }
    if (!lost_from_current.empty()) paragraphs.push_back(lost_paragraph(lost_from_current));

    std::ostringstream joined;
    for (std::size_t index = 0; index < paragraphs.size(); ++index) {
        if (index > 0) joined << "\n\n";
        joined << paragraphs[index];
    }
    return joined.str();
}

}  // namespace bamseek
