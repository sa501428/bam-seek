#include <bamseek/pileup_view.hpp>

#include <QPainter>
#include <QPaintEvent>
#include <QFontDatabase>

#include <algorithm>
#include <cctype>
#include <map>
#include <numeric>

namespace bamseek {
namespace {

constexpr int left_margin = 82;
constexpr int top_margin = 38;
constexpr int cell_width = 12;
constexpr int row_height = 17;

QColor strand_color(const bool reverse) {
    return reverse ? QColor(246, 183, 177) : QColor(171, 211, 238);
}

bool reverse_strand(const igv::Alignment& alignment) { return (alignment.flags & 0x10U) != 0; }

char reference_at(const PileupData& data, const std::int64_t position) {
    const auto offset = position - data.interval.start;
    if (!data.has_reference || offset < 0 || static_cast<std::size_t>(offset) >= data.reference_bases.size()) return 'N';
    return static_cast<char>(std::toupper(static_cast<unsigned char>(data.reference_bases[static_cast<std::size_t>(offset)])));
}

}  // namespace

PileupView::PileupView(QWidget* parent) : QWidget(parent) {
    setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    setAutoFillBackground(true);
}

void PileupView::set_data(PileupData data) {
    data_ = std::move(data);
    rebuild_layout();
    updateGeometry();
    update();
}

void PileupView::set_group_pairs(const bool enabled) {
    if (group_pairs_ == enabled) return;
    group_pairs_ = enabled;
    rebuild_layout();
    updateGeometry();
    update();
}

QSize PileupView::sizeHint() const {
    const auto columns = std::max<std::int64_t>(1, data_.interval.end - data_.interval.start);
    const auto row_count = rows_.empty() ? 1 : *std::max_element(rows_.begin(), rows_.end()) + 1;
    return {left_margin + static_cast<int>(columns) * cell_width + 20, top_margin + row_count * row_height + 30};
}

void PileupView::rebuild_layout() {
    order_.resize(data_.alignments.size());
    std::iota(order_.begin(), order_.end(), std::size_t{0});
    std::sort(order_.begin(), order_.end(), [&](const std::size_t first, const std::size_t second) {
        const auto& lhs = data_.alignments[first];
        const auto& rhs = data_.alignments[second];
        if (group_pairs_ && lhs.name != rhs.name) return lhs.name < rhs.name;
        if (lhs.interval.start != rhs.interval.start) return lhs.interval.start < rhs.interval.start;
        return lhs.name < rhs.name;
    });
    rows_.assign(data_.alignments.size(), 0);
    std::vector<std::int64_t> row_ends;
    for (const auto index : order_) {
        const auto& alignment = data_.alignments[index];
        int row = 0;
        while (row < static_cast<int>(row_ends.size()) && row_ends[static_cast<std::size_t>(row)] >= alignment.interval.start) ++row;
        if (row == static_cast<int>(row_ends.size())) row_ends.push_back(alignment.interval.end);
        else row_ends[static_cast<std::size_t>(row)] = alignment.interval.end;
        rows_[index] = row;
    }
    setMinimumSize(sizeHint());
}

void PileupView::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    painter.fillRect(rect(), Qt::white);
    painter.setFont(font());
    if (data_.interval.contig.empty()) {
        painter.drawText(12, 24, "Select a targeted variant and choose View pileup.");
        return;
    }

    painter.setPen(QColor(70, 70, 70));
    painter.drawText(8, 17, QString::fromStdString(data_.interval.contig));
    const auto length = data_.interval.end - data_.interval.start;
    for (std::int64_t offset = 0; offset < length; ++offset) {
        const auto position = data_.interval.start + offset;
        const int x = left_margin + static_cast<int>(offset) * cell_width;
        if ((position + 1) % 10 == 0 || position == data_.query.position) {
            painter.setPen(position == data_.query.position ? QColor(158, 42, 43) : QColor(90, 90, 90));
            painter.drawText(x - 9, 17, QString::number(position + 1));
            painter.drawLine(x + cell_width / 2, top_margin - 11, x + cell_width / 2, height());
        }
        painter.setPen(QColor(50, 50, 50));
        painter.drawText(x + 2, top_margin - 2, QString(reference_at(data_, position)));
    }

    if (group_pairs_) {
        std::map<std::string, std::vector<std::size_t>> pairs;
        for (const auto index : order_) pairs[data_.alignments[index].name].push_back(index);
        painter.setPen(QPen(QColor(120, 120, 120), 1, Qt::DashLine));
        for (const auto& [name, reads] : pairs) {
            (void)name;
            if (reads.size() < 2) continue;
            const auto first = reads[0];
            const auto second = reads[1];
            const int first_x = left_margin + static_cast<int>(data_.alignments[first].interval.end - data_.interval.start) * cell_width;
            const int second_x = left_margin + static_cast<int>(data_.alignments[second].interval.start - data_.interval.start) * cell_width;
            const int first_y = top_margin + rows_[first] * row_height + row_height / 2;
            const int second_y = top_margin + rows_[second] * row_height + row_height / 2;
            painter.drawLine(first_x, first_y, second_x, second_y);
        }
    }

    for (const auto index : order_) {
        const auto& alignment = data_.alignments[index];
        const int y = top_margin + rows_[index] * row_height;
        painter.setPen(QColor(70, 70, 70));
        painter.drawText(2, y + 12, QString::fromStdString(alignment.name).left(11));
        const auto fill = strand_color(reverse_strand(alignment));
        auto reference_position = alignment.interval.start;
        std::size_t read_position = 0;
        for (const auto& operation : alignment.cigar) {
            const auto operation_length = static_cast<std::size_t>(operation.length);
            if (operation.operation == 'M' || operation.operation == '=' || operation.operation == 'X') {
                for (std::size_t offset = 0; offset < operation_length && read_position + offset < alignment.sequence.size(); ++offset) {
                    const auto position = reference_position + static_cast<std::int64_t>(offset);
                    if (position < data_.interval.start || position >= data_.interval.end) continue;
                    const int x = left_margin + static_cast<int>(position - data_.interval.start) * cell_width;
                    const char base = static_cast<char>(std::toupper(static_cast<unsigned char>(alignment.sequence[read_position + offset])));
                    const int quality = read_position + offset < alignment.qualities.size()
                        ? static_cast<int>(static_cast<unsigned char>(alignment.qualities[read_position + offset])) - 33
                        : 0;
                    const bool low_quality = quality < data_.minimum_base_quality;
                    const bool mismatch = data_.has_reference && base != reference_at(data_, position);
                    painter.fillRect(x, y, cell_width - 1, row_height - 2,
                                     low_quality ? QColor(225, 225, 225) : mismatch ? QColor(255, 222, 89) : fill);
                    painter.setPen(low_quality ? QColor(130, 130, 130) : Qt::black);
                    painter.drawText(x + 2, y + 12, QString(base));
                }
                reference_position += static_cast<std::int64_t>(operation_length);
                read_position += operation_length;
            } else if (operation.operation == 'I') {
                const auto position = reference_position - 1;
                if (position >= data_.interval.start && position < data_.interval.end) {
                    const int x = left_margin + static_cast<int>(position - data_.interval.start) * cell_width + cell_width - 2;
                    painter.setPen(QPen(QColor(24, 132, 87), 2));
                    painter.drawLine(x, y + 2, x, y + row_height - 4);
                }
                read_position += operation_length;
            } else if (operation.operation == 'D' || operation.operation == 'N') {
                const auto start = std::max(reference_position, data_.interval.start);
                const auto end = std::min(reference_position + static_cast<std::int64_t>(operation_length), data_.interval.end);
                if (start < end) {
                    const int x = left_margin + static_cast<int>(start - data_.interval.start) * cell_width;
                    painter.setPen(QPen(QColor(24, 132, 87), 1));
                    painter.drawLine(x, y + row_height / 2, x + static_cast<int>(end - start) * cell_width, y + row_height / 2);
                }
                reference_position += static_cast<std::int64_t>(operation_length);
            } else if (operation.operation == 'S') {
                read_position += operation_length;
            }
        }
    }
}

}  // namespace bamseek
