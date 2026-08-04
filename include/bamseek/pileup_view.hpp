#pragma once

#include <bamseek/evidence.hpp>

#include <QWidget>

#include <vector>

namespace bamseek {

class PileupView final : public QWidget {
public:
    explicit PileupView(QWidget* parent = nullptr);

    void set_data(PileupData data);
    void set_group_pairs(bool enabled);
    [[nodiscard]] QSize sizeHint() const override;

protected:
    void paintEvent(QPaintEvent* event) override;

private:
    void rebuild_layout();

    PileupData data_;
    bool group_pairs_ = true;
    std::vector<std::size_t> order_;
    std::vector<int> rows_;
};

}  // namespace bamseek
