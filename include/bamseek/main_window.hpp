#pragma once

#include <bamseek/evidence.hpp>

#include <QFutureWatcher>
#include <QMainWindow>

class QComboBox;
class QCheckBox;
class QLineEdit;
class QPlainTextEdit;
class QPushButton;
class QTableWidget;
class QLabel;

namespace bamseek {

class MainWindow final : public QMainWindow {
public:
    MainWindow();

protected:
    void dragEnterEvent(QDragEnterEvent* event) override;
    void dropEvent(QDropEvent* event) override;

private:
    void choose_bam();
    void choose_reference();
    void run_queries();
    void show_results();
    void show_read_details(int row, int column);
    void export_audit();
    [[nodiscard]] FilterSettings filters() const;

    QLineEdit* bam_path_{};
    QLineEdit* reference_path_{};
    QComboBox* molecule_mode_{};
    QLineEdit* molecule_tag_{};
    QLineEdit* vaf_{};
    QLineEdit* minimum_alt_reads_{};
    QLineEdit* minimum_alt_molecules_{};
    QLineEdit* minimum_mapq_{};
    QLineEdit* minimum_baseq_{};
    QPlainTextEdit* query_text_{};
    QPlainTextEdit* read_details_{};
    QTableWidget* results_{};
    QLabel* status_{};
    QPushButton* run_button_{};
    QCheckBox* include_duplicates_{};
    QCheckBox* include_secondary_{};
    QCheckBox* include_supplementary_{};
    BatchEvidence last_batch_;
    FilterSettings last_filters_;
    QFutureWatcher<BatchEvidence> watcher_;
};

}  // namespace bamseek
