#pragma once

#include <bamseek/evidence.hpp>

#include <QFutureWatcher>
#include <QMainWindow>
#include <QString>

class QComboBox;
class QCheckBox;
class QLineEdit;
class QPlainTextEdit;
class QPushButton;
class QTableWidget;
class QLabel;
class QTabWidget;

namespace bamseek {

class PileupView;

class MainWindow final : public QMainWindow {
public:
    MainWindow();

protected:
    void dragEnterEvent(QDragEnterEvent* event) override;
    void dropEvent(QDropEvent* event) override;

private:
    void choose_bam();
    void choose_index();
    void choose_reference();
    void run_queries();
    void show_results();
    void show_read_details(int row, int column);
    void export_audit();
    void show_pileup();
    void pileup_loaded();
    void audit_saved();
    [[nodiscard]] FilterSettings filters() const;

    QLineEdit* bam_path_{};
    QLineEdit* index_path_{};
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
    QTabWidget* tabs_{};
    PileupView* pileup_view_{};
    QLabel* status_{};
    QPushButton* run_button_{};
    QPushButton* pileup_button_{};
    QPushButton* export_button_{};
    QCheckBox* include_duplicates_{};
    QCheckBox* include_secondary_{};
    QCheckBox* include_supplementary_{};
    QCheckBox* group_pairs_{};
    BatchEvidence last_batch_;
    FilterSettings last_filters_;
    QString last_bam_path_;
    QString last_index_path_;
    QString last_reference_path_;
    QString last_query_text_;
    QFutureWatcher<BatchEvidence> watcher_;
    struct PileupLoad {
        PileupData data;
        std::string error;
    };
    QFutureWatcher<PileupLoad> pileup_watcher_;
    struct AuditSave {
        QString path;
        std::string error;
    };
    QFutureWatcher<AuditSave> audit_watcher_;
};

}  // namespace bamseek
