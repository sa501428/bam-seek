#pragma once

#include <bamseek/evidence.hpp>

#include <QFutureWatcher>
#include <QMainWindow>
#include <QString>
#include <QStringList>

class QComboBox;
class QCheckBox;
class QLineEdit;
class QListWidget;
class QPlainTextEdit;
class QPushButton;
class QTableWidget;
class QLabel;
class QTabWidget;
class QSpinBox;
class QScrollArea;

namespace bamseek {

class PileupView;
class IgvCommandReceiver;

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
    void choose_clinical_mapping();
    void load_bams();
    void clear_loaded_bams();
    void refresh_loaded_bams();
    void run_queries();
    void show_results();
    void show_read_details(int row, int column);
    void export_audit();
    void show_pileup();
    void pileup_loaded();
    void audit_saved();
    void set_receiver_enabled(bool enabled);
    void append_received_command(const QString& description);
    [[nodiscard]] FilterSettings filters() const;

    QPlainTextEdit* bam_path_{};
    QPlainTextEdit* index_path_{};
    QListWidget* loaded_bams_{};
    QLineEdit* reference_path_{};
    QLineEdit* clinical_mapping_path_{};
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
    QScrollArea* pileup_scroll_{};
    QLabel* pileup_summary_{};
    QLabel* status_{};
    QPushButton* run_button_{};
    QPushButton* pileup_button_{};
    QPushButton* export_button_{};
    QPushButton* load_bams_button_{};
    QPushButton* clear_bams_button_{};
    QCheckBox* include_duplicates_{};
    QCheckBox* include_secondary_{};
    QCheckBox* include_supplementary_{};
    QCheckBox* group_pairs_{};
    QPlainTextEdit* broadcast_text_{};
    QCheckBox* receiver_enabled_{};
    QSpinBox* receiver_port_{};
    QLabel* receiver_status_{};
    IgvCommandReceiver* command_receiver_{};
    BatchEvidence last_batch_;
    FilterSettings last_filters_;
    QStringList loaded_bam_paths_;
    QStringList loaded_index_paths_;
    QStringList configured_bam_paths_;
    QStringList configured_index_paths_;
    QStringList result_bam_paths_;
    QStringList result_index_paths_;
    QString last_reference_path_;
    QString last_clinical_mapping_path_;
    QString last_query_text_;
    struct MultiBamBatch {
        BatchEvidence batch;
        QStringList result_bams;
        QStringList result_indexes;
    };
    QFutureWatcher<MultiBamBatch> watcher_;
    struct PileupLoad {
        PileupData data;
        QString summary;
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
