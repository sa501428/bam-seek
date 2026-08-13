#pragma once

#include <bamseek/evidence.hpp>

#include <QFutureWatcher>
#include <QMainWindow>
#include <QString>
#include <QStringList>

class QCheckBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QPlainTextEdit;
class QPushButton;
class QSpinBox;
class QTableWidget;
class QTabWidget;

namespace bamseek {

class IgvCommandReceiver;

class MainWindow final : public QMainWindow {
public:
    MainWindow();

protected:
    void dragEnterEvent(QDragEnterEvent* event) override;
    void dropEvent(QDropEvent* event) override;

private:
    void choose_bams();
    void load_bams();
    void clear_loaded_bams();
    void refresh_loaded_bams();
    void run_queries();
    void show_results();
    void show_result_details(int row, int column = 0);
    void copy_summary();
    void export_audit();
    void audit_saved();
    void set_receiver_enabled(bool enabled);
    void append_received_command(const QString& description);
    [[nodiscard]] FilterSettings filters() const;

    QPlainTextEdit* bam_path_{};
    QListWidget* loaded_bams_{};
    QLineEdit* minimum_mapq_{};
    QLineEdit* minimum_baseq_{};
    QPlainTextEdit* query_text_{};
    QPlainTextEdit* summary_text_{};
    QPlainTextEdit* read_details_{};
    QTableWidget* results_{};
    QTabWidget* tabs_{};
    QLabel* status_{};
    QPushButton* run_button_{};
    QPushButton* copy_summary_button_{};
    QPushButton* export_button_{};
    QPushButton* load_bams_button_{};
    QPushButton* clear_bams_button_{};
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
    QString last_query_text_;
    struct MultiBamBatch {
        BatchEvidence batch;
        QStringList result_bams;
    };
    QFutureWatcher<MultiBamBatch> watcher_;
    struct AuditSave {
        QString path;
        std::string error;
    };
    QFutureWatcher<AuditSave> audit_watcher_;
};

}  // namespace bamseek
