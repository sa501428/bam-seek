#pragma once

#include <bamseek/evidence.hpp>

#include <QFutureWatcher>
#include <QMainWindow>
#include <QString>
#include <QStringList>

#include <memory>
#include <vector>

class QCheckBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QPlainTextEdit;
class QPushButton;
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
    void bams_loaded();
    void remove_selected_bams();
    void clear_loaded_bams();
    void refresh_loaded_bams();
    void run_queries();
    void show_results();
    void show_result_summary(int row, int column = 0);
    void copy_summary();
    void toggle_theme();
    void apply_theme();
    void set_receiver_enabled(bool enabled);
    void append_received_command(const QString& description);
    [[nodiscard]] FilterSettings filters() const;
    [[nodiscard]] bool busy() const;

    QPlainTextEdit* bam_path_{};
    QListWidget* loaded_bams_{};
    QLabel* bam_count_{};
    QLineEdit* minimum_mapq_{};
    QLineEdit* minimum_baseq_{};
    QPlainTextEdit* query_text_{};
    QPlainTextEdit* summary_text_{};
    QTableWidget* results_{};
    QTabWidget* tabs_{};
    QLabel* status_{};
    QPushButton* run_button_{};
    QPushButton* copy_summary_button_{};
    QPushButton* theme_button_{};
    QPushButton* load_bams_button_{};
    QPushButton* remove_bams_button_{};
    QPushButton* clear_bams_button_{};
    QPlainTextEdit* broadcast_text_{};
    QCheckBox* receiver_enabled_{};
    QLineEdit* receiver_port_{};
    QLabel* receiver_status_{};
    IgvCommandReceiver* command_receiver_{};
    bool dark_mode_ = true;

    BatchEvidence last_batch_;
    FilterSettings last_filters_;
    QStringList loaded_bam_paths_;
    QStringList loaded_index_paths_;
    std::vector<std::shared_ptr<EvidenceEngine>> loaded_engines_;
    QStringList result_bam_paths_;

    struct BamLoadBatch {
        QStringList bams;
        QStringList indexes;
        std::vector<std::shared_ptr<EvidenceEngine>> engines;
        QStringList issues;
    };
    QFutureWatcher<BamLoadBatch> bam_load_watcher_;

    struct MultiBamBatch {
        BatchEvidence batch;
        QStringList result_bams;
    };
    QFutureWatcher<MultiBamBatch> watcher_;
};

}  // namespace bamseek
