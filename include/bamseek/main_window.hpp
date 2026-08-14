#pragma once

#include <bamseek/comparison.hpp>

#include <QFutureWatcher>
#include <QMainWindow>
#include <QString>
#include <QStringList>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class QCheckBox;
class QComboBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QListWidgetItem;
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
    void start_bam_preparation(const QStringList& pending, bool clear_manual_input);
    void bams_loaded();
    void enqueue_received_bams(const QStringList& paths);
    void load_queued_receiver_bams();
    void remove_selected_bams();
    void clear_loaded_bams();
    void refresh_loaded_bams();
    void set_current_bam(QListWidgetItem* changed_item);
    void clear_results();
    void run_queries();
    void show_results();
    void render_results();
    void copy_summary();
    void toggle_theme();
    void apply_theme();
    void set_receiver_enabled(bool enabled);
    void save_analysis_settings();
    void append_received_command(const QString& description);
    [[nodiscard]] FilterSettings filters() const;
    [[nodiscard]] bool busy() const;

    QPlainTextEdit* bam_path_{};
    QListWidget* loaded_bams_{};
    QLabel* bam_count_{};
    QLineEdit* minimum_mapq_{};
    QLineEdit* minimum_baseq_{};
    QCheckBox* include_duplicates_{};
    QCheckBox* include_secondary_{};
    QCheckBox* include_supplementary_{};
    QComboBox* molecule_mode_{};
    QLineEdit* molecule_tag_{};
    QPlainTextEdit* current_query_text_{};
    QPlainTextEdit* historical_query_text_{};
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
    QLineEdit* upstream_port_{};
    QLabel* receiver_status_{};
    IgvCommandReceiver* command_receiver_{};
    bool dark_mode_ = true;

    BatchEvidence last_batch_;
    FilterSettings last_filters_;
    QStringList loaded_bam_paths_;
    QStringList loaded_index_paths_;
    QStringList preparing_bam_paths_;
    QStringList queued_receiver_bams_;
    QString current_bam_path_;
    std::vector<std::shared_ptr<EvidenceEngine>> loaded_engines_;
    QStringList result_bam_paths_;
    std::vector<VariantOrigin> result_variant_origins_;
    std::unordered_map<std::string, std::unordered_map<std::string, VariantEvidence>> evidence_cache_;
    bool clear_bam_input_after_load_ = false;
    bool refreshing_bam_list_ = false;

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
        std::vector<VariantOrigin> result_variant_origins;
        std::unordered_map<std::string, std::unordered_map<std::string, VariantEvidence>> cache_updates;
        std::size_t cache_hits{};
        std::size_t calculations{};
    };
    QFutureWatcher<MultiBamBatch> watcher_;
};

}  // namespace bamseek
