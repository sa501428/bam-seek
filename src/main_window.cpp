#include <bamseek/main_window.hpp>

#include <bamseek/igv_command_receiver.hpp>
#include <bamseek/query.hpp>

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QDateTime>
#include <QDragEnterEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontDatabase>
#include <QFrame>
#include <QFuture>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QIntValidator>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QMimeData>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSettings>
#include <QSize>
#include <QSpinBox>
#include <QSplitter>
#include <QTabWidget>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QTabBar>
#include <QUrl>
#include <QVBoxLayout>
#include <QWidget>
#include <QtConcurrentRun>

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace bamseek {
namespace {

QString display_query(const VariantQuery& query) {
    const auto genomic = QString::fromStdString(query.contig) + ':' + QString::number(query.position + 1) + ' '
        + QString::fromStdString(query.reference) + '>' + QString::fromStdString(query.alternate);
    if (query.gene.empty()) return genomic;
    QString clinical = QString::fromStdString(query.gene);
    if (!query.coding_change.empty()) {
        clinical += ' ';
        if (!query.transcript.empty()) clinical += QString::fromStdString(query.transcript) + ':';
        clinical += QString::fromStdString(query.coding_change);
    }
    if (!query.protein_change.empty()) clinical += ' ' + QString::fromStdString(query.protein_change);
    return clinical + "  ·  " + genomic;
}

QString sanitized_resource_uri(const QString& path) {
    if (!path.startsWith("https://")) return path;
    QUrl sanitized(path);
    sanitized.setUserName({});
    sanitized.setPassword({});
    sanitized.setQuery({});
    sanitized.setFragment({});
    return sanitized.toString();
}

QString resource_label(const QString& path) {
    const auto name = path.startsWith("https://") ? QFileInfo(QUrl(path).path()).fileName() : QFileInfo(path).fileName();
    return name.isEmpty() ? sanitized_resource_uri(path) : name;
}

std::string sanitized_error(std::string message, const QString& resource) {
    auto sanitized = QString::fromStdString(message);
    if (resource.startsWith("https://")) sanitized.replace(resource, sanitized_resource_uri(resource));
    return sanitized.toStdString();
}

QString index_path_for(const QString& bam) {
    const QFileInfo info(bam);
    const QStringList candidates{
        bam + ".bai",
        bam + ".csi",
        info.absolutePath() + '/' + info.completeBaseName() + ".bai",
        info.absolutePath() + '/' + info.completeBaseName() + ".csi",
    };
    for (const auto& candidate : candidates) {
        if (QFileInfo::exists(candidate)) return candidate;
    }
    return {};
}

QStringList nonempty_lines(const QPlainTextEdit* edit) {
    QStringList result;
    for (auto line : edit->toPlainText().split('\n')) {
        line = line.trimmed();
        if (!line.isEmpty()) result.append(line);
    }
    return result;
}

void append_resource_lines(QPlainTextEdit* edit, const QStringList& paths) {
    auto existing = edit->toPlainText();
    if (!existing.isEmpty() && !existing.endsWith('\n')) existing += '\n';
    edit->setPlainText(existing + paths.join('\n'));
    auto cursor = edit->textCursor();
    cursor.movePosition(QTextCursor::End);
    edit->setTextCursor(cursor);
}

QString percent(const double fraction) {
    return QString::number(fraction * 100.0, 'f', 3) + '%';
}

QString evidence_summary(const VariantEvidence& evidence, const QString& bam_label) {
    const auto& counts = evidence.counts;
    return QString("In %1, %2 was supported by %3 of %4 informative reads (VAF %5; %6 REF and %7 OTHER/N reads). "
                   "After collapsing reads with the same read name into paired fragments, %8 of %9 informative molecules supported the variant "
                   "(molecule VAF %10; %11 ambiguous fragments). ALT support included %12 forward and %13 reverse reads.")
        .arg(bam_label, display_query(evidence.query))
        .arg(counts.alternate_reads)
        .arg(counts.informative_read_depth())
        .arg(percent(counts.allele_fraction()))
        .arg(counts.reference_reads)
        .arg(counts.other_reads)
        .arg(evidence.molecule_counts_available ? QString::number(counts.alternate_molecules) : "N/A")
        .arg(evidence.molecule_counts_available ? QString::number(counts.molecule_depth()) : "N/A")
        .arg(evidence.molecule_counts_available ? percent(counts.molecule_allele_fraction()) : "unavailable")
        .arg(evidence.molecule_counts_available ? QString::number(counts.other_molecules) : "N/A")
        .arg(counts.alternate_forward_reads)
        .arg(counts.alternate_reverse_reads);
}

QLabel* section_title(const QString& text, QWidget* parent) {
    auto* label = new QLabel(text, parent);
    label->setObjectName("SectionTitle");
    return label;
}

QLabel* muted_label(const QString& text, QWidget* parent) {
    auto* label = new QLabel(text, parent);
    label->setObjectName("MutedLabel");
    label->setWordWrap(true);
    return label;
}

QString premium_style(const bool dark) {
    QString style = R"QSS(
QMainWindow, QWidget#AppRoot {
    background: @BG@;
    color: @TEXT@;
}
QWidget {
    font-size: 13px;
    color: @TEXT@;
}
QTabWidget#WorkspaceTabs::pane {
    border: 0;
    background: @BG@;
}
QTabWidget#WorkspaceTabs QTabBar::tab {
    background: transparent;
    color: @MUTED@;
    border: 0;
    border-bottom: 2px solid transparent;
    padding: 13px 24px;
    margin-right: 4px;
    font-weight: 600;
}
QTabWidget#WorkspaceTabs QTabBar::tab:selected {
    color: @TAB_SELECTED@;
    border-bottom: 2px solid #7c6cff;
}
QTabWidget#WorkspaceTabs QTabBar::tab:hover { color: @TAB_HOVER@; }
QFrame#Card {
    background: @CARD@;
    border: 1px solid @BORDER@;
    border-radius: 14px;
}
QLabel#SectionTitle { font-size: 16px; font-weight: 650; color: @TAB_SELECTED@; }
QLabel#MutedLabel { color: @MUTED@; }
QLabel#Pill {
    color: @PILL_TEXT@;
    background: @PILL_BG@;
    border: 1px solid @PILL_BORDER@;
    border-radius: 10px;
    padding: 3px 10px;
    font-weight: 600;
}
QPlainTextEdit, QLineEdit, QListWidget, QTableWidget, QSpinBox {
    background: @FIELD@;
    color: @TEXT@;
    border: 1px solid @FIELD_BORDER@;
    border-radius: 9px;
    selection-background-color: @SELECT@;
    selection-color: @SELECT_TEXT@;
}
QPlainTextEdit, QLineEdit, QListWidget { padding: 8px; }
QPlainTextEdit:focus, QLineEdit:focus, QListWidget:focus, QTableWidget:focus, QSpinBox:focus {
    border: 1px solid #796cff;
}
QPlainTextEdit#SummaryBox {
    background: @SUMMARY@;
    border: 1px solid @SUMMARY_BORDER@;
    padding: 12px;
}
QListWidget::item { border-radius: 7px; padding: 6px; margin: 2px 0; }
QListWidget::item:selected { background: @LIST_SELECTED@; color: @SELECT_TEXT@; }
QPushButton {
    background: @BUTTON@;
    color: @BUTTON_TEXT@;
    border: 1px solid @BUTTON_BORDER@;
    border-radius: 9px;
    padding: 8px 14px;
    font-weight: 600;
}
QPushButton:hover { background: @BUTTON_HOVER@; border-color: @BUTTON_HOVER_BORDER@; }
QPushButton:pressed { background: @BUTTON_PRESSED@; }
QPushButton:disabled { color: @DISABLED_TEXT@; background: @DISABLED_BG@; border-color: @DISABLED_BORDER@; }
QPushButton#PrimaryButton {
    color: #ffffff;
    background: #6f5cf6;
    border: 1px solid #8575ff;
    padding: 10px 20px;
}
QPushButton#PrimaryButton:hover { background: #7d6bff; }
QPushButton#PrimaryButton:disabled {
    color: @DISABLED_TEXT@;
    background: @DISABLED_BG@;
    border-color: @DISABLED_BORDER@;
}
QPushButton#DangerButton { color: @DANGER@; }
QPushButton#DangerButton:disabled {
    color: @DISABLED_TEXT@;
    background: @DISABLED_BG@;
    border-color: @DISABLED_BORDER@;
}
QPushButton#ThemeButton { padding: 5px 11px; margin: 5px 2px 5px 8px; }
QHeaderView::section {
    background: @HEADER@;
    color: @HEADER_TEXT@;
    border: 0;
    border-bottom: 1px solid @HEADER_BORDER@;
    padding: 9px 8px;
    font-weight: 600;
}
QTableWidget {
    gridline-color: @HEADER_BORDER@;
    alternate-background-color: @ALT_ROW@;
}
QTableWidget::item { padding: 6px; }
QTableWidget::item:selected { background: @TABLE_SELECT@; color: @SELECT_TEXT@; }
QScrollBar:vertical { background: transparent; width: 10px; margin: 2px; }
QScrollBar::handle:vertical { background: @SCROLL@; border-radius: 5px; min-height: 28px; }
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QSplitter::handle { background: transparent; width: 10px; height: 10px; }
QCheckBox { spacing: 8px; color: @CHECK_TEXT@; }
QCheckBox::indicator { width: 16px; height: 16px; }
QToolTip { color: @TOOLTIP_TEXT@; background: @TOOLTIP_BG@; border: 1px solid @TOOLTIP_BORDER@; padding: 5px; }
)QSS";
    const auto replace = [&style](const char* token, const char* dark_color, const char* light_color, const bool use_dark) {
        style.replace(token, use_dark ? dark_color : light_color);
    };
    replace("@BG@", "#0a0e16", "#f4f6fa", dark);
    replace("@TEXT@", "#e8edf6", "#172033", dark);
    replace("@MUTED@", "#8f9bb0", "#667085", dark);
    replace("@TAB_SELECTED@", "#f4f7fb", "#172033", dark);
    replace("@TAB_HOVER@", "#cfd7e6", "#344054", dark);
    replace("@CARD@", "#111824", "#ffffff", dark);
    replace("@BORDER@", "#222d3d", "#d9e0ea", dark);
    replace("@PILL_TEXT@", "#b9afff", "#5b48d6", dark);
    replace("@PILL_BG@", "#211d43", "#eeebff", dark);
    replace("@PILL_BORDER@", "#39316d", "#d5ceff", dark);
    replace("@FIELD@", "#0c121c", "#f8fafc", dark);
    replace("@FIELD_BORDER@", "#273244", "#ccd5e1", dark);
    replace("@SELECT@", "#5146a8", "#dcd8ff", dark);
    replace("@SELECT_TEXT@", "#ffffff", "#172033", dark);
    replace("@SUMMARY@", "#0e1521", "#fbfcfe", dark);
    replace("@SUMMARY_BORDER@", "#313d52", "#d6dee9", dark);
    replace("@LIST_SELECTED@", "#242d47", "#e9e6ff", dark);
    replace("@BUTTON@", "#1a2331", "#ffffff", dark);
    replace("@BUTTON_TEXT@", "#dce3ef", "#344054", dark);
    replace("@BUTTON_BORDER@", "#303c4f", "#cfd7e3", dark);
    replace("@BUTTON_HOVER@", "#222d3d", "#f5f7fb", dark);
    replace("@BUTTON_HOVER_BORDER@", "#46546b", "#aeb9c8", dark);
    replace("@BUTTON_PRESSED@", "#151d29", "#e9edf3", dark);
    replace("@DISABLED_TEXT@", "#586477", "#98a2b3", dark);
    replace("@DISABLED_BG@", "#111722", "#f1f3f6", dark);
    replace("@DISABLED_BORDER@", "#202938", "#e1e5eb", dark);
    replace("@DANGER@", "#f1a5ad", "#c24155", dark);
    replace("@HEADER@", "#151e2b", "#f2f5f9", dark);
    replace("@HEADER_TEXT@", "#9eabc0", "#667085", dark);
    replace("@HEADER_BORDER@", "#2a3547", "#dce2eb", dark);
    replace("@ALT_ROW@", "#0f1621", "#fafbfc", dark);
    replace("@TABLE_SELECT@", "#262e52", "#e9e6ff", dark);
    replace("@SCROLL@", "#344055", "#c3cad5", dark);
    replace("@CHECK_TEXT@", "#cbd3e1", "#344054", dark);
    replace("@TOOLTIP_TEXT@", "#edf2fa", "#ffffff", dark);
    replace("@TOOLTIP_BG@", "#17202d", "#172033", dark);
    replace("@TOOLTIP_BORDER@", "#354156", "#344054", dark);
    return style;
}

}  // namespace

MainWindow::MainWindow() {
    setWindowTitle("BAM Seek — variant allele frequency");
    setAcceptDrops(true);
    resize(1320, 900);
    setMinimumSize(1040, 720);
    setFont(QFontDatabase::systemFont(QFontDatabase::GeneralFont));
    dark_mode_ = QSettings().value("appearance/dark_mode", true).toBool();
    setStyleSheet(premium_style(dark_mode_));

    auto* root = new QWidget(this);
    root->setObjectName("AppRoot");
    auto* root_layout = new QVBoxLayout(root);
    root_layout->setContentsMargins(20, 14, 20, 18);
    root_layout->setSpacing(0);

    tabs_ = new QTabWidget(root);
    tabs_->setObjectName("WorkspaceTabs");
    tabs_->setDocumentMode(true);
    tabs_->tabBar()->setDrawBase(false);
    theme_button_ = new QPushButton(tabs_);
    theme_button_->setObjectName("ThemeButton");
    tabs_->setCornerWidget(theme_button_, Qt::TopRightCorner);
    apply_theme();

    auto* analysis_page = new QWidget(tabs_);
    auto* analysis_layout = new QVBoxLayout(analysis_page);
    analysis_layout->setContentsMargins(4, 14, 4, 4);
    analysis_layout->setSpacing(14);

    auto* input_splitter = new QSplitter(Qt::Horizontal, analysis_page);
    input_splitter->setChildrenCollapsible(false);

    auto* bam_card = new QFrame(input_splitter);
    bam_card->setObjectName("Card");
    auto* bam_layout = new QVBoxLayout(bam_card);
    bam_layout->setContentsMargins(18, 16, 18, 18);
    bam_layout->setSpacing(11);
    auto* bam_header = new QHBoxLayout();
    bam_header->addWidget(section_title("BAM sources", bam_card));
    bam_header->addStretch(1);
    bam_count_ = new QLabel("0 ready", bam_card);
    bam_count_->setObjectName("Pill");
    bam_header->addWidget(bam_count_);
    bam_layout->addLayout(bam_header);
    bam_path_ = new QPlainTextEdit(bam_card);
    bam_path_->setPlaceholderText("/data/sample.bam\nhttps://example.org/sample.bam");
    bam_path_->setMaximumHeight(72);
    bam_layout->addWidget(bam_path_);
    auto* bam_add_row = new QHBoxLayout();
    auto* browse_bam = new QPushButton("Choose files…", bam_card);
    load_bams_button_ = new QPushButton("Add BAMs", bam_card);
    load_bams_button_->setObjectName("PrimaryButton");
    bam_add_row->addWidget(browse_bam);
    bam_add_row->addStretch(1);
    bam_add_row->addWidget(load_bams_button_);
    bam_layout->addLayout(bam_add_row);
    auto* bam_manage_row = new QHBoxLayout();
    bam_manage_row->addWidget(muted_label("Prepared BAMs", bam_card));
    bam_manage_row->addStretch(1);
    remove_bams_button_ = new QPushButton("Remove", bam_card);
    remove_bams_button_->setEnabled(false);
    clear_bams_button_ = new QPushButton("Clear", bam_card);
    clear_bams_button_->setObjectName("DangerButton");
    clear_bams_button_->setEnabled(false);
    bam_manage_row->addWidget(remove_bams_button_);
    bam_manage_row->addWidget(clear_bams_button_);
    bam_layout->addLayout(bam_manage_row);
    loaded_bams_ = new QListWidget(bam_card);
    loaded_bams_->setAlternatingRowColors(false);
    loaded_bams_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    loaded_bams_->setMinimumHeight(105);
    loaded_bams_->setToolTip("Prepared BAMs are reused when VAF calculation is requested.");
    bam_layout->addWidget(loaded_bams_, 1);

    auto* variant_card = new QFrame(input_splitter);
    variant_card->setObjectName("Card");
    auto* variant_layout = new QVBoxLayout(variant_card);
    variant_layout->setContentsMargins(18, 16, 18, 18);
    variant_layout->setSpacing(11);
    variant_layout->addWidget(section_title("Variants", variant_card));
    query_text_ = new QPlainTextEdit(variant_card);
    query_text_->setPlaceholderText(
        "chr7:140453136 A>T\n"
        "chr7:g.140453136A>T\n"
        "BRAF c.1799T>A p.V600E chr7:140453136 A>T");
    variant_layout->addWidget(query_text_, 1);
    auto* filter_row = new QHBoxLayout();
    filter_row->addWidget(muted_label("mapQ ≥", variant_card));
    minimum_mapq_ = new QLineEdit("20", variant_card);
    minimum_mapq_->setMaximumWidth(58);
    minimum_mapq_->setAlignment(Qt::AlignCenter);
    minimum_mapq_->setValidator(new QIntValidator(0, 255, minimum_mapq_));
    filter_row->addWidget(minimum_mapq_);
    filter_row->addSpacing(12);
    filter_row->addWidget(muted_label("baseQ ≥", variant_card));
    minimum_baseq_ = new QLineEdit("20", variant_card);
    minimum_baseq_->setMaximumWidth(58);
    minimum_baseq_->setAlignment(Qt::AlignCenter);
    minimum_baseq_->setValidator(new QIntValidator(0, 255, minimum_baseq_));
    filter_row->addWidget(minimum_baseq_);
    filter_row->addStretch(1);
    variant_layout->addLayout(filter_row);

    input_splitter->addWidget(bam_card);
    input_splitter->addWidget(variant_card);
    input_splitter->setSizes({580, 700});
    analysis_layout->addWidget(input_splitter, 1);

    auto* action_row = new QHBoxLayout();
    status_ = new QLabel("Add BAMs and variants to begin.", analysis_page);
    status_->setObjectName("MutedLabel");
    run_button_ = new QPushButton("Calculate VAFs", analysis_page);
    run_button_->setObjectName("PrimaryButton");
    run_button_->setEnabled(false);
    action_row->addWidget(status_, 1);
    action_row->addWidget(run_button_);
    analysis_layout->addLayout(action_row);

    auto* results_card = new QFrame(analysis_page);
    results_card->setObjectName("Card");
    auto* results_layout = new QVBoxLayout(results_card);
    results_layout->setContentsMargins(16, 14, 16, 16);
    results_layout->setSpacing(10);
    results_layout->addWidget(section_title("VAF results", results_card));
    results_ = new QTableWidget(results_card);
    results_->setColumnCount(11);
    results_->setHorizontalHeaderLabels({"BAM", "Variant", "ALT reads", "Read depth", "VAF",
        "ALT molecules", "Molecule depth", "Molecule VAF", "OTHER/N", "Ambiguous", "ALT F / R"});
    results_->setAlternatingRowColors(true);
    results_->setSelectionBehavior(QAbstractItemView::SelectRows);
    results_->setSelectionMode(QAbstractItemView::SingleSelection);
    results_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    results_->setShowGrid(false);
    results_->verticalHeader()->setVisible(false);
    results_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    results_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    results_->setMinimumHeight(175);
    results_layout->addWidget(results_, 1);
    auto* summary_header = new QHBoxLayout();
    summary_header->addWidget(muted_label("Summary", results_card));
    summary_header->addStretch(1);
    copy_summary_button_ = new QPushButton("Copy summary", results_card);
    copy_summary_button_->setEnabled(false);
    summary_header->addWidget(copy_summary_button_);
    results_layout->addLayout(summary_header);
    summary_text_ = new QPlainTextEdit(results_card);
    summary_text_->setObjectName("SummaryBox");
    summary_text_->setReadOnly(true);
    summary_text_->setMaximumHeight(90);
    summary_text_->setPlaceholderText("Select a result to generate a copy-ready paragraph.");
    results_layout->addWidget(summary_text_);
    analysis_layout->addWidget(results_card, 2);

    auto* broadcast_page = new QWidget(tabs_);
    auto* broadcast_layout = new QVBoxLayout(broadcast_page);
    broadcast_layout->setContentsMargins(4, 14, 4, 4);
    auto* receiver_card = new QFrame(broadcast_page);
    receiver_card->setObjectName("Card");
    auto* receiver_layout = new QVBoxLayout(receiver_card);
    receiver_layout->setContentsMargins(18, 16, 18, 18);
    auto* receiver_controls = new QHBoxLayout();
    receiver_enabled_ = new QCheckBox("Receiver enabled", receiver_card);
    receiver_enabled_->setChecked(true);
    receiver_port_ = new QSpinBox(receiver_card);
    receiver_port_->setRange(1024, 65535);
    receiver_port_->setValue(60151);
    receiver_port_->setMaximumWidth(100);
    receiver_status_ = new QLabel(receiver_card);
    receiver_status_->setObjectName("MutedLabel");
    receiver_controls->addWidget(receiver_enabled_);
    receiver_controls->addSpacing(16);
    receiver_controls->addWidget(muted_label("Port", receiver_card));
    receiver_controls->addWidget(receiver_port_);
    receiver_controls->addSpacing(16);
    receiver_controls->addWidget(receiver_status_, 1);
    receiver_layout->addLayout(receiver_controls);
    broadcast_text_ = new QPlainTextEdit(receiver_card);
    broadcast_text_->setPlaceholderText("Received /load and /goto commands will appear here…");
    receiver_layout->addWidget(broadcast_text_, 1);
    broadcast_layout->addWidget(receiver_card, 1);

    tabs_->addTab(analysis_page, "VAF workspace");
    tabs_->addTab(broadcast_page, "Broadcast receiver");
    root_layout->addWidget(tabs_);
    setCentralWidget(root);

    connect(browse_bam, &QPushButton::clicked, this, [this] { choose_bams(); });
    connect(load_bams_button_, &QPushButton::clicked, this, [this] { load_bams(); });
    connect(remove_bams_button_, &QPushButton::clicked, this, [this] { remove_selected_bams(); });
    connect(clear_bams_button_, &QPushButton::clicked, this, [this] { clear_loaded_bams(); });
    connect(loaded_bams_, &QListWidget::itemSelectionChanged, this, [this] {
        remove_bams_button_->setEnabled(!busy() && !loaded_bams_->selectedItems().isEmpty());
    });
    connect(run_button_, &QPushButton::clicked, this, [this] { run_queries(); });
    connect(copy_summary_button_, &QPushButton::clicked, this, [this] { copy_summary(); });
    connect(theme_button_, &QPushButton::clicked, this, [this] { toggle_theme(); });
    connect(results_, &QTableWidget::cellClicked, this, [this](const int row, const int column) { show_result_summary(row, column); });
    connect(&bam_load_watcher_, &QFutureWatcher<BamLoadBatch>::finished, this, [this] { bams_loaded(); });
    connect(&watcher_, &QFutureWatcher<MultiBamBatch>::finished, this, [this] { show_results(); });

    command_receiver_ = new IgvCommandReceiver(this);
    connect(receiver_enabled_, &QCheckBox::toggled, this, [this](const bool enabled) { set_receiver_enabled(enabled); });
    connect(receiver_port_, &QSpinBox::valueChanged, this, [this] {
        if (receiver_enabled_->isChecked()) set_receiver_enabled(true);
    });
    connect(command_receiver_, &IgvCommandReceiver::request_received, this, [this](const QString& description) {
        append_received_command(description);
    });
    connect(command_receiver_, &IgvCommandReceiver::listening_changed, this,
        [this](const bool listening, const quint16 port, const QString& error) {
            receiver_status_->setText(listening ? QString("Listening on 127.0.0.1:%1").arg(port)
                                                : (error.isEmpty() ? "Receiver stopped" : "Unavailable · " + error));
        });
    set_receiver_enabled(true);
}

bool MainWindow::busy() const {
    return bam_load_watcher_.isRunning() || watcher_.isRunning();
}

void MainWindow::toggle_theme() {
    dark_mode_ = !dark_mode_;
    QSettings().setValue("appearance/dark_mode", dark_mode_);
    apply_theme();
}

void MainWindow::apply_theme() {
    setStyleSheet(premium_style(dark_mode_));
    if (theme_button_ != nullptr) {
        theme_button_->setText(dark_mode_ ? "Light" : "Dark");
        theme_button_->setToolTip(dark_mode_ ? "Switch to light mode" : "Switch to dark mode");
    }
    if (loaded_bams_ != nullptr) {
        const QColor ready_color(dark_mode_ ? "#b7efcf" : "#157347");
        for (int row = 0; row < loaded_bams_->count(); ++row) loaded_bams_->item(row)->setForeground(ready_color);
    }
    if (results_ != nullptr) {
        const QColor accent(dark_mode_ ? "#b8adff" : "#5b48d6");
        for (int row = 0; row < results_->rowCount(); ++row) {
            if (results_->item(row, 4) != nullptr) results_->item(row, 4)->setForeground(accent);
            if (results_->item(row, 7) != nullptr) results_->item(row, 7)->setForeground(accent);
        }
    }
}

void MainWindow::set_receiver_enabled(const bool enabled) {
    if (!enabled) {
        command_receiver_->close();
        return;
    }
    (void)command_receiver_->listen(static_cast<quint16>(receiver_port_->value()));
}

void MainWindow::append_received_command(const QString& description) {
    auto cursor = broadcast_text_->textCursor();
    cursor.movePosition(QTextCursor::End);
    if (!broadcast_text_->document()->isEmpty()) cursor.insertText("\n\n");
    cursor.insertText(QDateTime::currentDateTime().toString(Qt::ISODate) + "\n" + description);
    broadcast_text_->setTextCursor(cursor);
    broadcast_text_->ensureCursorVisible();
}

void MainWindow::dragEnterEvent(QDragEnterEvent* event) {
    if (event->mimeData()->hasUrls()) event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent* event) {
    QStringList bams;
    for (const auto& url : event->mimeData()->urls()) {
        if (url.isLocalFile() && url.toLocalFile().endsWith(".bam", Qt::CaseInsensitive)) bams.append(url.toLocalFile());
        else if (url.scheme() == "https") bams.append(url.toString());
    }
    if (bams.isEmpty()) return;
    append_resource_lines(bam_path_, bams);
    tabs_->setCurrentIndex(0);
    status_->setText(QString("%1 BAM source(s) ready to add.").arg(bams.size()));
    event->acceptProposedAction();
}

void MainWindow::choose_bams() {
    const auto files = QFileDialog::getOpenFileNames(this, "Choose BAM files", {}, "BAM files (*.bam)");
    if (!files.isEmpty()) append_resource_lines(bam_path_, files);
}

void MainWindow::load_bams() {
    if (busy()) return;
    const auto pending = nonempty_lines(bam_path_);
    if (pending.isEmpty()) {
        QMessageBox::information(this, "No BAM sources", "Add local BAM paths or HTTPS BAM URLs first.");
        return;
    }

    QStringList candidates;
    QStringList candidate_indexes;
    QStringList issues;
    for (const auto& bam : pending) {
        const bool remote = bam.startsWith("https://");
        const QUrl remote_url(bam);
        if ((bam.contains("://") && !remote)
            || (remote && (!remote_url.isValid() || remote_url.host().isEmpty()))
            || (!remote && !bam.endsWith(".bam", Qt::CaseInsensitive))) {
            issues.append(sanitized_resource_uri(bam) + ": expected a local .bam path or valid HTTPS BAM URL.");
            continue;
        }
        if (!remote && (!QFileInfo(bam).exists() || !QFileInfo(bam).isFile())) {
            issues.append(bam + ": BAM file not found.");
            continue;
        }
        const auto index = remote ? QString{} : index_path_for(bam);
        if (!remote && index.isEmpty()) {
            issues.append(resource_label(bam) + ": no adjacent .bai or .csi index was found.");
            continue;
        }
        if (loaded_bam_paths_.contains(bam) || candidates.contains(bam)) {
            issues.append(resource_label(bam) + ": already in the BAM panel.");
            continue;
        }
        candidates.append(bam);
        candidate_indexes.append(index);
    }
    if (candidates.isEmpty()) {
        QMessageBox::warning(this, "No BAMs added", issues.join('\n'));
        return;
    }

    bam_path_->setEnabled(false);
    load_bams_button_->setEnabled(false);
    run_button_->setEnabled(false);
    remove_bams_button_->setEnabled(false);
    clear_bams_button_->setEnabled(false);
    status_->setText(QString("Preparing %1 BAM source(s) in the background…").arg(candidates.size()));
    bam_load_watcher_.setFuture(QtConcurrent::run([candidates, candidate_indexes, issues = std::move(issues)]() mutable {
        BamLoadBatch loaded;
        loaded.issues = std::move(issues);
        for (qsizetype i = 0; i < candidates.size(); ++i) {
            const auto bam = candidates[i];
            const auto index = candidate_indexes[i];
            try {
                igv::Resource resource{.uri = bam.toStdString()};
                if (!index.isEmpty()) resource.index_uri = index.toStdString();
                auto engine = std::make_shared<EvidenceEngine>(std::move(resource));
                loaded.bams.append(bam);
                loaded.indexes.append(index);
                loaded.engines.push_back(std::move(engine));
            } catch (const std::exception& error) {
                loaded.issues.append(resource_label(bam) + ": " + QString::fromStdString(sanitized_error(error.what(), bam)));
            }
        }
        return loaded;
    }));
}

void MainWindow::bams_loaded() {
    const auto loaded = bam_load_watcher_.result();
    for (qsizetype i = 0; i < loaded.bams.size(); ++i) {
        loaded_bam_paths_.append(loaded.bams[i]);
        loaded_index_paths_.append(loaded.indexes[i]);
        loaded_engines_.push_back(loaded.engines[static_cast<std::size_t>(i)]);
    }
    if (!loaded.bams.isEmpty()) {
        bam_path_->clear();
        last_batch_ = {};
        result_bam_paths_.clear();
        results_->setRowCount(0);
        summary_text_->clear();
        copy_summary_button_->setEnabled(false);
    }
    bam_path_->setEnabled(true);
    load_bams_button_->setEnabled(true);
    run_button_->setEnabled(true);
    refresh_loaded_bams();
    status_->setText(QString("%1 BAM(s) ready. No VAF calculation has run yet.").arg(loaded_bam_paths_.size()));
    if (!loaded.issues.isEmpty()) QMessageBox::warning(this, "BAM preparation notes", loaded.issues.join('\n'));
}

void MainWindow::remove_selected_bams() {
    if (busy()) return;
    QList<int> rows;
    for (auto* item : loaded_bams_->selectedItems()) rows.append(loaded_bams_->row(item));
    std::sort(rows.begin(), rows.end(), std::greater<>());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    for (const int row : rows) {
        loaded_bam_paths_.removeAt(row);
        loaded_index_paths_.removeAt(row);
        loaded_engines_.erase(loaded_engines_.begin() + row);
    }
    last_batch_ = {};
    result_bam_paths_.clear();
    results_->setRowCount(0);
    summary_text_->clear();
    copy_summary_button_->setEnabled(false);
    refresh_loaded_bams();
    status_->setText(QString("%1 BAM(s) ready. Results cleared after source changes.").arg(loaded_bam_paths_.size()));
}

void MainWindow::clear_loaded_bams() {
    if (busy()) return;
    loaded_bam_paths_.clear();
    loaded_index_paths_.clear();
    loaded_engines_.clear();
    result_bam_paths_.clear();
    last_batch_ = {};
    results_->setRowCount(0);
    summary_text_->clear();
    copy_summary_button_->setEnabled(false);
    refresh_loaded_bams();
    status_->setText("BAM panel cleared.");
}

void MainWindow::refresh_loaded_bams() {
    loaded_bams_->clear();
    for (qsizetype i = 0; i < loaded_bam_paths_.size(); ++i) {
        const auto& bam = loaded_bam_paths_[i];
        const bool remote = bam.startsWith("https://");
        auto* item = new QListWidgetItem(QString("●  %1\n    %2 · Ready")
            .arg(resource_label(bam), remote ? "HTTPS" : "Local"), loaded_bams_);
        item->setForeground(QColor(dark_mode_ ? "#b7efcf" : "#157347"));
        item->setSizeHint(QSize(0, 48));
        item->setToolTip("BAM: " + sanitized_resource_uri(bam) + "\nIndex: "
            + (loaded_index_paths_[i].isEmpty() ? "automatic remote discovery" : loaded_index_paths_[i]));
    }
    bam_count_->setText(QString("%1 ready").arg(loaded_bam_paths_.size()));
    const bool has_bams = !loaded_bam_paths_.isEmpty();
    run_button_->setEnabled(has_bams && !busy());
    remove_bams_button_->setEnabled(has_bams && !loaded_bams_->selectedItems().isEmpty() && !busy());
    clear_bams_button_->setEnabled(has_bams && !busy());
}

FilterSettings MainWindow::filters() const {
    FilterSettings values;
    values.minimum_mapping_quality = std::max(0, minimum_mapq_->text().toInt());
    values.minimum_base_quality = std::max(0, minimum_baseq_->text().toInt());
    values.minimum_variant_allele_fraction = 0.0;
    values.minimum_alternate_reads = 1;
    values.minimum_alternate_molecules = 1;
    values.molecule_mode = MoleculeMode::raw_reads;
    return values;
}

void MainWindow::run_queries() {
    if (busy()) return;
    if (loaded_engines_.empty()) {
        QMessageBox::warning(this, "No prepared BAMs", "Add at least one BAM to the BAM panel first.");
        return;
    }
    if (!minimum_mapq_->hasAcceptableInput() || !minimum_baseq_->hasAcceptableInput()) {
        QMessageBox::warning(this, "Invalid quality filters", "MapQ and baseQ must be integers from 0 to 255.");
        return;
    }
    auto parsed = parse_queries(query_text_->toPlainText().toStdString());
    std::vector<Query> variants;
    variants.reserve(parsed.queries.size());
    for (auto& query : parsed.queries) {
        if (std::holds_alternative<VariantQuery>(query)) variants.push_back(std::move(query));
        else parsed.errors.push_back(std::get<RegionQuery>(query).source_text
            + ": enter a specific REF>ALT variant rather than a region");
    }
    if (variants.empty()) {
        QMessageBox::warning(this, "No valid variants", QString::fromStdString(parsed.errors.empty()
            ? "Enter at least one variant." : parsed.errors.front()));
        return;
    }

    last_filters_ = filters();
    const auto bams = loaded_bam_paths_;
    const auto engines = loaded_engines_;
    run_button_->setEnabled(false);
    load_bams_button_->setEnabled(false);
    remove_bams_button_->setEnabled(false);
    clear_bams_button_->setEnabled(false);
    status_->setText(QString("Calculating %1 variant(s) across %2 prepared BAM(s)…")
        .arg(variants.size()).arg(bams.size()));
    watcher_.setFuture(QtConcurrent::run([queries = std::move(variants), errors = std::move(parsed.errors),
                                          filter_values = last_filters_, bams, engines] {
        MultiBamBatch combined;
        combined.batch.errors = errors;
        for (std::size_t source_index = 0; source_index < engines.size(); ++source_index) {
            const auto& bam = bams[static_cast<qsizetype>(source_index)];
            try {
                auto evaluated = engines[source_index]->evaluate(queries, filter_values);
                for (auto& result : evaluated.results) {
                    combined.batch.results.push_back(std::move(result));
                    combined.result_bams.append(bam);
                }
                for (auto& error : evaluated.errors) {
                    combined.batch.errors.push_back("[" + resource_label(bam).toStdString() + "] " + error);
                }
            } catch (const std::exception& error) {
                combined.batch.errors.push_back("[" + resource_label(bam).toStdString() + "] "
                    + sanitized_error(error.what(), bam));
            }
        }
        return combined;
    }));
}

void MainWindow::show_results() {
    const auto combined = watcher_.result();
    last_batch_ = combined.batch;
    result_bam_paths_ = combined.result_bams;
    results_->setRowCount(0);
    for (int row = 0; row < static_cast<int>(last_batch_.results.size()); ++row) {
        const auto* evidence = std::get_if<VariantEvidence>(&last_batch_.results[static_cast<std::size_t>(row)]);
        if (evidence == nullptr) continue;
        results_->insertRow(row);
        const auto& count = evidence->counts;
        const QList<QString> values{resource_label(result_bam_paths_[row]), display_query(evidence->query),
            QString::number(count.alternate_reads), QString::number(count.informative_read_depth()), percent(count.allele_fraction()),
            evidence->molecule_counts_available ? QString::number(count.alternate_molecules) : "N/A",
            evidence->molecule_counts_available ? QString::number(count.molecule_depth()) : "N/A",
            evidence->molecule_counts_available ? percent(count.molecule_allele_fraction()) : "N/A",
            QString::number(count.other_reads),
            evidence->molecule_counts_available ? QString::number(count.other_molecules) : "N/A",
            QString("%1 / %2").arg(count.alternate_forward_reads).arg(count.alternate_reverse_reads)};
        for (int column = 0; column < values.size(); ++column) {
            auto* item = new QTableWidgetItem(values[column]);
            item->setTextAlignment(column >= 2 ? Qt::AlignCenter : Qt::AlignVCenter | Qt::AlignLeft);
            if (column == 0) item->setToolTip(sanitized_resource_uri(result_bam_paths_[row]));
            if ((column == 4 || column == 7) && count.alternate_reads > 0) {
                item->setForeground(QColor(dark_mode_ ? "#b8adff" : "#5b48d6"));
            }
            results_->setItem(row, column, item);
        }
    }
    run_button_->setEnabled(true);
    load_bams_button_->setEnabled(true);
    refresh_loaded_bams();
    status_->setText(QString("%1 result(s) ready · %2 issue(s)")
        .arg(last_batch_.results.size()).arg(last_batch_.errors.size()));
    QString issue_text;
    for (const auto& error : last_batch_.errors) issue_text += QString::fromStdString(error) + '\n';
    status_->setToolTip(issue_text.trimmed());
    if (results_->rowCount() > 0) {
        results_->selectRow(0);
        show_result_summary(0);
    } else {
        summary_text_->setPlainText(issue_text.trimmed());
        copy_summary_button_->setEnabled(!issue_text.trimmed().isEmpty());
    }
}

void MainWindow::show_result_summary(const int row, const int) {
    if (row < 0 || static_cast<std::size_t>(row) >= last_batch_.results.size()) return;
    const auto* evidence = std::get_if<VariantEvidence>(&last_batch_.results[static_cast<std::size_t>(row)]);
    if (evidence == nullptr) return;
    summary_text_->setPlainText(evidence_summary(*evidence, resource_label(result_bam_paths_[row])));
    copy_summary_button_->setEnabled(true);
}

void MainWindow::copy_summary() {
    const auto summary = summary_text_->toPlainText().trimmed();
    if (summary.isEmpty()) return;
    QApplication::clipboard()->setText(summary);
    status_->setText("Summary copied to clipboard.");
}

}  // namespace bamseek
