#include <bamseek/main_window.hpp>

#include <bamseek/igv_command_receiver.hpp>
#include <bamseek/query.hpp>

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QDateTime>
#include <QDragEnterEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontDatabase>
#include <QFontMetrics>
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
#include <QScrollArea>
#include <QSettings>
#include <QSize>
#include <QSplitter>
#include <QTabWidget>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QTabBar>
#include <QTimer>
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
    if (!query.variant_type.empty()) clinical += ' ' + QString::fromStdString(query.variant_type);
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

QString variant_origin_label(const VariantOrigin origin) {
    if (origin == VariantOrigin::current) return "Current";
    if (origin == VariantOrigin::historical) return "Historical";
    return "Both";
}

QString compact_resource_label(const QString& path, const QFont& font) {
    return QFontMetrics(font).elidedText(resource_label(path), Qt::ElideMiddle, 176);
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

QLabel* role_label(const QString& text, const char* object_name, QWidget* parent) {
    auto* label = new QLabel(text, parent);
    label->setObjectName(object_name);
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
QFrame#SummaryPanel {
    background: @SUMMARY@;
    border: 1px solid @SUMMARY_BORDER@;
    border-radius: 10px;
}
QScrollArea { border: 0; background: @BG@; }
QScrollArea > QWidget > QWidget { background: @BG@; }
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
QLabel#CurrentRole, QLabel#HistoricalRole {
    border-radius: 8px;
    padding: 3px 9px;
    font-weight: 650;
}
QLabel#CurrentRole { color: @CURRENT_TEXT@; background: @CURRENT_BG@; border: 1px solid @CURRENT_BORDER@; }
QLabel#HistoricalRole { color: @HISTORY_TEXT@; background: @HISTORY_BG@; border: 1px solid @HISTORY_BORDER@; }
QPlainTextEdit, QLineEdit, QComboBox, QListWidget, QTableWidget {
    background: @FIELD@;
    color: @TEXT@;
    border: 1px solid @FIELD_BORDER@;
    border-radius: 9px;
    selection-background-color: @SELECT@;
    selection-color: @SELECT_TEXT@;
}
QPlainTextEdit, QLineEdit, QComboBox, QListWidget { padding: 8px; }
QPlainTextEdit:focus, QLineEdit:focus, QComboBox:focus, QListWidget:focus, QTableWidget:focus {
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
    replace("@CURRENT_TEXT@", "#c8c1ff", "#5541c8", dark);
    replace("@CURRENT_BG@", "#28234f", "#efedff", dark);
    replace("@CURRENT_BORDER@", "#4a4084", "#d5ceff", dark);
    replace("@HISTORY_TEXT@", "#8fddd2", "#087568", dark);
    replace("@HISTORY_BG@", "#102f31", "#e5f7f4", dark);
    replace("@HISTORY_BORDER@", "#245254", "#b9e4de", dark);
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
    resize(1440, 980);
    setMinimumSize(1080, 760);
    setFont(QFontDatabase::systemFont(QFontDatabase::GeneralFont));
    QSettings stored_settings;
    dark_mode_ = stored_settings.value("appearance/dark_mode", true).toBool();
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

    auto* analysis_page = new QScrollArea(tabs_);
    analysis_page->setWidgetResizable(true);
    analysis_page->setFrameShape(QFrame::NoFrame);
    analysis_page->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    analysis_page->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    auto* analysis_content = new QWidget(analysis_page);
    analysis_content->setObjectName("AppRoot");
    analysis_content->setMinimumSize(1020, 680);
    analysis_page->setWidget(analysis_content);
    auto* analysis_layout = new QVBoxLayout(analysis_content);
    analysis_layout->setContentsMargins(4, 14, 8, 8);
    analysis_layout->setSpacing(14);

    auto* input_splitter = new QSplitter(Qt::Horizontal, analysis_content);
    input_splitter->setChildrenCollapsible(false);
    input_splitter->setFixedHeight(340);

    auto* bam_card = new QFrame(input_splitter);
    bam_card->setObjectName("Card");
    auto* bam_layout = new QVBoxLayout(bam_card);
    bam_layout->setContentsMargins(14, 12, 14, 14);
    bam_layout->setSpacing(8);
    auto* bam_header = new QHBoxLayout();
    bam_header->addWidget(section_title("BAM sources", bam_card));
    bam_header->addStretch(1);
    bam_count_ = new QLabel("0 ready", bam_card);
    bam_count_->setObjectName("Pill");
    bam_header->addWidget(bam_count_);
    bam_layout->addLayout(bam_header);
    bam_path_ = new QPlainTextEdit(bam_card);
    bam_path_->setPlaceholderText("/data/sample.bam\nhttps://example.org/sample.bam");
    bam_path_->setMaximumHeight(58);
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
    loaded_bams_->setMinimumHeight(70);
    loaded_bams_->setToolTip("Check one BAM to designate the current sample. Unchecked BAMs are historical.");
    bam_layout->addWidget(loaded_bams_, 1);

    auto* variant_card = new QFrame(input_splitter);
    variant_card->setObjectName("Card");
    auto* variant_layout = new QVBoxLayout(variant_card);
    variant_layout->setContentsMargins(14, 12, 14, 14);
    variant_layout->setSpacing(8);
    variant_layout->addWidget(section_title("Variants", variant_card));
    auto* current_header = new QHBoxLayout();
    current_header->addWidget(role_label("Current variants", "CurrentRole", variant_card));
    current_header->addStretch(1);
    variant_layout->addLayout(current_header);
    current_query_text_ = new QPlainTextEdit(variant_card);
    current_query_text_->setPlaceholderText(
        "IGV 9 139390861 NOTCH1 stop_gained c.7330C>T p.Q2444* G A 3.1");
    variant_layout->addWidget(current_query_text_, 1);
    auto* historical_header = new QHBoxLayout();
    historical_header->addWidget(role_label("Historical variants", "HistoricalRole", variant_card));
    historical_header->addStretch(1);
    variant_layout->addLayout(historical_header);
    historical_query_text_ = new QPlainTextEdit(variant_card);
    historical_query_text_->setPlaceholderText(
        "5 176943930 DDX41 ENST00000507955.1 MISSENSE c.17C>T p.P6L G A");
    variant_layout->addWidget(historical_query_text_, 1);

    input_splitter->addWidget(bam_card);
    input_splitter->addWidget(variant_card);
    input_splitter->setSizes({580, 700});
    analysis_layout->addWidget(input_splitter, 1);

    auto* action_row = new QHBoxLayout();
    status_ = new QLabel("Add BAMs and variants to begin.", analysis_content);
    status_->setObjectName("MutedLabel");
    run_button_ = new QPushButton("Calculate VAFs", analysis_content);
    run_button_->setObjectName("PrimaryButton");
    run_button_->setEnabled(false);
    action_row->addWidget(status_, 1);
    action_row->addWidget(run_button_);
    analysis_layout->addLayout(action_row);

    auto* results_card = new QFrame(analysis_content);
    results_card->setObjectName("Card");
    auto* results_layout = new QVBoxLayout(results_card);
    results_layout->setContentsMargins(14, 12, 14, 14);
    results_layout->setSpacing(8);
    results_layout->addWidget(section_title("VAF results", results_card));
    results_ = new QTableWidget(results_card);
    results_->setColumnCount(13);
    results_->setHorizontalHeaderLabels({"Sample", "BAM", "Variant set", "Variant", "ALT reads", "REF reads", "VAF",
        "ALT mol.", "REF mol.", "Mol. VAF", "Other/N", "Ambig.", "F / R"});
    const QStringList result_header_tips{"Current or historical sample", "BAM source", "Variant input set", "Variant",
        "Alternate-supporting reads", "Reference-supporting reads", "Read variant allele frequency",
        "Alternate-supporting molecules", "Reference-supporting molecules", "Molecule variant allele frequency",
        "OTHER/N reads", "Ambiguous molecules", "ALT forward / reverse reads"};
    for (int column = 0; column < result_header_tips.size(); ++column) {
        results_->horizontalHeaderItem(column)->setToolTip(result_header_tips[column]);
    }
    results_->setAlternatingRowColors(true);
    results_->setSelectionBehavior(QAbstractItemView::SelectRows);
    results_->setSelectionMode(QAbstractItemView::SingleSelection);
    results_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    results_->setShowGrid(false);
    results_->verticalHeader()->setVisible(false);
    results_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    results_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Fixed);
    results_->horizontalHeader()->resizeSection(1, 190);
    results_->horizontalHeader()->setSectionResizeMode(3, QHeaderView::Stretch);
    results_->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    results_->setMinimumHeight(100);
    results_layout->addWidget(results_, 1);
    auto* summary_panel = new QFrame(results_card);
    summary_panel->setObjectName("SummaryPanel");
    auto* summary_layout = new QVBoxLayout(summary_panel);
    summary_layout->setContentsMargins(10, 6, 10, 8);
    summary_layout->setSpacing(5);
    auto* summary_header = new QHBoxLayout();
    summary_header->addWidget(muted_label("Comparison narrative", summary_panel));
    summary_header->addStretch(1);
    copy_summary_button_ = new QPushButton("Copy summary", summary_panel);
    copy_summary_button_->setEnabled(false);
    summary_header->addWidget(copy_summary_button_);
    summary_layout->addLayout(summary_header);
    summary_text_ = new QPlainTextEdit(summary_panel);
    summary_text_->setObjectName("SummaryBox");
    summary_text_->setReadOnly(true);
    summary_text_->setMinimumHeight(120);
    summary_text_->setMaximumHeight(190);
    summary_text_->setPlaceholderText("Run the analysis to generate copy-ready comparison paragraphs.");
    summary_layout->addWidget(summary_text_);
    results_layout->addWidget(summary_panel);
    analysis_layout->addWidget(results_card, 1);

    auto* settings_page = new QScrollArea(tabs_);
    settings_page->setWidgetResizable(true);
    settings_page->setFrameShape(QFrame::NoFrame);
    auto* settings_content = new QWidget(settings_page);
    settings_content->setObjectName("AppRoot");
    settings_page->setWidget(settings_content);
    auto* settings_layout = new QVBoxLayout(settings_content);
    settings_layout->setContentsMargins(4, 14, 8, 8);
    settings_layout->setSpacing(14);

    auto* evidence_settings_card = new QFrame(settings_content);
    evidence_settings_card->setObjectName("Card");
    auto* evidence_settings_layout = new QVBoxLayout(evidence_settings_card);
    evidence_settings_layout->setContentsMargins(18, 16, 18, 18);
    evidence_settings_layout->setSpacing(12);
    evidence_settings_layout->addWidget(section_title("Evidence calculation", evidence_settings_card));
    evidence_settings_layout->addWidget(muted_label(
        "These settings apply to the next VAF calculation. Changes are saved automatically.", evidence_settings_card));

    auto* quality_controls = new QHBoxLayout();
    quality_controls->addWidget(muted_label("Minimum mapping quality", evidence_settings_card));
    minimum_mapq_ = new QLineEdit(
        QString::number(stored_settings.value("analysis/minimum_mapping_quality", 20).toInt()), evidence_settings_card);
    minimum_mapq_->setMaximumWidth(72);
    minimum_mapq_->setAlignment(Qt::AlignCenter);
    minimum_mapq_->setValidator(new QIntValidator(0, 255, minimum_mapq_));
    quality_controls->addWidget(minimum_mapq_);
    quality_controls->addSpacing(22);
    quality_controls->addWidget(muted_label("Minimum base quality", evidence_settings_card));
    minimum_baseq_ = new QLineEdit(
        QString::number(stored_settings.value("analysis/minimum_base_quality", 20).toInt()), evidence_settings_card);
    minimum_baseq_->setMaximumWidth(72);
    minimum_baseq_->setAlignment(Qt::AlignCenter);
    minimum_baseq_->setValidator(new QIntValidator(0, 255, minimum_baseq_));
    quality_controls->addWidget(minimum_baseq_);
    quality_controls->addStretch(1);
    evidence_settings_layout->addLayout(quality_controls);

    auto* alignment_flags = new QHBoxLayout();
    alignment_flags->addWidget(muted_label("Include alignments marked as", evidence_settings_card));
    include_duplicates_ = new QCheckBox("Duplicates", evidence_settings_card);
    include_duplicates_->setChecked(stored_settings.value("analysis/include_duplicates", false).toBool());
    alignment_flags->addWidget(include_duplicates_);
    include_secondary_ = new QCheckBox("Secondary", evidence_settings_card);
    include_secondary_->setChecked(stored_settings.value("analysis/include_secondary", false).toBool());
    alignment_flags->addWidget(include_secondary_);
    include_supplementary_ = new QCheckBox("Supplementary", evidence_settings_card);
    include_supplementary_->setChecked(stored_settings.value("analysis/include_supplementary", false).toBool());
    alignment_flags->addWidget(include_supplementary_);
    alignment_flags->addStretch(1);
    evidence_settings_layout->addLayout(alignment_flags);

    auto* molecule_controls = new QHBoxLayout();
    molecule_controls->addWidget(muted_label("Molecule grouping", evidence_settings_card));
    molecule_mode_ = new QComboBox(evidence_settings_card);
    molecule_mode_->addItem("Paired fragments (read name)", static_cast<int>(MoleculeMode::raw_reads));
    molecule_mode_->addItem("Auto-detect MI / RX / UB", static_cast<int>(MoleculeMode::auto_detect));
    molecule_mode_->addItem("Specific BAM tag", static_cast<int>(MoleculeMode::selected_tag));
    const auto stored_molecule_mode = stored_settings.value(
        "analysis/molecule_mode", static_cast<int>(MoleculeMode::raw_reads)).toInt();
    const auto molecule_mode_index = molecule_mode_->findData(stored_molecule_mode);
    molecule_mode_->setCurrentIndex(molecule_mode_index < 0 ? 0 : molecule_mode_index);
    molecule_controls->addWidget(molecule_mode_);
    molecule_controls->addSpacing(12);
    molecule_controls->addWidget(muted_label("Tag", evidence_settings_card));
    molecule_tag_ = new QLineEdit(stored_settings.value("analysis/molecule_tag", "MI").toString().toUpper(),
                                  evidence_settings_card);
    molecule_tag_->setPlaceholderText("MI");
    molecule_tag_->setMaxLength(2);
    molecule_tag_->setMaximumWidth(64);
    molecule_tag_->setAlignment(Qt::AlignCenter);
    molecule_tag_->setEnabled(molecule_mode_->currentData().toInt() == static_cast<int>(MoleculeMode::selected_tag));
    molecule_controls->addWidget(molecule_tag_);
    molecule_controls->addStretch(1);
    evidence_settings_layout->addLayout(molecule_controls);
    evidence_settings_layout->addWidget(muted_label(
        "Auto-detect uses MI, RX, or UB when at least 90% of callable alignments carry the tag; otherwise reads are grouped into paired fragments.",
        evidence_settings_card));
    settings_layout->addWidget(evidence_settings_card);

    auto* receiver_card = new QFrame(settings_content);
    receiver_card->setObjectName("Card");
    auto* receiver_layout = new QVBoxLayout(receiver_card);
    receiver_layout->setContentsMargins(18, 16, 18, 18);
    receiver_layout->setSpacing(10);
    receiver_layout->addWidget(section_title("IGV command relay", receiver_card));
    receiver_layout->addWidget(muted_label(
        "BAM Seek listens only on localhost, records broadcast commands below, and forwards their bytes unchanged to IGV.",
        receiver_card));
    auto* receiver_controls = new QHBoxLayout();
    receiver_enabled_ = new QCheckBox("Receiver enabled", receiver_card);
    receiver_enabled_->setChecked(stored_settings.value("receiver/enabled", true).toBool());
    receiver_port_ = new QLineEdit(
        QString::number(stored_settings.value("receiver/listen_port", 60151).toInt()), receiver_card);
    receiver_port_->setValidator(new QIntValidator(1, 65535, receiver_port_));
    receiver_port_->setAlignment(Qt::AlignCenter);
    receiver_port_->setMaximumWidth(82);
    upstream_port_ = new QLineEdit(
        QString::number(stored_settings.value("receiver/igv_port", 60152).toInt()), receiver_card);
    upstream_port_->setValidator(new QIntValidator(1, 65535, upstream_port_));
    upstream_port_->setAlignment(Qt::AlignCenter);
    upstream_port_->setMaximumWidth(82);
    auto* apply_receiver_ports = new QPushButton("Apply ports", receiver_card);
    receiver_status_ = new QLabel(receiver_card);
    receiver_status_->setObjectName("MutedLabel");
    receiver_controls->addWidget(receiver_enabled_);
    receiver_controls->addSpacing(16);
    receiver_controls->addWidget(muted_label("Receive on 127.0.0.1", receiver_card));
    receiver_controls->addWidget(receiver_port_);
    receiver_controls->addSpacing(16);
    receiver_controls->addWidget(muted_label("Send to IGV on 127.0.0.1", receiver_card));
    receiver_controls->addWidget(upstream_port_);
    receiver_controls->addSpacing(16);
    receiver_controls->addWidget(apply_receiver_ports);
    receiver_controls->addStretch(1);
    receiver_layout->addLayout(receiver_controls);
    receiver_layout->addWidget(receiver_status_);
    receiver_layout->addWidget(muted_label("Broadcast activity", receiver_card));
    broadcast_text_ = new QPlainTextEdit(receiver_card);
    broadcast_text_->setReadOnly(true);
    broadcast_text_->setPlaceholderText("Received /load and /goto commands will appear here…");
    broadcast_text_->setMinimumHeight(240);
    receiver_layout->addWidget(broadcast_text_, 1);
    settings_layout->addWidget(receiver_card, 1);

    tabs_->addTab(analysis_page, "VAF workspace");
    tabs_->addTab(settings_page, "Settings");
    root_layout->addWidget(tabs_);
    setCentralWidget(root);

    connect(browse_bam, &QPushButton::clicked, this, [this] { choose_bams(); });
    connect(load_bams_button_, &QPushButton::clicked, this, [this] { load_bams(); });
    connect(remove_bams_button_, &QPushButton::clicked, this, [this] { remove_selected_bams(); });
    connect(clear_bams_button_, &QPushButton::clicked, this, [this] { clear_loaded_bams(); });
    connect(loaded_bams_, &QListWidget::itemSelectionChanged, this, [this] {
        remove_bams_button_->setEnabled(!busy() && !loaded_bams_->selectedItems().isEmpty());
    });
    connect(loaded_bams_, &QListWidget::itemChanged, this, [this](QListWidgetItem* item) { set_current_bam(item); });
    connect(run_button_, &QPushButton::clicked, this, [this] { run_queries(); });
    connect(copy_summary_button_, &QPushButton::clicked, this, [this] { copy_summary(); });
    connect(theme_button_, &QPushButton::clicked, this, [this] { toggle_theme(); });
    connect(&bam_load_watcher_, &QFutureWatcher<BamLoadBatch>::finished, this, [this] { bams_loaded(); });
    connect(&watcher_, &QFutureWatcher<MultiBamBatch>::finished, this, [this] { show_results(); });

    const auto save_analysis = [this] { save_analysis_settings(); };
    connect(minimum_mapq_, &QLineEdit::editingFinished, this, save_analysis);
    connect(minimum_baseq_, &QLineEdit::editingFinished, this, save_analysis);
    connect(include_duplicates_, &QCheckBox::toggled, this, save_analysis);
    connect(include_secondary_, &QCheckBox::toggled, this, save_analysis);
    connect(include_supplementary_, &QCheckBox::toggled, this, save_analysis);
    connect(molecule_mode_, &QComboBox::currentIndexChanged, this, [this](const int) {
        molecule_tag_->setEnabled(molecule_mode_->currentData().toInt() == static_cast<int>(MoleculeMode::selected_tag));
        save_analysis_settings();
    });
    connect(molecule_tag_, &QLineEdit::editingFinished, this, [this] {
        molecule_tag_->setText(molecule_tag_->text().trimmed().toUpper());
        save_analysis_settings();
    });

    command_receiver_ = new IgvCommandReceiver(this);
    connect(receiver_enabled_, &QCheckBox::toggled, this, [this](const bool enabled) { set_receiver_enabled(enabled); });
    connect(apply_receiver_ports, &QPushButton::clicked, this,
        [this] { set_receiver_enabled(receiver_enabled_->isChecked()); });
    connect(command_receiver_, &IgvCommandReceiver::request_received, this, [this](const QString& description) {
        append_received_command(description);
    });
    connect(command_receiver_, &IgvCommandReceiver::bam_load_requested, this,
        [this](const QStringList& paths) { enqueue_received_bams(paths); });
    connect(command_receiver_, &IgvCommandReceiver::listening_changed, this,
        [this](const bool listening, const quint16 port, const QString& error) {
            receiver_status_->setText(listening
                ? QString("Listening on 127.0.0.1:%1 → IGV 127.0.0.1:%2").arg(port).arg(upstream_port_->text())
                                                : (error.isEmpty() ? "Receiver stopped" : "Unavailable · " + error));
        });
    connect(command_receiver_, &IgvCommandReceiver::forwarding_changed, this,
        [this](const bool connected, const quint16 port, const QString& error) {
            if (connected) {
                receiver_status_->setText(QString("Listening on 127.0.0.1:%1 → IGV 127.0.0.1:%2")
                    .arg(command_receiver_->port()).arg(port));
            } else {
                receiver_status_->setText(QString("Listening on 127.0.0.1:%1 · IGV 127.0.0.1:%2 unavailable")
                    .arg(command_receiver_->port()).arg(port));
                append_received_command(QString("IGV forwarding failed on 127.0.0.1:%1\n%2").arg(port).arg(error));
            }
        });
    set_receiver_enabled(receiver_enabled_->isChecked());
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
        const QColor current_color(dark_mode_ ? "#c8c1ff" : "#5541c8");
        const QColor historical_color(dark_mode_ ? "#8fddd2" : "#087568");
        for (int row = 0; row < loaded_bams_->count(); ++row) {
            auto* item = loaded_bams_->item(row);
            item->setForeground(item->checkState() == Qt::Checked ? current_color : historical_color);
        }
    }
    if (results_ != nullptr) {
        const QColor accent(dark_mode_ ? "#b8adff" : "#5b48d6");
        const QColor current_color(dark_mode_ ? "#c8c1ff" : "#5541c8");
        const QColor historical_color(dark_mode_ ? "#8fddd2" : "#087568");
        for (int row = 0; row < results_->rowCount(); ++row) {
            if (results_->item(row, 0) != nullptr) {
                results_->item(row, 0)->setForeground(results_->item(row, 0)->data(Qt::UserRole).toBool()
                    ? current_color : historical_color);
            }
            if (results_->item(row, 2) != nullptr) {
                const auto origin = static_cast<VariantOrigin>(results_->item(row, 2)->data(Qt::UserRole).toInt());
                results_->item(row, 2)->setForeground(origin == VariantOrigin::current ? current_color
                    : origin == VariantOrigin::historical ? historical_color : accent);
            }
            if (results_->item(row, 6) != nullptr) results_->item(row, 6)->setForeground(accent);
            if (results_->item(row, 9) != nullptr) results_->item(row, 9)->setForeground(accent);
        }
    }
}

void MainWindow::set_receiver_enabled(const bool enabled) {
    QSettings settings;
    settings.setValue("receiver/enabled", enabled);
    if (!enabled) {
        if (receiver_port_->hasAcceptableInput() && upstream_port_->hasAcceptableInput()
            && receiver_port_->text().toUShort() != upstream_port_->text().toUShort()) {
            settings.setValue("receiver/listen_port", receiver_port_->text().toUShort());
            settings.setValue("receiver/igv_port", upstream_port_->text().toUShort());
        }
        command_receiver_->close();
        receiver_status_->setText("Receiver stopped");
        return;
    }
    if (!receiver_port_->hasAcceptableInput() || !upstream_port_->hasAcceptableInput()) {
        receiver_status_->setText("Port changes not applied: enter receiving and IGV ports from 1 to 65535.");
        return;
    }
    const auto receive_port = receiver_port_->text().toUShort();
    const auto send_port = upstream_port_->text().toUShort();
    if (receive_port == send_port) {
        receiver_status_->setText("Port changes not applied: receiving and IGV ports must be different.");
        return;
    }
    settings.setValue("receiver/listen_port", receive_port);
    settings.setValue("receiver/igv_port", send_port);
    (void)command_receiver_->listen(receive_port, send_port);
}

void MainWindow::save_analysis_settings() {
    QSettings settings;
    if (minimum_mapq_->hasAcceptableInput()) {
        settings.setValue("analysis/minimum_mapping_quality", minimum_mapq_->text().toInt());
    }
    if (minimum_baseq_->hasAcceptableInput()) {
        settings.setValue("analysis/minimum_base_quality", minimum_baseq_->text().toInt());
    }
    settings.setValue("analysis/include_duplicates", include_duplicates_->isChecked());
    settings.setValue("analysis/include_secondary", include_secondary_->isChecked());
    settings.setValue("analysis/include_supplementary", include_supplementary_->isChecked());
    settings.setValue("analysis/molecule_mode", molecule_mode_->currentData().toInt());
    settings.setValue("analysis/molecule_tag", molecule_tag_->text().trimmed().toUpper());
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
    start_bam_preparation(pending, true);
}

void MainWindow::start_bam_preparation(const QStringList& pending, const bool clear_manual_input) {
    if (busy() || pending.isEmpty()) return;

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
    loaded_bams_->setEnabled(false);
    clear_bam_input_after_load_ = clear_manual_input;
    preparing_bam_paths_ = candidates;
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
    preparing_bam_paths_.clear();
    for (qsizetype i = 0; i < loaded.bams.size(); ++i) {
        loaded_bam_paths_.append(loaded.bams[i]);
        loaded_index_paths_.append(loaded.indexes[i]);
        loaded_engines_.push_back(loaded.engines[static_cast<std::size_t>(i)]);
    }
    if (!loaded.bams.isEmpty()) {
        if (clear_bam_input_after_load_) bam_path_->clear();
        clear_results();
    }
    bam_path_->setEnabled(true);
    load_bams_button_->setEnabled(true);
    run_button_->setEnabled(true);
    loaded_bams_->setEnabled(true);
    refresh_loaded_bams();
    status_->setText(QString("%1 BAM(s) ready. No VAF calculation has run yet.").arg(loaded_bam_paths_.size()));
    if (!loaded.issues.isEmpty()) QMessageBox::warning(this, "BAM preparation notes", loaded.issues.join('\n'));
    clear_bam_input_after_load_ = false;
    if (!queued_receiver_bams_.isEmpty()) {
        QTimer::singleShot(0, this, [this] { load_queued_receiver_bams(); });
    }
}

void MainWindow::enqueue_received_bams(const QStringList& paths) {
    for (auto path : paths) {
        path = path.trimmed();
        if (path.isEmpty() || loaded_bam_paths_.contains(path) || preparing_bam_paths_.contains(path)
            || queued_receiver_bams_.contains(path)) continue;
        queued_receiver_bams_.append(path);
    }
    if (queued_receiver_bams_.isEmpty()) return;
    if (busy()) {
        status_->setText(QString("%1 receiver BAM(s) queued for preparation.").arg(queued_receiver_bams_.size()));
        return;
    }
    load_queued_receiver_bams();
}

void MainWindow::load_queued_receiver_bams() {
    if (busy() || queued_receiver_bams_.isEmpty()) return;
    const auto pending = queued_receiver_bams_;
    queued_receiver_bams_.clear();
    start_bam_preparation(pending, false);
}

void MainWindow::remove_selected_bams() {
    if (busy()) return;
    QList<int> rows;
    for (auto* item : loaded_bams_->selectedItems()) rows.append(loaded_bams_->row(item));
    std::sort(rows.begin(), rows.end(), std::greater<>());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    for (const int row : rows) {
        evidence_cache_.erase(loaded_bam_paths_[row].toStdString());
        if (loaded_bam_paths_[row] == current_bam_path_) current_bam_path_.clear();
        loaded_bam_paths_.removeAt(row);
        loaded_index_paths_.removeAt(row);
        loaded_engines_.erase(loaded_engines_.begin() + row);
    }
    clear_results();
    refresh_loaded_bams();
    status_->setText(QString("%1 BAM(s) ready. Results cleared after source changes.").arg(loaded_bam_paths_.size()));
}

void MainWindow::clear_loaded_bams() {
    if (busy()) return;
    loaded_bam_paths_.clear();
    loaded_index_paths_.clear();
    loaded_engines_.clear();
    evidence_cache_.clear();
    current_bam_path_.clear();
    clear_results();
    refresh_loaded_bams();
    status_->setText("BAM panel cleared.");
}

void MainWindow::refresh_loaded_bams() {
    if (!current_bam_path_.isEmpty() && !loaded_bam_paths_.contains(current_bam_path_)) current_bam_path_.clear();
    refreshing_bam_list_ = true;
    loaded_bams_->clear();
    for (qsizetype i = 0; i < loaded_bam_paths_.size(); ++i) {
        const auto& bam = loaded_bam_paths_[i];
        const bool remote = bam.startsWith("https://");
        const bool current = bam == current_bam_path_;
        auto* item = new QListWidgetItem(QString("%1\n%2 · %3 · Ready")
            .arg(resource_label(bam), current ? "Current" : "Historical", remote ? "HTTPS" : "Local"), loaded_bams_);
        item->setData(Qt::UserRole, bam);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(current ? Qt::Checked : Qt::Unchecked);
        item->setForeground(QColor(current ? (dark_mode_ ? "#c8c1ff" : "#5541c8")
                                           : (dark_mode_ ? "#8fddd2" : "#087568")));
        item->setSizeHint(QSize(0, 48));
        item->setToolTip(QString("%1 sample · check to make this the only current BAM\nBAM: ")
            .arg(current ? "Current" : "Historical") + sanitized_resource_uri(bam) + "\nIndex: "
            + (loaded_index_paths_[i].isEmpty() ? "automatic remote discovery" : loaded_index_paths_[i]));
    }
    refreshing_bam_list_ = false;
    bam_count_->setText(QString("%1 ready").arg(loaded_bam_paths_.size()));
    const bool has_bams = !loaded_bam_paths_.isEmpty();
    run_button_->setEnabled(has_bams && !busy());
    remove_bams_button_->setEnabled(has_bams && !loaded_bams_->selectedItems().isEmpty() && !busy());
    clear_bams_button_->setEnabled(has_bams && !busy());
}

void MainWindow::set_current_bam(QListWidgetItem* changed_item) {
    if (refreshing_bam_list_ || changed_item == nullptr || busy()) return;
    const auto path = changed_item->data(Qt::UserRole).toString();
    const auto previous = current_bam_path_;
    if (changed_item->checkState() == Qt::Checked) current_bam_path_ = path;
    else if (current_bam_path_ == path) current_bam_path_.clear();
    if (current_bam_path_ == previous) return;
    refreshing_bam_list_ = true;
    for (int row = 0; row < loaded_bams_->count(); ++row) {
        auto* item = loaded_bams_->item(row);
        const auto item_path = item->data(Qt::UserRole).toString();
        const bool current = item_path == current_bam_path_;
        const bool remote = item_path.startsWith("https://");
        item->setCheckState(current ? Qt::Checked : Qt::Unchecked);
        item->setText(QString("%1\n%2 · %3 · Ready")
            .arg(resource_label(item_path), current ? "Current" : "Historical", remote ? "HTTPS" : "Local"));
    }
    refreshing_bam_list_ = false;
    apply_theme();
    const bool has_results = !last_batch_.results.empty();
    if (has_results) render_results();
    status_->setText(current_bam_path_.isEmpty()
        ? QString("No current BAM selected; all BAMs are historical.%1").arg(has_results ? " Existing VAF evidence reused." : "")
        : resource_label(current_bam_path_) + " designated as current."
            + (has_results ? " Existing VAF evidence reused." : ""));
}

void MainWindow::clear_results() {
    last_batch_ = {};
    result_bam_paths_.clear();
    result_variant_origins_.clear();
    results_->setRowCount(0);
    summary_text_->clear();
    copy_summary_button_->setEnabled(false);
}

FilterSettings MainWindow::filters() const {
    FilterSettings values;
    values.minimum_mapping_quality = std::max(0, minimum_mapq_->text().toInt());
    values.minimum_base_quality = std::max(0, minimum_baseq_->text().toInt());
    values.include_duplicates = include_duplicates_->isChecked();
    values.include_secondary = include_secondary_->isChecked();
    values.include_supplementary = include_supplementary_->isChecked();
    values.minimum_variant_allele_fraction = 0.0;
    values.minimum_alternate_reads = 1;
    values.minimum_alternate_molecules = 1;
    values.molecule_mode = static_cast<MoleculeMode>(molecule_mode_->currentData().toInt());
    values.molecule_tag = molecule_tag_->text().trimmed().toUpper().toStdString();
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
    const auto selected_tag_mode = molecule_mode_->currentData().toInt() == static_cast<int>(MoleculeMode::selected_tag);
    const auto molecule_tag = molecule_tag_->text().trimmed();
    if (selected_tag_mode && (molecule_tag.size() != 2 || !molecule_tag.front().isLetter()
        || !molecule_tag.back().isLetterOrNumber())) {
        QMessageBox::warning(this, "Invalid molecule tag",
            "A specific BAM molecule tag must contain exactly two characters: a letter followed by a letter or number.");
        tabs_->setCurrentIndex(1);
        molecule_tag_->setFocus();
        return;
    }
    save_analysis_settings();
    auto current_parsed = parse_queries(current_query_text_->toPlainText().toStdString());
    auto historical_parsed = parse_queries(historical_query_text_->toPlainText().toStdString());
    std::vector<VariantQuery> current_variants;
    std::vector<VariantQuery> historical_variants;
    std::vector<std::string> errors;
    const auto collect = [&errors](ParsedQueries& parsed, std::vector<VariantQuery>& destination, const std::string& label) {
        for (auto& error : parsed.errors) errors.push_back('[' + label + "] " + std::move(error));
        for (auto& query : parsed.queries) {
            if (auto* variant = std::get_if<VariantQuery>(&query)) destination.push_back(std::move(*variant));
            else errors.push_back('[' + label + "] " + std::get<RegionQuery>(query).source_text
                + ": enter a specific REF>ALT variant rather than a region");
        }
    };
    collect(current_parsed, current_variants, "Current variants");
    collect(historical_parsed, historical_variants, "Historical variants");
    auto classified = classify_variants(current_variants, historical_variants);
    if (classified.empty()) {
        QMessageBox::warning(this, "No valid variants", QString::fromStdString(errors.empty()
            ? "Enter at least one current or historical variant." : errors.front()));
        return;
    }
    last_filters_ = filters();
    const auto bams = loaded_bam_paths_;
    const auto engines = loaded_engines_;
    auto cache_snapshot = evidence_cache_;
    run_button_->setEnabled(false);
    load_bams_button_->setEnabled(false);
    remove_bams_button_->setEnabled(false);
    clear_bams_button_->setEnabled(false);
    loaded_bams_->setEnabled(false);
    current_query_text_->setEnabled(false);
    historical_query_text_->setEnabled(false);
    status_->setText(QString("Calculating %1 distinct variant(s) across %2 prepared BAM(s)…")
        .arg(classified.size()).arg(bams.size()));
    watcher_.setFuture(QtConcurrent::run([classified = std::move(classified), errors = std::move(errors),
                                          filter_values = last_filters_, bams, engines,
                                          cache_snapshot = std::move(cache_snapshot)] {
        MultiBamBatch combined;
        combined.batch.errors = errors;
        for (std::size_t source_index = 0; source_index < engines.size(); ++source_index) {
            const auto& bam = bams[static_cast<qsizetype>(source_index)];
            const auto bam_key = bam.toStdString();
            const auto cached_bam = cache_snapshot.find(bam_key);
            for (const auto& item : classified) {
                const auto pair_key = evidence_cache_key(item.query, filter_values);
                if (cached_bam != cache_snapshot.end()) {
                    const auto cached = cached_bam->second.find(pair_key);
                    if (cached != cached_bam->second.end()) {
                        auto evidence = cached->second;
                        evidence.query = item.query;
                        combined.batch.results.emplace_back(std::move(evidence));
                        combined.result_bams.append(bam);
                        combined.result_variant_origins.push_back(item.origin);
                        ++combined.cache_hits;
                        continue;
                    }
                }
                ++combined.calculations;
                try {
                    auto evaluated = engines[source_index]->evaluate({Query{item.query}}, filter_values);
                    for (auto& error : evaluated.errors) {
                        combined.batch.errors.push_back("[" + resource_label(bam).toStdString() + "] " + error);
                    }
                    if (evaluated.results.empty()) continue;
                    auto* result = std::get_if<VariantEvidence>(&evaluated.results.front());
                    if (result == nullptr) continue;
                    result->reads.clear();
                    result->reads.shrink_to_fit();
                    combined.cache_updates[bam_key][pair_key] = *result;
                    combined.batch.results.emplace_back(std::move(*result));
                    combined.result_bams.append(bam);
                    combined.result_variant_origins.push_back(item.origin);
                } catch (const std::exception& error) {
                    combined.batch.errors.push_back("[" + resource_label(bam).toStdString() + "] "
                        + sanitized_error(error.what(), bam));
                }
            }
        }
        return combined;
    }));
}

void MainWindow::show_results() {
    auto combined = watcher_.result();
    for (auto& [bam, entries] : combined.cache_updates) {
        auto& destination = evidence_cache_[bam];
        for (auto& [key, evidence] : entries) destination[key] = std::move(evidence);
    }
    last_batch_ = std::move(combined.batch);
    result_bam_paths_ = std::move(combined.result_bams);
    result_variant_origins_ = std::move(combined.result_variant_origins);
    run_button_->setEnabled(true);
    load_bams_button_->setEnabled(true);
    loaded_bams_->setEnabled(true);
    current_query_text_->setEnabled(true);
    historical_query_text_->setEnabled(true);
    refresh_loaded_bams();
    render_results();
    status_->setText(QString("%1 result(s) ready · %2 reused · %3 calculated · %4 issue(s)")
        .arg(last_batch_.results.size()).arg(combined.cache_hits).arg(combined.calculations).arg(last_batch_.errors.size()));
    if (!queued_receiver_bams_.isEmpty()) {
        QTimer::singleShot(0, this, [this] { load_queued_receiver_bams(); });
    }
}

void MainWindow::render_results() {
    results_->setRowCount(0);
    std::vector<std::size_t> order(last_batch_.results.size());
    for (std::size_t index = 0; index < order.size(); ++index) order[index] = index;
    std::sort(order.begin(), order.end(), [this](const std::size_t left, const std::size_t right) {
        const bool left_current = !current_bam_path_.isEmpty() && result_bam_paths_[static_cast<qsizetype>(left)] == current_bam_path_;
        const bool right_current = !current_bam_path_.isEmpty() && result_bam_paths_[static_cast<qsizetype>(right)] == current_bam_path_;
        if (left_current != right_current) return left_current;
        const auto left_bam = resource_label(result_bam_paths_[static_cast<qsizetype>(left)]).toCaseFolded();
        const auto right_bam = resource_label(result_bam_paths_[static_cast<qsizetype>(right)]).toCaseFolded();
        if (left_bam != right_bam) return left_bam < right_bam;
        const auto* left_evidence = std::get_if<VariantEvidence>(&last_batch_.results[left]);
        const auto* right_evidence = std::get_if<VariantEvidence>(&last_batch_.results[right]);
        if (left_evidence == nullptr || right_evidence == nullptr) return left < right;
        return display_query(left_evidence->query).toCaseFolded() < display_query(right_evidence->query).toCaseFolded();
    });

    std::vector<ComparativeEvidence> narrative_evidence;
    narrative_evidence.reserve(last_batch_.results.size());
    for (const auto source_row : order) {
        const auto* evidence = std::get_if<VariantEvidence>(&last_batch_.results[source_row]);
        if (evidence == nullptr) continue;
        const auto bam_index = static_cast<qsizetype>(source_row);
        const auto& bam_path = result_bam_paths_[bam_index];
        const bool current_bam = !current_bam_path_.isEmpty() && bam_path == current_bam_path_;
        const auto origin = result_variant_origins_[source_row];
        const int row = results_->rowCount();
        results_->insertRow(row);
        const auto& count = evidence->counts;
        const QList<QString> values{current_bam ? "Current" : "Historical", compact_resource_label(bam_path, results_->font()),
            variant_origin_label(origin), display_query(evidence->query),
            QString::number(count.alternate_reads), QString::number(count.reference_reads), percent(count.allele_fraction()),
            evidence->molecule_counts_available ? QString::number(count.alternate_molecules) : "N/A",
            evidence->molecule_counts_available ? QString::number(count.reference_molecules) : "N/A",
            evidence->molecule_counts_available ? percent(count.molecule_allele_fraction()) : "N/A",
            QString::number(count.other_reads),
            evidence->molecule_counts_available ? QString::number(count.other_molecules) : "N/A",
            QString("%1 / %2").arg(count.alternate_forward_reads).arg(count.alternate_reverse_reads)};
        const auto diagnostic = QString("%1 overlapping alignment(s)\n%2 excluded by mapping quality or alignment flags\n"
                                        "%3 without a callable allele\n%4 below baseQ %5\n%6 counted as REF, ALT, or OTHER/N")
            .arg(evidence->overlapping_alignments)
            .arg(evidence->filtered_alignments)
            .arg(evidence->uncallable_alignments)
            .arg(evidence->low_base_quality_alignments)
            .arg(last_filters_.minimum_base_quality)
            .arg(count.depth());
        for (int column = 0; column < values.size(); ++column) {
            auto* item = new QTableWidgetItem(values[column]);
            item->setTextAlignment(column >= 4 ? Qt::AlignCenter : Qt::AlignVCenter | Qt::AlignLeft);
            if (column == 0) item->setData(Qt::UserRole, current_bam);
            if (column == 1) item->setToolTip(sanitized_resource_uri(bam_path));
            if (column == 2) item->setData(Qt::UserRole, static_cast<int>(origin));
            if (column == 3 || (column >= 4 && column <= 10)) item->setToolTip(diagnostic);
            results_->setItem(row, column, item);
        }
        narrative_evidence.push_back({resource_label(bam_path).toStdString(), current_bam, origin, *evidence});
    }
    apply_theme();
    QString issue_text;
    for (const auto& error : last_batch_.errors) issue_text += QString::fromStdString(error) + '\n';
    status_->setToolTip(issue_text.trimmed());
    const auto narrative = QString::fromStdString(comparison_narrative(narrative_evidence)).trimmed();
    summary_text_->setPlainText(narrative.isEmpty()
        ? "No comparison-specific narrative was generated. Variants listed in both sets are excluded, and historical-only variants require a current BAM."
        : narrative);
    copy_summary_button_->setEnabled(!narrative.isEmpty());
}

void MainWindow::copy_summary() {
    const auto summary = summary_text_->toPlainText().trimmed();
    if (summary.isEmpty()) return;
    QApplication::clipboard()->setText(summary);
    status_->setText("Summary copied to clipboard.");
}

}  // namespace bamseek
