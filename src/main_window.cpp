#include <bamseek/main_window.hpp>

#include <bamseek/igv_command_receiver.hpp>
#include <bamseek/pileup_view.hpp>
#include <bamseek/query.hpp>

#include <QComboBox>
#include <QCheckBox>
#include <QCryptographicHash>
#include <QDateTime>
#include <QDoubleValidator>
#include <QDragEnterEvent>
#include <QFileDialog>
#include <QFile>
#include <QFileInfo>
#include <QFormLayout>
#include <QFuture>
#include <QtConcurrentRun>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QIntValidator>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QMimeData>
#include <QIODevice>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QRegularExpression>
#include <QRegularExpressionValidator>
#include <QSaveFile>
#include <QScrollArea>
#include <QSpinBox>
#include <QSplitter>
#include <QTabWidget>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QVBoxLayout>
#include <QWidget>
#include <QUrl>

#include <algorithm>
#include <stdexcept>

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
    return clinical + " | " + genomic;
}

QString mode_name(const MoleculeMode mode) {
    switch (mode) {
        case MoleculeMode::raw_reads: return "Read pairs/fragments (no UMI)";
        case MoleculeMode::auto_detect: return "Auto-detect UMI; pair fallback";
        case MoleculeMode::selected_tag: return "Selected UMI tag; pair fallback";
    }
    return "Unknown";
}

QString evidence_summary(const VariantEvidence& evidence) {
    const auto& counts = evidence.counts;
    const auto query = display_query(evidence.query);
    const auto read_vaf = QString::number(counts.allele_fraction() * 100.0, 'f', 4) + '%';
    if (!evidence.molecule_counts_available) {
        return QString("For %1, ALT (%2) was seen in %3/%4 individual reads (read VAF %5); molecule consensus was unavailable. %6 OTHER/N read(s) were excluded from VAF.")
            .arg(query, QString::fromStdString(evidence.query.alternate))
            .arg(counts.alternate_reads).arg(counts.informative_read_depth()).arg(read_vaf).arg(counts.other_reads);
    }
    const auto molecule_vaf = QString::number(counts.molecule_allele_fraction() * 100.0, 'f', 4) + '%';
    return QString("For %1, ALT (%2) was seen in %3/%4 consensus reads (unique molecules; molecule VAF %5) and %6/%7 individual reads (read VAF %8). %9 OTHER/N read(s) and %10 ambiguous molecule(s) were excluded from VAF.")
        .arg(query, QString::fromStdString(evidence.query.alternate))
        .arg(counts.alternate_molecules).arg(counts.molecule_depth()).arg(molecule_vaf)
        .arg(counts.alternate_reads).arg(counts.informative_read_depth()).arg(read_vaf)
        .arg(counts.other_reads).arg(counts.other_molecules);
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

std::string sanitized_error(std::string message, const std::vector<std::string>& resources) {
    auto sanitized_message = QString::fromStdString(message);
    for (const auto& resource : resources) {
        const auto original = QString::fromStdString(resource);
        if (original.startsWith("https://")) sanitized_message.replace(original, sanitized_resource_uri(original));
    }
    return sanitized_message.toStdString();
}

QJsonObject file_identity(const QString& path) {
    if (path.startsWith("https://")) {
        return QJsonObject{{"uri", sanitized_resource_uri(path)}, {"remote", true}, {"sha256", QJsonValue::Null}};
    }
    QFileInfo info(path);
    QJsonObject identity{{"path", path}, {"exists", info.exists()}};
    if (!info.exists() || !info.isFile()) return identity;
    const auto original_size = info.size();
    const auto original_modified = info.lastModified();
    identity["bytes"] = static_cast<qint64>(original_size);
    identity["modified_utc"] = original_modified.toUTC().toString(Qt::ISODateWithMs);
    QFile input(path);
    if (!input.open(QIODevice::ReadOnly)) return identity;
    QCryptographicHash digest(QCryptographicHash::Sha256);
    while (!input.atEnd()) digest.addData(input.read(4 * 1024 * 1024));
    identity["sha256"] = QString::fromLatin1(digest.result().toHex());
    const QFileInfo after(path);
    identity["stable_during_hash"] = after.size() == original_size && after.lastModified() == original_modified;
    return identity;
}

QString index_path_for(const QString& alignment) {
    const QStringList candidates{alignment + ".bai", alignment + ".csi", alignment + ".crai",
        QFileInfo(alignment).absolutePath() + '/' + QFileInfo(alignment).completeBaseName() + ".bai",
        QFileInfo(alignment).absolutePath() + '/' + QFileInfo(alignment).completeBaseName() + ".csi",
        QFileInfo(alignment).absolutePath() + '/' + QFileInfo(alignment).completeBaseName() + ".crai"};
    for (const auto& candidate : candidates) if (QFileInfo::exists(candidate)) return candidate;
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

QStringList index_lines_for(const QPlainTextEdit* edit, const qsizetype count) {
    auto lines = edit->toPlainText().split('\n', Qt::KeepEmptyParts);
    for (auto& line : lines) line = line.trimmed();
    while (!lines.isEmpty() && lines.back().isEmpty() && lines.size() > count) lines.removeLast();
    lines.resize(count);
    return lines;
}

void append_resource_lines(QPlainTextEdit* edit, const QStringList& paths) {
    auto existing = edit->toPlainText();
    if (!existing.isEmpty() && !existing.endsWith('\n')) existing += '\n';
    edit->setPlainText(existing + paths.join('\n'));
    auto cursor = edit->textCursor();
    cursor.movePosition(QTextCursor::End);
    edit->setTextCursor(cursor);
}

QString resource_label(const QString& path) {
    const QUrl url(path);
    const auto name = path.startsWith("https://") ? QFileInfo(url.path()).fileName() : QFileInfo(path).fileName();
    return name.isEmpty() ? sanitized_resource_uri(path) : name;
}

}  // namespace

MainWindow::MainWindow() {
    setWindowTitle("BAM Seek — targeted BAM evidence");
    setAcceptDrops(true);
    resize(1240, 820);

    auto* root = new QWidget(this);
    auto* layout = new QVBoxLayout(root);
    auto* input_form = new QFormLayout();
    auto* bam_row = new QHBoxLayout();
    bam_path_ = new QPlainTextEdit(root);
    bam_path_->setPlaceholderText("One BAM/CRAM path or HTTPS URL per line");
    bam_path_->setMaximumHeight(70);
    bam_path_->setToolTip("Editable list: use one local path or HTTPS URL per line. Each BAM is evaluated separately.");
    auto* browse_bam = new QPushButton("Browse…", root);
    load_bams_button_ = new QPushButton("Load BAM(s)", root);
    load_bams_button_->setToolTip("Add the pending paths to the active BAM set used for VAF computation.");
    bam_row->addWidget(bam_path_);
    bam_row->addWidget(browse_bam);
    bam_row->addWidget(load_bams_button_);
    input_form->addRow("BAM / CRAM paths or URLs", bam_row);
    auto* index_row = new QHBoxLayout();
    index_path_ = new QPlainTextEdit(root);
    index_path_->setPlaceholderText("Optional: one index per BAM line; leave a line blank for automatic discovery");
    index_path_->setMaximumHeight(70);
    index_path_->setToolTip("Line N is the explicit index for BAM line N. A blank line uses automatic discovery.");
    auto* browse_index = new QPushButton("Browse…", root);
    index_row->addWidget(index_path_);
    index_row->addWidget(browse_index);
    input_form->addRow("Index", index_row);
    auto* loaded_row = new QHBoxLayout();
    loaded_bams_ = new QListWidget(root);
    loaded_bams_->setMaximumHeight(92);
    loaded_bams_->setAlternatingRowColors(true);
    loaded_bams_->setSelectionMode(QAbstractItemView::NoSelection);
    loaded_bams_->setToolTip("Only BAMs shown in this active list are included when evidence and VAF are computed.");
    clear_bams_button_ = new QPushButton("Clear all BAMs", root);
    clear_bams_button_->setEnabled(false);
    loaded_row->addWidget(loaded_bams_);
    loaded_row->addWidget(clear_bams_button_, 0, Qt::AlignTop);
    input_form->addRow("Actively loaded BAMs", loaded_row);
    auto* reference_row = new QHBoxLayout();
    reference_path_ = new QLineEdit(root);
    reference_path_->setPlaceholderText("Optional for BAM; required for CRAM (hg19 FASTA)");
    reference_path_->setClearButtonEnabled(true);
    reference_path_->setToolTip("Editable: type or paste the indexed reference FASTA path/URL. Browse is optional.");
    auto* browse_reference = new QPushButton("Browse…", root);
    reference_row->addWidget(reference_path_);
    reference_row->addWidget(browse_reference);
    input_form->addRow("Reference", reference_row);
    auto* clinical_mapping_row = new QHBoxLayout();
    clinical_mapping_path_ = new QLineEdit(root);
    clinical_mapping_path_->setPlaceholderText("Optional local TSV for resolving GENE c.… p.… to hg19 genomic alleles");
    clinical_mapping_path_->setClearButtonEnabled(true);
    clinical_mapping_path_->setToolTip("Editable: type or paste a local clinical-variant mapping TSV path. Browse is optional.");
    auto* browse_clinical_mapping = new QPushButton("Browse…", root);
    clinical_mapping_row->addWidget(clinical_mapping_path_);
    clinical_mapping_row->addWidget(browse_clinical_mapping);
    input_form->addRow("Clinical mapping TSV", clinical_mapping_row);
    layout->addLayout(input_form);

    auto* settings = new QHBoxLayout();
    molecule_mode_ = new QComboBox(root);
    molecule_mode_->addItem("Auto UMI; pair fallback", static_cast<int>(MoleculeMode::auto_detect));
    molecule_mode_->addItem("Read pairs / fragments", static_cast<int>(MoleculeMode::raw_reads));
    molecule_mode_->addItem("Selected UMI tag; pair fallback", static_cast<int>(MoleculeMode::selected_tag));
    molecule_tag_ = new QLineEdit(root);
    molecule_tag_->setPlaceholderText("e.g. MI");
    molecule_tag_->setEnabled(false);
    vaf_ = new QLineEdit("0.05", root);
    minimum_alt_reads_ = new QLineEdit("1", root);
    minimum_alt_molecules_ = new QLineEdit("1", root);
    minimum_mapq_ = new QLineEdit("20", root);
    minimum_baseq_ = new QLineEdit("20", root);
    vaf_->setValidator(new QDoubleValidator(0.0, 100.0, 6, vaf_));
    minimum_alt_reads_->setValidator(new QIntValidator(1, 1000000000, minimum_alt_reads_));
    minimum_alt_molecules_->setValidator(new QIntValidator(1, 1000000000, minimum_alt_molecules_));
    minimum_mapq_->setValidator(new QIntValidator(0, 255, minimum_mapq_));
    minimum_baseq_->setValidator(new QIntValidator(0, 255, minimum_baseq_));
    molecule_tag_->setMaxLength(2);
    molecule_tag_->setValidator(new QRegularExpressionValidator(QRegularExpression("[A-Za-z][A-Za-z0-9]"), molecule_tag_));
    settings->addWidget(new QLabel("Molecule grouping", root)); settings->addWidget(molecule_mode_);
    settings->addWidget(molecule_tag_);
    settings->addWidget(new QLabel("Min VAF (%)", root)); settings->addWidget(vaf_);
    settings->addWidget(new QLabel("Min alt reads", root)); settings->addWidget(minimum_alt_reads_);
    settings->addWidget(new QLabel("Min alt molecules", root)); settings->addWidget(minimum_alt_molecules_);
    settings->addWidget(new QLabel("Min mapQ", root)); settings->addWidget(minimum_mapq_);
    settings->addWidget(new QLabel("Min baseQ", root)); settings->addWidget(minimum_baseq_);
    layout->addLayout(settings);
    auto* flags = new QHBoxLayout();
    include_duplicates_ = new QCheckBox("Include duplicates", root);
    include_secondary_ = new QCheckBox("Include secondary alignments", root);
    include_supplementary_ = new QCheckBox("Include supplementary alignments", root);
    flags->addWidget(include_duplicates_); flags->addWidget(include_secondary_); flags->addWidget(include_supplementary_); flags->addStretch(1);
    layout->addLayout(flags);

    query_text_ = new QPlainTextEdit(root);
    query_text_->setPlaceholderText("One query per line:\nchr7:140453136 A>T\nchr1:100000-101000");
    query_text_->setPlainText("# One-based genomic or clinical notation\n# chr7:140453136 A>T\n# BRAF c.1799T>A p.V600E chr7:140453136 A>T");
    layout->addWidget(query_text_, 1);

    auto* buttons = new QHBoxLayout();
    run_button_ = new QPushButton("Query evidence", root);
    pileup_button_ = new QPushButton("View pileup", root);
    export_button_ = new QPushButton("Export audit JSON…", root);
    status_ = new QLabel("Enter an indexed BAM/CRAM, choose Load BAM(s), then enter queries.", root);
    buttons->addWidget(run_button_); buttons->addWidget(pileup_button_); buttons->addWidget(export_button_); buttons->addWidget(status_, 1);
    layout->addLayout(buttons);

    tabs_ = new QTabWidget(root);
    auto* evidence_tab = new QWidget(tabs_);
    auto* evidence_layout = new QVBoxLayout(evidence_tab);
    evidence_layout->setContentsMargins(0, 0, 0, 0);
    auto* splitter = new QSplitter(Qt::Vertical, evidence_tab);
    results_ = new QTableWidget(splitter);
    results_->setColumnCount(19);
    results_->setHorizontalHeaderLabels({"BAM / CRAM", "Query", "Status", "REF reads", "ALT reads", "REF+ALT reads", "Read VAF",
        "REF molecules", "ALT molecules", "REF+ALT molecules", "Molecule VAF", "OTHER/N reads", "Ambiguous molecules",
        "ALT Fwd", "ALT Rev", "Strand P", "Molecule grouping", "Missing tags", "Evidence summary"});
    results_->setSelectionBehavior(QAbstractItemView::SelectRows);
    results_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    results_->horizontalHeader()->setStretchLastSection(true);
    read_details_ = new QPlainTextEdit(splitter);
    read_details_->setReadOnly(true);
    read_details_->setPlaceholderText("Select a variant result to inspect supporting reads.");
    splitter->addWidget(results_);
    splitter->addWidget(read_details_);
    splitter->setSizes({330, 220});
    evidence_layout->addWidget(splitter);
    tabs_->addTab(evidence_tab, "Evidence");
    auto* pileup_tab = new QWidget(tabs_);
    auto* pileup_layout = new QVBoxLayout(pileup_tab);
    auto* pileup_controls = new QHBoxLayout();
    group_pairs_ = new QCheckBox("Group and link read pairs", pileup_tab);
    group_pairs_->setChecked(true);
    pileup_controls->addWidget(group_pairs_);
    pileup_controls->addWidget(new QLabel("Blue: forward   Red: reverse   Yellow: mismatch   Green: indel   Gray: low base quality", pileup_tab));
    pileup_controls->addStretch(1);
    pileup_layout->addLayout(pileup_controls);
    pileup_summary_ = new QLabel(pileup_tab);
    pileup_summary_->setWordWrap(true);
    pileup_summary_->setTextInteractionFlags(Qt::TextSelectableByMouse);
    pileup_summary_->setText("Select a targeted variant result and choose View pileup.");
    pileup_layout->addWidget(pileup_summary_);
    pileup_scroll_ = new QScrollArea(pileup_tab);
    pileup_scroll_->setAlignment(Qt::AlignLeft | Qt::AlignTop);
    pileup_view_ = new PileupView(pileup_scroll_);
    pileup_scroll_->setWidget(pileup_view_);
    pileup_scroll_->setWidgetResizable(false);
    pileup_layout->addWidget(pileup_scroll_);
    tabs_->addTab(pileup_tab, "Pileup");
    auto* broadcast_tab = new QWidget(tabs_);
    auto* broadcast_layout = new QVBoxLayout(broadcast_tab);
    auto* receiver_controls = new QHBoxLayout();
    receiver_enabled_ = new QCheckBox("Listen for IGV commands", broadcast_tab);
    receiver_enabled_->setChecked(true);
    receiver_port_ = new QSpinBox(broadcast_tab);
    receiver_port_->setRange(1024, 65535);
    receiver_port_->setValue(60151);
    receiver_port_->setToolTip("IGV's default command port is 60151. Only one application can listen on a port at a time.");
    receiver_status_ = new QLabel(broadcast_tab);
    receiver_controls->addWidget(receiver_enabled_);
    receiver_controls->addWidget(new QLabel("Local port", broadcast_tab));
    receiver_controls->addWidget(receiver_port_);
    receiver_controls->addWidget(receiver_status_, 1);
    broadcast_layout->addLayout(receiver_controls);
    broadcast_text_ = new QPlainTextEdit(broadcast_tab);
    broadcast_text_->setPlaceholderText("Received IGV /load and /goto information will be appended here. This box remains editable; received commands never load files or change the current query.");
    broadcast_text_->setToolTip("Passive inbox. Editing this text does not trigger any BAM Seek action.");
    broadcast_layout->addWidget(broadcast_text_);
    tabs_->addTab(broadcast_tab, "Broadcast inbox");
    layout->addWidget(tabs_, 2);
    setCentralWidget(root);

    connect(browse_bam, &QPushButton::clicked, this, [this] { choose_bam(); });
    connect(browse_index, &QPushButton::clicked, this, [this] { choose_index(); });
    connect(load_bams_button_, &QPushButton::clicked, this, [this] { load_bams(); });
    connect(clear_bams_button_, &QPushButton::clicked, this, [this] { clear_loaded_bams(); });
    connect(browse_reference, &QPushButton::clicked, this, [this] { choose_reference(); });
    connect(browse_clinical_mapping, &QPushButton::clicked, this, [this] { choose_clinical_mapping(); });
    connect(run_button_, &QPushButton::clicked, this, [this] { run_queries(); });
    connect(pileup_button_, &QPushButton::clicked, this, [this] { show_pileup(); });
    connect(export_button_, &QPushButton::clicked, this, [this] { export_audit(); });
    connect(molecule_mode_, &QComboBox::currentIndexChanged, this, [this] {
        molecule_tag_->setEnabled(molecule_mode_->currentData().toInt() == static_cast<int>(MoleculeMode::selected_tag));
    });
    connect(results_, &QTableWidget::cellClicked, this, [this](const int row, const int column) { show_read_details(row, column); });
    connect(&watcher_, &QFutureWatcher<MultiBamBatch>::finished, this, [this] { show_results(); });
    connect(&pileup_watcher_, &QFutureWatcher<PileupLoad>::finished, this, [this] { pileup_loaded(); });
    connect(&audit_watcher_, &QFutureWatcher<AuditSave>::finished, this, [this] { audit_saved(); });
    connect(group_pairs_, &QCheckBox::toggled, this, [this](const bool enabled) { pileup_view_->set_group_pairs(enabled); });
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
            receiver_status_->setText(listening ? QString("Listening passively on 127.0.0.1:%1").arg(port)
                                                : (error.isEmpty() ? "Receiver stopped" : "Receiver unavailable: " + error));
        });
    set_receiver_enabled(true);
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
    const auto urls = event->mimeData()->urls();
    QStringList files;
    for (const auto& url : urls) if (url.isLocalFile()) files.append(url.toLocalFile());
    if (files.isEmpty()) return;
    append_resource_lines(bam_path_, files);
    status_->setText(QString("%1 pending BAM/CRAM path(s) added. Choose Load BAM(s) to activate them.").arg(files.size()));
    event->acceptProposedAction();
}

void MainWindow::choose_bam() {
    const auto files = QFileDialog::getOpenFileNames(this, "Choose indexed BAM or CRAM files", {}, "Alignment files (*.bam *.cram *.sam);;All files (*)");
    if (!files.isEmpty()) append_resource_lines(bam_path_, files);
}

void MainWindow::choose_index() {
    const auto files = QFileDialog::getOpenFileNames(this, "Choose alignment indexes in BAM order", {}, "Alignment indexes (*.bai *.csi *.crai);;All files (*)");
    if (!files.isEmpty()) append_resource_lines(index_path_, files);
}

void MainWindow::load_bams() {
    if (watcher_.isRunning() || pileup_watcher_.isRunning() || audit_watcher_.isRunning()) return;
    const auto pending_bams = nonempty_lines(bam_path_);
    const auto pending_indexes = index_lines_for(index_path_, pending_bams.size());
    if (pending_bams.isEmpty()) {
        QMessageBox::information(this, "No pending BAM", "Type, paste, browse, or drop a BAM/CRAM path, then choose Load BAM(s).");
        return;
    }

    QStringList issues;
    int added = 0;
    int updated = 0;
    for (qsizetype pending_index = 0; pending_index < pending_bams.size(); ++pending_index) {
        const auto bam = pending_bams[pending_index];
        auto index = pending_indexes[pending_index];
        const auto remote_bam = bam.startsWith("https://");
        if (bam.startsWith("http://") || (!remote_bam && bam.contains("://"))) {
            issues.append(resource_label(bam) + ": only local paths and HTTPS URLs are accepted.");
            continue;
        }
        if (!remote_bam && !QFileInfo::exists(bam)) {
            issues.append(bam + ": alignment file not found.");
            continue;
        }
        if (!index.isEmpty()) {
            const auto remote_index = index.startsWith("https://");
            if (index.startsWith("http://") || (!remote_index && index.contains("://"))) {
                issues.append(resource_label(bam) + ": index must be a local path or HTTPS URL.");
                continue;
            }
            if (!remote_index && !QFileInfo::exists(index)) {
                issues.append(resource_label(bam) + ": explicit index file not found.");
                continue;
            }
        } else if (!remote_bam) {
            index = index_path_for(bam);
            if (index.isEmpty()) {
                issues.append(resource_label(bam) + ": no accompanying .bai, .csi, or .crai index was found.");
                continue;
            }
        }

        const auto existing = loaded_bam_paths_.indexOf(bam);
        if (existing >= 0) {
            if (!index.isEmpty() || loaded_index_paths_[existing].isEmpty()) loaded_index_paths_[existing] = index;
            ++updated;
        } else {
            loaded_bam_paths_.append(bam);
            loaded_index_paths_.append(index);
            ++added;
        }
    }

    refresh_loaded_bams();
    if (added > 0 || updated > 0) {
        bam_path_->clear();
        index_path_->clear();
    }
    status_->setText(QString("%1 active BAM(s); %2 added, %3 updated. VAF queries use this active set.")
        .arg(loaded_bam_paths_.size()).arg(added).arg(updated));
    if (!issues.isEmpty()) QMessageBox::warning(this, "Some BAMs were not loaded", issues.join('\n'));
}

void MainWindow::clear_loaded_bams() {
    if (watcher_.isRunning() || pileup_watcher_.isRunning() || audit_watcher_.isRunning()) return;
    loaded_bam_paths_.clear();
    loaded_index_paths_.clear();
    configured_bam_paths_.clear();
    configured_index_paths_.clear();
    result_bam_paths_.clear();
    result_index_paths_.clear();
    last_batch_ = {};
    results_->setRowCount(0);
    read_details_->clear();
    pileup_view_->set_data({});
    refresh_loaded_bams();
    status_->setText("All BAMs cleared. Load at least one BAM before computing VAF.");
}

void MainWindow::refresh_loaded_bams() {
    loaded_bams_->clear();
    for (qsizetype index = 0; index < loaded_bam_paths_.size(); ++index) {
        const auto& bam = loaded_bam_paths_[index];
        const auto& bam_index = loaded_index_paths_[index];
        auto* item = new QListWidgetItem(resource_label(bam)
            + (bam_index.isEmpty() ? "  —  index: automatic" : "  —  index: " + resource_label(bam_index)), loaded_bams_);
        item->setToolTip("BAM/CRAM: " + sanitized_resource_uri(bam) + "\nIndex: "
            + (bam_index.isEmpty() ? "automatic" : sanitized_resource_uri(bam_index)));
    }
    clear_bams_button_->setEnabled(!loaded_bam_paths_.isEmpty() && !watcher_.isRunning());
}

void MainWindow::choose_reference() {
    const auto file = QFileDialog::getOpenFileName(this, "Choose hg19 reference FASTA", {}, "FASTA files (*.fa *.fasta *.fna);;All files (*)");
    if (!file.isEmpty()) reference_path_->setText(file);
}

void MainWindow::choose_clinical_mapping() {
    const auto file = QFileDialog::getOpenFileName(this, "Choose local clinical mapping TSV", {}, "TSV files (*.tsv *.txt);;All files (*)");
    if (!file.isEmpty()) clinical_mapping_path_->setText(file);
}

FilterSettings MainWindow::filters() const {
    FilterSettings values;
    values.molecule_mode = static_cast<MoleculeMode>(molecule_mode_->currentData().toInt());
    values.molecule_tag = molecule_tag_->text().trimmed().toStdString();
    values.minimum_variant_allele_fraction = std::max(0.0, vaf_->text().toDouble() / 100.0);
    values.minimum_alternate_reads = std::max(1, minimum_alt_reads_->text().toInt());
    values.minimum_alternate_molecules = std::max(1, minimum_alt_molecules_->text().toInt());
    values.minimum_mapping_quality = std::max(0, minimum_mapq_->text().toInt());
    values.minimum_base_quality = std::max(0, minimum_baseq_->text().toInt());
    values.include_duplicates = include_duplicates_->isChecked();
    values.include_secondary = include_secondary_->isChecked();
    values.include_supplementary = include_supplementary_->isChecked();
    return values;
}

void MainWindow::run_queries() {
    const auto bams = loaded_bam_paths_;
    const auto indexes = loaded_index_paths_;
    if (bams.isEmpty()) {
        QMessageBox::warning(this, "No active alignment file", "Type or browse to a BAM/CRAM, choose Load BAM(s), and then run the query.");
        return;
    }
    const auto has_insecure_resource = [](const QStringList& resources) {
        return std::any_of(resources.cbegin(), resources.cend(), [](const QString& value) { return value.startsWith("http://"); });
    };
    if (has_insecure_resource(bams) || has_insecure_resource(indexes) || reference_path_->text().startsWith("http://")) {
        QMessageBox::warning(this, "Secure remote access", "BAM Seek accepts local paths or HTTPS resources only. HTTP URLs are not permitted.");
        return;
    }
    if (!vaf_->hasAcceptableInput() || !minimum_alt_reads_->hasAcceptableInput() || !minimum_alt_molecules_->hasAcceptableInput()
        || !minimum_mapq_->hasAcceptableInput() || !minimum_baseq_->hasAcceptableInput()
        || (molecule_mode_->currentData().toInt() == static_cast<int>(MoleculeMode::selected_tag) && !molecule_tag_->hasAcceptableInput())) {
        QMessageBox::warning(this, "Invalid filters", "Correct the highlighted numeric filters and enter a valid two-character BAM tag.");
        return;
    }
    const auto loaded_mappings = load_clinical_mappings(clinical_mapping_path_->text().trimmed().toStdString());
    auto parsed = parse_queries(query_text_->toPlainText().toStdString(), loaded_mappings.mappings);
    parsed.errors.insert(parsed.errors.begin(), loaded_mappings.errors.begin(), loaded_mappings.errors.end());
    if (parsed.queries.empty()) {
        QMessageBox::warning(this, "No valid queries", QString::fromStdString(parsed.errors.empty() ? "Enter at least one query." : parsed.errors.front()));
        return;
    }
    last_filters_ = filters();
    configured_bam_paths_ = bams;
    configured_index_paths_ = indexes;
    last_reference_path_ = reference_path_->text().trimmed();
    last_clinical_mapping_path_ = clinical_mapping_path_->text().trimmed();
    last_query_text_ = query_text_->toPlainText();
    const auto reference = last_reference_path_.toStdString();
    run_button_->setEnabled(false);
    load_bams_button_->setEnabled(false);
    clear_bams_button_->setEnabled(false);
    status_->setText(QString("Querying evidence in %1 alignment(s)…").arg(bams.size()));
    watcher_.setFuture(QtConcurrent::run([queries = parsed.queries, errors = parsed.errors, filter_values = last_filters_, bams, indexes, reference] {
        MultiBamBatch combined;
        combined.batch.errors = errors;
        for (qsizetype source_index = 0; source_index < bams.size(); ++source_index) {
            const auto bam_qt = bams[source_index];
            const auto index_qt = indexes[source_index];
            const auto bam = bam_qt.toStdString();
            const auto index = index_qt.toStdString();
            try {
                igv::Resource resource{.uri = bam};
                if (!index.empty()) resource.index_uri = index;
                if (!reference.empty()) resource.reference_uri = reference;
                EvidenceEngine engine(std::move(resource));
                auto evaluated = engine.evaluate(queries, filter_values);
                for (auto& result : evaluated.results) {
                    combined.batch.results.push_back(std::move(result));
                    combined.result_bams.append(bam_qt);
                    combined.result_indexes.append(index_qt);
                }
                for (auto& error : evaluated.errors) {
                    combined.batch.errors.push_back("[" + resource_label(bam_qt).toStdString() + "] " + error);
                }
            } catch (const std::exception& error) {
                combined.batch.errors.push_back("[" + resource_label(bam_qt).toStdString() + "] "
                    + sanitized_error(error.what(), {bam, index, reference}));
            }
        }
        return combined;
    }));
}

void MainWindow::show_results() {
    const auto combined = watcher_.result();
    last_batch_ = combined.batch;
    result_bam_paths_ = combined.result_bams;
    result_index_paths_ = combined.result_indexes;
    results_->setRowCount(0);
    for (int index = 0; index < static_cast<int>(last_batch_.results.size()); ++index) {
        results_->insertRow(index);
        const auto& result = last_batch_.results[static_cast<std::size_t>(index)];
        if (const auto* evidence = std::get_if<VariantEvidence>(&result)) {
            const auto& count = evidence->counts;
            const auto strand_p = count.strand_bias_p_value();
            const QList<QString> values{resource_label(result_bam_paths_[index]), display_query(evidence->query), evidence->passes_thresholds ? "PRESENT" : "not detected",
                QString::number(count.reference_reads), QString::number(count.alternate_reads), QString::number(count.informative_read_depth()),
                QString::number(count.allele_fraction() * 100.0, 'f', 4) + '%',
                evidence->molecule_counts_available ? QString::number(count.reference_molecules) : "N/A",
                evidence->molecule_counts_available ? QString::number(count.alternate_molecules) : "N/A",
                evidence->molecule_counts_available ? QString::number(count.molecule_depth()) : "N/A",
                evidence->molecule_counts_available ? QString::number(count.molecule_allele_fraction() * 100.0, 'f', 4) + '%' : "N/A",
                QString::number(count.other_reads), evidence->molecule_counts_available ? QString::number(count.other_molecules) : "N/A",
                QString::number(count.alternate_forward_reads), QString::number(count.alternate_reverse_reads),
                strand_p ? QString::number(*strand_p, 'g', 3) : "N/A",
                QString::fromStdString(evidence->molecule_counts_available ? evidence->molecule_tag_used : "Unavailable"),
                QString::number(evidence->reads_missing_molecule_tag), evidence_summary(*evidence)};
            for (int column = 0; column < values.size(); ++column) {
                auto* item = new QTableWidgetItem(values[column]);
                if (column == 0) item->setToolTip(sanitized_resource_uri(result_bam_paths_[index]));
                results_->setItem(index, column, item);
            }
        } else {
            const auto& region = std::get<RegionEvidence>(result);
            QList<QString> values(19);
            values[0] = resource_label(result_bam_paths_[index]);
            values[1] = QString::fromStdString(region.query.source_text);
            values[2] = region.candidates.empty() ? "NO CANDIDATES" : "CANDIDATES";
            values[18] = QString::fromStdString(region.note);
            for (int column = 0; column < values.size(); ++column) {
                auto* item = new QTableWidgetItem(values[column]);
                if (column == 0) item->setToolTip(sanitized_resource_uri(result_bam_paths_[index]));
                results_->setItem(index, column, item);
            }
        }
    }
    run_button_->setEnabled(true);
    load_bams_button_->setEnabled(true);
    clear_bams_button_->setEnabled(!loaded_bam_paths_.isEmpty());
    status_->setText(QString("%1 result(s) from %2 active BAM(s), %3 issue(s). Select a row for read evidence.")
        .arg(last_batch_.results.size()).arg(loaded_bam_paths_.size()).arg(last_batch_.errors.size()));
    if (!last_batch_.errors.empty()) read_details_->setPlainText("Issues:\n" + QString::fromStdString([&] {
        std::string combined;
        for (const auto& error : last_batch_.errors) combined += error + '\n';
        return combined;
    }()));
}

void MainWindow::show_read_details(const int row, const int) {
    if (row < 0 || static_cast<std::size_t>(row) >= last_batch_.results.size()) return;
    const auto* evidence = std::get_if<VariantEvidence>(&last_batch_.results[static_cast<std::size_t>(row)]);
    if (evidence == nullptr) {
        const auto& region = std::get<RegionEvidence>(last_batch_.results[static_cast<std::size_t>(row)]);
        QString text = "Alignment: " + sanitized_resource_uri(result_bam_paths_[row]) + "\n" + QString::fromStdString(region.note) + "\n\n";
        for (const auto& candidate : region.candidates) {
            text += evidence_summary(candidate) + '\n';
        }
        read_details_->setPlainText(text);
        return;
    }
    QString text = "Alignment: " + sanitized_resource_uri(result_bam_paths_[row]) + "\n" + evidence_summary(*evidence)
        + "\nMolecule grouping: " + QString::fromStdString(evidence->molecule_counts_available ? evidence->molecule_tag_used : "unavailable") + "\n\n";
    for (const auto& read : evidence->reads) {
        text += QString::fromStdString(read.read_name) + (read.reverse_strand ? "  reverse  " : "  forward  ")
            + QString::fromStdString(read.summary);
        if (!read.molecule_id.empty()) text += "  molecule=" + QString::fromStdString(read.molecule_id);
        text += '\n';
    }
    read_details_->setPlainText(text);
}

void MainWindow::show_pileup() {
    if (pileup_watcher_.isRunning()) return;
    const int row = results_->currentRow();
    if (row < 0 || static_cast<std::size_t>(row) >= last_batch_.results.size()) {
        QMessageBox::information(this, "Select a variant", "Select a targeted variant result, then choose View pileup.");
        return;
    }
    const auto* evidence = std::get_if<VariantEvidence>(&last_batch_.results[static_cast<std::size_t>(row)]);
    if (evidence == nullptr) {
        QMessageBox::information(this, "Select a targeted variant", "Pileup is available for targeted variant rows. Region candidates can be pasted into the query box.");
        return;
    }
    const auto bam = result_bam_paths_[row].toStdString();
    const auto index = result_index_paths_[row].toStdString();
    const auto reference = last_reference_path_.toStdString();
    const auto query = evidence->query;
    const auto filter_values = last_filters_;
    const auto summary = evidence_summary(*evidence);
    pileup_button_->setEnabled(false);
    load_bams_button_->setEnabled(false);
    clear_bams_button_->setEnabled(false);
    status_->setText("Loading local alignment pileup…");
    pileup_watcher_.setFuture(QtConcurrent::run([bam, index, reference, query, filter_values, summary] {
        PileupLoad loaded;
        loaded.summary = summary;
        try {
            igv::Resource resource{.uri = bam};
            if (!index.empty()) resource.index_uri = index;
            if (!reference.empty()) resource.reference_uri = reference;
            EvidenceEngine engine(std::move(resource));
            loaded.data = engine.pileup(query, filter_values, 40);
        } catch (const std::exception& error) {
            loaded.error = sanitized_error(error.what(), {bam, index, reference});
        }
        return loaded;
    }));
}

void MainWindow::pileup_loaded() {
    const auto loaded = pileup_watcher_.result();
    pileup_button_->setEnabled(true);
    load_bams_button_->setEnabled(true);
    clear_bams_button_->setEnabled(!loaded_bam_paths_.isEmpty());
    if (!loaded.error.empty()) {
        QMessageBox::warning(this, "Pileup unavailable", QString::fromStdString(loaded.error));
        status_->setText("Could not load pileup.");
        pileup_summary_->setText("Pileup unavailable: " + QString::fromStdString(loaded.error));
        return;
    }
    pileup_view_->set_data(loaded.data);
    pileup_summary_->setText(loaded.summary);
    tabs_->setCurrentIndex(1);
    pileup_scroll_->ensureVisible(pileup_view_->variant_x(), 0, 120, 20);
    status_->setText(QString("Rendered %1 of %2 filtered alignments locally%3.")
        .arg(loaded.data.alignments.size()).arg(loaded.data.total_alignments).arg(loaded.data.truncated ? " (display limit reached)" : ""));
}

void MainWindow::export_audit() {
    if (audit_watcher_.isRunning()) return;
    if (last_batch_.results.empty() && last_batch_.errors.empty()) {
        QMessageBox::information(this, "Nothing to export", "Run an evidence query before exporting an audit record.");
        return;
    }
    const auto path = QFileDialog::getSaveFileName(this, "Export audit JSON", "bam-seek-audit.json", "JSON files (*.json)");
    if (path.isEmpty()) return;
    export_button_->setEnabled(false);
    load_bams_button_->setEnabled(false);
    clear_bams_button_->setEnabled(false);
    status_->setText("Creating audit fingerprints in the background…");
    QJsonObject root{{"application", "BAM Seek"}, {"version", QStringLiteral(BAM_SEEK_VERSION)}, {"generated_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs)},
        {"query_input", last_query_text_},
        {"molecule_mode", mode_name(last_filters_.molecule_mode)},
        {"molecule_tag", QString::fromStdString(last_filters_.molecule_tag)},
        {"minimum_vaf", last_filters_.minimum_variant_allele_fraction}, {"minimum_alt_reads", last_filters_.minimum_alternate_reads},
        {"minimum_alt_molecules", last_filters_.minimum_alternate_molecules}, {"minimum_mapq", last_filters_.minimum_mapping_quality}, {"minimum_baseq", last_filters_.minimum_base_quality},
        {"include_duplicates", last_filters_.include_duplicates}, {"include_secondary", last_filters_.include_secondary}, {"include_supplementary", last_filters_.include_supplementary}};
    QJsonArray result_array;
    const auto variant_json = [](const VariantEvidence& evidence) {
        QJsonObject row{{"query", display_query(evidence.query)}, {"present", evidence.passes_thresholds},
            {"depth", evidence.counts.depth()}, {"total_called_reads", evidence.counts.depth()},
            {"informative_read_depth", evidence.counts.informative_read_depth()},
            {"gene", QString::fromStdString(evidence.query.gene)}, {"transcript", QString::fromStdString(evidence.query.transcript)},
            {"coding_change", QString::fromStdString(evidence.query.coding_change)}, {"protein_change", QString::fromStdString(evidence.query.protein_change)},
            {"contig", QString::fromStdString(evidence.query.contig)}, {"position_one_based", evidence.query.position + 1},
            {"reference_allele", QString::fromStdString(evidence.query.reference)}, {"alternate_allele", QString::fromStdString(evidence.query.alternate)},
            {"ref_reads", evidence.counts.reference_reads}, {"other_reads", evidence.counts.other_reads},
            {"ref_forward_reads", evidence.counts.reference_forward_reads}, {"ref_reverse_reads", evidence.counts.reference_reverse_reads},
            {"alt_reads", evidence.counts.alternate_reads}, {"alt_forward_reads", evidence.counts.alternate_forward_reads},
            {"alt_reverse_reads", evidence.counts.alternate_reverse_reads}, {"read_vaf", evidence.counts.allele_fraction()},
            {"vaf", evidence.counts.allele_fraction()},
            {"molecule_counts_available", evidence.molecule_counts_available}, {"reads_missing_molecule_tag", evidence.reads_missing_molecule_tag},
            {"molecule_tag", QString::fromStdString(evidence.molecule_tag_used)},
            {"molecule_grouping", QString::fromStdString(evidence.molecule_tag_used)}, {"evidence_summary", evidence_summary(evidence)}};
        row["alt_molecules"] = evidence.molecule_counts_available ? QJsonValue(evidence.counts.alternate_molecules) : QJsonValue::Null;
        row["ref_molecules"] = evidence.molecule_counts_available ? QJsonValue(evidence.counts.reference_molecules) : QJsonValue::Null;
        row["other_molecules"] = evidence.molecule_counts_available ? QJsonValue(evidence.counts.other_molecules) : QJsonValue::Null;
        row["informative_molecule_depth"] = evidence.molecule_counts_available ? QJsonValue(evidence.counts.molecule_depth()) : QJsonValue::Null;
        row["molecule_vaf"] = evidence.molecule_counts_available ? QJsonValue(evidence.counts.molecule_allele_fraction()) : QJsonValue::Null;
        const auto strand_p = evidence.counts.strand_bias_p_value();
        row["strand_bias_fisher_p"] = strand_p ? QJsonValue(*strand_p) : QJsonValue::Null;
        QJsonArray reads;
        for (const auto& read : evidence.reads) reads.append(QJsonObject{{"name", QString::fromStdString(read.read_name)}, {"allele", QString::fromStdString(allele_name(read.allele))}, {"reverse", read.reverse_strand}, {"mapq", read.mapping_quality}, {"minimum_baseq", read.minimum_base_quality}, {"molecule", QString::fromStdString(read.molecule_id)}, {"summary", QString::fromStdString(read.summary)}});
        row["reads"] = reads;
        return row;
    };
    for (std::size_t result_index = 0; result_index < last_batch_.results.size(); ++result_index) {
        const auto& item = last_batch_.results[result_index];
        const auto source = sanitized_resource_uri(result_bam_paths_[static_cast<qsizetype>(result_index)]);
        if (const auto* evidence = std::get_if<VariantEvidence>(&item)) {
            auto row = variant_json(*evidence);
            row["alignment"] = source;
            result_array.append(row);
        } else {
            const auto& region = std::get<RegionEvidence>(item);
            QJsonObject row{{"alignment", source}, {"region", QString::fromStdString(region.query.source_text)}, {"note", QString::fromStdString(region.note)}};
            QJsonArray candidates;
            for (const auto& candidate : region.candidates) candidates.append(variant_json(candidate));
            row["candidates"] = candidates;
            result_array.append(row);
        }
    }
    root["results"] = result_array;
    QJsonArray errors;
    for (const auto& error : last_batch_.errors) errors.append(QString::fromStdString(error));
    root["errors"] = errors;
    const auto alignments = configured_bam_paths_;
    auto indexes = configured_index_paths_;
    for (qsizetype source_index = 0; source_index < alignments.size(); ++source_index) {
        if (indexes[source_index].isEmpty() && !alignments[source_index].startsWith("https://")) {
            indexes[source_index] = index_path_for(alignments[source_index]);
        }
    }
    const auto reference = last_reference_path_;
    const auto clinical_mapping = last_clinical_mapping_path_;
    audit_watcher_.setFuture(QtConcurrent::run([root = std::move(root), alignments, indexes, reference, clinical_mapping, path]() mutable {
        AuditSave result{path, {}};
        QJsonArray alignment_array;
        for (qsizetype source_index = 0; source_index < alignments.size(); ++source_index) {
            QJsonObject resource{{"alignment", file_identity(alignments[source_index])}};
            resource["index"] = indexes[source_index].isEmpty()
                ? QJsonObject{{"auto_detected", true}}
                : file_identity(indexes[source_index]);
            alignment_array.append(resource);
        }
        root["alignments"] = alignment_array;
        root["reference"] = file_identity(reference);
        root["clinical_mapping"] = file_identity(clinical_mapping);
        QSaveFile output(path);
        if (!output.open(QIODevice::WriteOnly)) {
            result.error = "Could not open the selected audit file for writing.";
            return result;
        }
        const auto document = QJsonDocument(root).toJson(QJsonDocument::Indented);
        if (output.write(document) != document.size() || !output.commit()) {
            result.error = "Could not atomically save the complete audit file.";
        }
        return result;
    }));
}

void MainWindow::audit_saved() {
    const auto result = audit_watcher_.result();
    export_button_->setEnabled(true);
    load_bams_button_->setEnabled(true);
    clear_bams_button_->setEnabled(!loaded_bam_paths_.isEmpty());
    if (!result.error.empty()) {
        QMessageBox::critical(this, "Export failed", QString::fromStdString(result.error));
        status_->setText("Audit export failed.");
        return;
    }
    status_->setText("Audit export saved: " + result.path);
}

}  // namespace bamseek
