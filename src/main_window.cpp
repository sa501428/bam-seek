#include <bamseek/main_window.hpp>

#include <bamseek/igv_command_receiver.hpp>
#include <bamseek/query.hpp>

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QCryptographicHash>
#include <QDateTime>
#include <QDragEnterEvent>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFuture>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QIntValidator>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QMimeData>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSaveFile>
#include <QSpinBox>
#include <QSplitter>
#include <QTabWidget>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QUrl>
#include <QVBoxLayout>
#include <QWidget>
#include <QtConcurrentRun>

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
    const auto query = display_query(evidence.query);
    const auto molecule_vaf = evidence.molecule_counts_available
        ? percent(counts.molecule_allele_fraction())
        : QString("unavailable");
    return QString("In %1, %2 was supported by %3 of %4 informative reads (VAF %5; %6 REF and %7 OTHER/N reads). "
                   "After collapsing reads with the same read name into paired fragments, %8 of %9 informative molecules supported the variant "
                   "(molecule VAF %10; %11 ambiguous fragments). ALT support included %12 forward and %13 reverse reads.")
        .arg(bam_label, query)
        .arg(counts.alternate_reads)
        .arg(counts.informative_read_depth())
        .arg(percent(counts.allele_fraction()))
        .arg(counts.reference_reads)
        .arg(counts.other_reads)
        .arg(evidence.molecule_counts_available ? QString::number(counts.alternate_molecules) : "N/A")
        .arg(evidence.molecule_counts_available ? QString::number(counts.molecule_depth()) : "N/A")
        .arg(molecule_vaf)
        .arg(evidence.molecule_counts_available ? QString::number(counts.other_molecules) : "N/A")
        .arg(counts.alternate_forward_reads)
        .arg(counts.alternate_reverse_reads);
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

}  // namespace

MainWindow::MainWindow() {
    setWindowTitle("BAM Seek — variant allele frequency");
    setAcceptDrops(true);
    resize(1220, 820);

    auto* root = new QWidget(this);
    auto* layout = new QVBoxLayout(root);

    auto* bam_group = new QGroupBox("1. Add BAMs", root);
    auto* bam_group_layout = new QVBoxLayout(bam_group);
    auto* bam_row = new QHBoxLayout();
    bam_path_ = new QPlainTextEdit(bam_group);
    bam_path_->setPlaceholderText("One local path or HTTPS BAM URL per line, or drag local BAMs here");
    bam_path_->setMaximumHeight(62);
    auto* browse_bam = new QPushButton("Choose BAMs…", bam_group);
    load_bams_button_ = new QPushButton("Add to analysis", bam_group);
    bam_row->addWidget(bam_path_, 1);
    bam_row->addWidget(browse_bam);
    bam_row->addWidget(load_bams_button_);
    bam_group_layout->addLayout(bam_row);
    auto* loaded_row = new QHBoxLayout();
    loaded_bams_ = new QListWidget(bam_group);
    loaded_bams_->setMaximumHeight(82);
    loaded_bams_->setAlternatingRowColors(true);
    loaded_bams_->setSelectionMode(QAbstractItemView::NoSelection);
    clear_bams_button_ = new QPushButton("Clear BAMs", bam_group);
    clear_bams_button_->setEnabled(false);
    loaded_row->addWidget(loaded_bams_, 1);
    loaded_row->addWidget(clear_bams_button_, 0, Qt::AlignTop);
    bam_group_layout->addLayout(loaded_row);
    layout->addWidget(bam_group);

    auto* query_group = new QGroupBox("2. Enter variants", root);
    auto* query_layout = new QVBoxLayout(query_group);
    query_text_ = new QPlainTextEdit(query_group);
    query_text_->setPlaceholderText("One variant per line");
    query_text_->setPlainText(
        "# Accepted examples (one-based coordinates)\n"
        "# chr7:140453136 A>T\n"
        "# chr7:g.140453136A>T\n"
        "# chr7 140453136 A T\n"
        "# chr7 140453136 . A T\n"
        "# BRAF c.1799T>A p.V600E chr7:140453136 A>T");
    query_text_->setMaximumHeight(150);
    query_layout->addWidget(query_text_);
    auto* filter_row = new QHBoxLayout();
    minimum_mapq_ = new QLineEdit("20", query_group);
    minimum_baseq_ = new QLineEdit("20", query_group);
    minimum_mapq_->setMaximumWidth(60);
    minimum_baseq_->setMaximumWidth(60);
    minimum_mapq_->setValidator(new QIntValidator(0, 255, minimum_mapq_));
    minimum_baseq_->setValidator(new QIntValidator(0, 255, minimum_baseq_));
    filter_row->addWidget(new QLabel("Minimum mapQ", query_group));
    filter_row->addWidget(minimum_mapq_);
    filter_row->addSpacing(16);
    filter_row->addWidget(new QLabel("Minimum baseQ", query_group));
    filter_row->addWidget(minimum_baseq_);
    filter_row->addStretch(1);
    filter_row->addWidget(new QLabel("Molecules are paired fragments grouped by read name (no UMI).", query_group));
    query_layout->addLayout(filter_row);
    layout->addWidget(query_group);

    auto* action_row = new QHBoxLayout();
    run_button_ = new QPushButton("Calculate VAFs", root);
    export_button_ = new QPushButton("Export audit JSON…", root);
    status_ = new QLabel("Add one or more BAMs, enter variants, then calculate VAFs.", root);
    action_row->addWidget(run_button_);
    action_row->addWidget(export_button_);
    action_row->addWidget(status_, 1);
    layout->addLayout(action_row);

    tabs_ = new QTabWidget(root);
    auto* results_tab = new QWidget(tabs_);
    auto* results_layout = new QVBoxLayout(results_tab);
    results_layout->setContentsMargins(0, 0, 0, 0);
    auto* splitter = new QSplitter(Qt::Vertical, results_tab);
    results_ = new QTableWidget(splitter);
    results_->setColumnCount(12);
    results_->setHorizontalHeaderLabels({"BAM", "Variant", "Observed", "ALT reads", "Read depth", "VAF",
        "ALT molecules", "Molecule depth", "Molecule VAF", "OTHER/N reads", "Ambiguous molecules", "Summary"});
    results_->setSelectionBehavior(QAbstractItemView::SelectRows);
    results_->setSelectionMode(QAbstractItemView::SingleSelection);
    results_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    results_->horizontalHeader()->setStretchLastSection(true);
    results_->verticalHeader()->setVisible(false);

    auto* details = new QWidget(splitter);
    auto* details_layout = new QVBoxLayout(details);
    details_layout->setContentsMargins(0, 0, 0, 0);
    auto* summary_row = new QHBoxLayout();
    summary_row->addWidget(new QLabel("Copyable result summary", details));
    summary_row->addStretch(1);
    copy_summary_button_ = new QPushButton("Copy summary", details);
    copy_summary_button_->setEnabled(false);
    summary_row->addWidget(copy_summary_button_);
    details_layout->addLayout(summary_row);
    summary_text_ = new QPlainTextEdit(details);
    summary_text_->setReadOnly(true);
    summary_text_->setMaximumHeight(88);
    summary_text_->setPlaceholderText("Select a result to create a concise summary paragraph.");
    details_layout->addWidget(summary_text_);
    read_details_ = new QPlainTextEdit(details);
    read_details_->setReadOnly(true);
    read_details_->setPlaceholderText("Per-read evidence for the selected result appears here.");
    details_layout->addWidget(read_details_);
    splitter->addWidget(results_);
    splitter->addWidget(details);
    splitter->setSizes({360, 230});
    results_layout->addWidget(splitter);
    tabs_->addTab(results_tab, "VAF results");

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
    broadcast_text_->setPlaceholderText("Received IGV /load and /goto information is appended here. The receiver remains passive.");
    broadcast_layout->addWidget(broadcast_text_);
    tabs_->addTab(broadcast_tab, "Broadcast inbox");
    layout->addWidget(tabs_, 1);
    setCentralWidget(root);

    connect(browse_bam, &QPushButton::clicked, this, [this] { choose_bams(); });
    connect(load_bams_button_, &QPushButton::clicked, this, [this] { load_bams(); });
    connect(clear_bams_button_, &QPushButton::clicked, this, [this] { clear_loaded_bams(); });
    connect(run_button_, &QPushButton::clicked, this, [this] { run_queries(); });
    connect(copy_summary_button_, &QPushButton::clicked, this, [this] { copy_summary(); });
    connect(export_button_, &QPushButton::clicked, this, [this] { export_audit(); });
    connect(results_, &QTableWidget::cellClicked, this, [this](const int row, const int column) { show_result_details(row, column); });
    connect(&watcher_, &QFutureWatcher<MultiBamBatch>::finished, this, [this] { show_results(); });
    connect(&audit_watcher_, &QFutureWatcher<AuditSave>::finished, this, [this] { audit_saved(); });

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
    QStringList bams;
    for (const auto& url : event->mimeData()->urls()) {
        if (url.isLocalFile() && url.toLocalFile().endsWith(".bam", Qt::CaseInsensitive)) bams.append(url.toLocalFile());
    }
    if (bams.isEmpty()) return;
    append_resource_lines(bam_path_, bams);
    status_->setText(QString("%1 BAM path(s) added. Choose Add to analysis.").arg(bams.size()));
    event->acceptProposedAction();
}

void MainWindow::choose_bams() {
    const auto files = QFileDialog::getOpenFileNames(this, "Choose BAM files", {}, "BAM files (*.bam)");
    if (!files.isEmpty()) append_resource_lines(bam_path_, files);
}

void MainWindow::load_bams() {
    if (watcher_.isRunning() || audit_watcher_.isRunning()) return;
    const auto pending_bams = nonempty_lines(bam_path_);
    if (pending_bams.isEmpty()) {
        QMessageBox::information(this, "No pending BAM", "Type, paste, browse, or drop one or more BAM paths or HTTPS URLs first.");
        return;
    }

    QStringList issues;
    int added = 0;
    for (const auto& bam : pending_bams) {
        const bool remote = bam.startsWith("https://");
        const QFileInfo info(bam);
        const QUrl remote_url(bam);
        if ((bam.contains("://") && !remote)
            || (remote && (!remote_url.isValid() || remote_url.host().isEmpty()))
            || (!remote && !bam.endsWith(".bam", Qt::CaseInsensitive))) {
            issues.append(sanitized_resource_uri(bam) + ": expected a local .bam path or valid HTTPS BAM URL.");
            continue;
        }
        if (!remote && (!info.exists() || !info.isFile())) {
            issues.append(bam + ": BAM file not found.");
            continue;
        }
        const auto index = remote ? QString{} : index_path_for(bam);
        if (!remote && index.isEmpty()) {
            issues.append(resource_label(bam) + ": no .bai or .csi index was found beside the BAM.");
            continue;
        }
        if (loaded_bam_paths_.contains(bam)) continue;
        loaded_bam_paths_.append(bam);
        loaded_index_paths_.append(index);
        ++added;
    }
    refresh_loaded_bams();
    if (added > 0) bam_path_->clear();
    status_->setText(QString("%1 BAM(s) active; %2 added. Each BAM will receive a separate VAF result.")
        .arg(loaded_bam_paths_.size()).arg(added));
    if (!issues.isEmpty()) QMessageBox::warning(this, "Some BAMs were not added", issues.join('\n'));
}

void MainWindow::clear_loaded_bams() {
    if (watcher_.isRunning() || audit_watcher_.isRunning()) return;
    loaded_bam_paths_.clear();
    loaded_index_paths_.clear();
    configured_bam_paths_.clear();
    configured_index_paths_.clear();
    result_bam_paths_.clear();
    last_batch_ = {};
    results_->setRowCount(0);
    summary_text_->clear();
    read_details_->clear();
    copy_summary_button_->setEnabled(false);
    refresh_loaded_bams();
    status_->setText("All BAMs cleared.");
}

void MainWindow::refresh_loaded_bams() {
    loaded_bams_->clear();
    for (qsizetype index = 0; index < loaded_bam_paths_.size(); ++index) {
        const auto& bam = loaded_bam_paths_[index];
        const auto& bam_index = loaded_index_paths_[index];
        const auto index_status = bam_index.isEmpty() ? "remote index: automatic" : "index found automatically";
        auto* item = new QListWidgetItem(resource_label(bam) + "  —  " + index_status, loaded_bams_);
        item->setToolTip("BAM: " + sanitized_resource_uri(bam) + "\nIndex: "
            + (bam_index.isEmpty() ? "automatic remote discovery" : bam_index));
    }
    clear_bams_button_->setEnabled(!loaded_bam_paths_.isEmpty() && !watcher_.isRunning());
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
    const auto bams = loaded_bam_paths_;
    const auto indexes = loaded_index_paths_;
    if (bams.isEmpty()) {
        QMessageBox::warning(this, "No active BAM", "Add at least one indexed BAM before calculating VAFs.");
        return;
    }
    if (!minimum_mapq_->hasAcceptableInput() || !minimum_baseq_->hasAcceptableInput()) {
        QMessageBox::warning(this, "Invalid read filters", "MapQ and baseQ must be integers from 0 to 255.");
        return;
    }
    auto parsed = parse_queries(query_text_->toPlainText().toStdString());
    std::vector<Query> variants;
    variants.reserve(parsed.queries.size());
    for (auto& query : parsed.queries) {
        if (std::holds_alternative<VariantQuery>(query)) variants.push_back(std::move(query));
        else parsed.errors.push_back(std::get<RegionQuery>(query).source_text + ": region scans are not part of the VAF workflow; enter a specific REF>ALT variant");
    }
    if (variants.empty()) {
        QMessageBox::warning(this, "No valid variants", QString::fromStdString(parsed.errors.empty() ? "Enter at least one variant." : parsed.errors.front()));
        return;
    }

    last_filters_ = filters();
    configured_bam_paths_ = bams;
    configured_index_paths_ = indexes;
    last_query_text_ = query_text_->toPlainText();
    run_button_->setEnabled(false);
    load_bams_button_->setEnabled(false);
    clear_bams_button_->setEnabled(false);
    status_->setText(QString("Calculating VAFs in %1 BAM(s)…").arg(bams.size()));
    watcher_.setFuture(QtConcurrent::run([queries = std::move(variants), errors = std::move(parsed.errors),
                                          filter_values = last_filters_, bams, indexes] {
        MultiBamBatch combined;
        combined.batch.errors = errors;
        for (qsizetype source_index = 0; source_index < bams.size(); ++source_index) {
            const auto bam_qt = bams[source_index];
            const auto index_qt = indexes[source_index];
            try {
                igv::Resource resource{.uri = bam_qt.toStdString()};
                if (!index_qt.isEmpty()) resource.index_uri = index_qt.toStdString();
                EvidenceEngine engine(std::move(resource));
                auto evaluated = engine.evaluate(queries, filter_values);
                for (auto& result : evaluated.results) {
                    combined.batch.results.push_back(std::move(result));
                    combined.result_bams.append(bam_qt);
                }
                for (auto& error : evaluated.errors) {
                    combined.batch.errors.push_back("[" + resource_label(bam_qt).toStdString() + "] " + error);
                }
            } catch (const std::exception& error) {
                combined.batch.errors.push_back("[" + resource_label(bam_qt).toStdString() + "] "
                    + sanitized_error(error.what(), bam_qt));
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
    for (int index = 0; index < static_cast<int>(last_batch_.results.size()); ++index) {
        const auto* evidence = std::get_if<VariantEvidence>(&last_batch_.results[static_cast<std::size_t>(index)]);
        if (evidence == nullptr) continue;
        results_->insertRow(index);
        const auto& count = evidence->counts;
        const auto bam_label = resource_label(result_bam_paths_[index]);
        const QList<QString> values{bam_label, display_query(evidence->query), count.alternate_molecules > 0 ? "Yes" : "No",
            QString::number(count.alternate_reads), QString::number(count.informative_read_depth()), percent(count.allele_fraction()),
            QString::number(count.alternate_molecules), QString::number(count.molecule_depth()), percent(count.molecule_allele_fraction()),
            QString::number(count.other_reads), QString::number(count.other_molecules), evidence_summary(*evidence, bam_label)};
        for (int column = 0; column < values.size(); ++column) {
            auto* item = new QTableWidgetItem(values[column]);
            if (column == 0) item->setToolTip(sanitized_resource_uri(result_bam_paths_[index]));
            results_->setItem(index, column, item);
        }
    }
    run_button_->setEnabled(true);
    load_bams_button_->setEnabled(true);
    clear_bams_button_->setEnabled(!loaded_bam_paths_.isEmpty());
    status_->setText(QString("%1 VAF result(s) from %2 BAM(s); %3 issue(s).")
        .arg(last_batch_.results.size()).arg(loaded_bam_paths_.size()).arg(last_batch_.errors.size()));
    if (results_->rowCount() > 0) {
        results_->selectRow(0);
        show_result_details(0);
    } else {
        summary_text_->clear();
        copy_summary_button_->setEnabled(false);
    }
    if (!last_batch_.errors.empty()) {
        QString issues = "Issues:\n";
        for (const auto& error : last_batch_.errors) issues += QString::fromStdString(error) + '\n';
        read_details_->setPlainText(issues);
    }
}

void MainWindow::show_result_details(const int row, const int) {
    if (row < 0 || static_cast<std::size_t>(row) >= last_batch_.results.size()) return;
    const auto* evidence = std::get_if<VariantEvidence>(&last_batch_.results[static_cast<std::size_t>(row)]);
    if (evidence == nullptr) return;
    const auto summary = evidence_summary(*evidence, resource_label(result_bam_paths_[row]));
    summary_text_->setPlainText(summary);
    copy_summary_button_->setEnabled(true);
    QString text = "BAM: " + sanitized_resource_uri(result_bam_paths_[row])
        + "\nGrouping: paired fragments by read name (no UMI)\n\n";
    for (const auto& read : evidence->reads) {
        text += QString::fromStdString(read.read_name) + (read.reverse_strand ? "  reverse  " : "  forward  ")
            + QString::fromStdString(read.summary);
        if (!read.molecule_id.empty()) text += "  fragment=" + QString::fromStdString(read.molecule_id);
        text += '\n';
    }
    read_details_->setPlainText(text);
}

void MainWindow::copy_summary() {
    const auto summary = summary_text_->toPlainText().trimmed();
    if (summary.isEmpty()) return;
    QApplication::clipboard()->setText(summary);
    status_->setText("Result summary copied to the clipboard.");
}

void MainWindow::export_audit() {
    if (audit_watcher_.isRunning()) return;
    if (last_batch_.results.empty() && last_batch_.errors.empty()) {
        QMessageBox::information(this, "Nothing to export", "Calculate VAFs before exporting an audit record.");
        return;
    }
    const auto path = QFileDialog::getSaveFileName(this, "Export audit JSON", "bam-seek-audit.json", "JSON files (*.json)");
    if (path.isEmpty()) return;
    export_button_->setEnabled(false);
    load_bams_button_->setEnabled(false);
    clear_bams_button_->setEnabled(false);
    status_->setText("Creating audit fingerprints in the background…");

    QJsonObject root{{"application", "BAM Seek"}, {"version", QStringLiteral(BAM_SEEK_VERSION)},
        {"generated_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs)}, {"query_input", last_query_text_},
        {"molecule_grouping", "paired fragments by read name (no UMI)"},
        {"minimum_mapq", last_filters_.minimum_mapping_quality}, {"minimum_baseq", last_filters_.minimum_base_quality}};
    QJsonArray result_array;
    for (std::size_t result_index = 0; result_index < last_batch_.results.size(); ++result_index) {
        const auto* evidence = std::get_if<VariantEvidence>(&last_batch_.results[result_index]);
        if (evidence == nullptr) continue;
        const auto bam = result_bam_paths_[static_cast<qsizetype>(result_index)];
        const auto& counts = evidence->counts;
        QJsonObject row{{"alignment", sanitized_resource_uri(bam)}, {"query", display_query(evidence->query)},
            {"contig", QString::fromStdString(evidence->query.contig)}, {"position_one_based", evidence->query.position + 1},
            {"reference_allele", QString::fromStdString(evidence->query.reference)},
            {"alternate_allele", QString::fromStdString(evidence->query.alternate)},
            {"ref_reads", counts.reference_reads}, {"alt_reads", counts.alternate_reads}, {"other_reads", counts.other_reads},
            {"informative_read_depth", counts.informative_read_depth()}, {"vaf", counts.allele_fraction()},
            {"ref_molecules", counts.reference_molecules}, {"alt_molecules", counts.alternate_molecules},
            {"ambiguous_molecules", counts.other_molecules}, {"informative_molecule_depth", counts.molecule_depth()},
            {"molecule_vaf", counts.molecule_allele_fraction()},
            {"summary", evidence_summary(*evidence, resource_label(bam))}};
        QJsonArray reads;
        for (const auto& read : evidence->reads) {
            reads.append(QJsonObject{{"name", QString::fromStdString(read.read_name)},
                {"allele", QString::fromStdString(allele_name(read.allele))}, {"reverse", read.reverse_strand},
                {"mapq", read.mapping_quality}, {"minimum_baseq", read.minimum_base_quality},
                {"fragment", QString::fromStdString(read.molecule_id)}});
        }
        row["reads"] = reads;
        result_array.append(row);
    }
    root["results"] = result_array;
    QJsonArray errors;
    for (const auto& error : last_batch_.errors) errors.append(QString::fromStdString(error));
    root["errors"] = errors;

    const auto alignments = configured_bam_paths_;
    const auto indexes = configured_index_paths_;
    audit_watcher_.setFuture(QtConcurrent::run([root = std::move(root), alignments, indexes, path]() mutable {
        AuditSave result{path, {}};
        QJsonArray alignment_array;
        for (qsizetype source_index = 0; source_index < alignments.size(); ++source_index) {
            QJsonObject resource{{"bam", file_identity(alignments[source_index])}};
            resource["index"] = indexes[source_index].isEmpty()
                ? QJsonObject{{"auto_discovery", true}, {"remote", true}}
                : file_identity(indexes[source_index]);
            alignment_array.append(resource);
        }
        root["alignments"] = alignment_array;
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
