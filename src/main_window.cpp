#include <bamseek/main_window.hpp>

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
#include <QMessageBox>
#include <QMimeData>
#include <QIODevice>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QRegularExpression>
#include <QRegularExpressionValidator>
#include <QSaveFile>
#include <QScrollArea>
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
    return QString::fromStdString(query.contig) + ':' + QString::number(query.position + 1) + ' '
        + QString::fromStdString(query.reference) + '>' + QString::fromStdString(query.alternate);
}

QString mode_name(const MoleculeMode mode) {
    switch (mode) {
        case MoleculeMode::raw_reads: return "Raw reads";
        case MoleculeMode::auto_detect: return "Auto-detect (MI, RX, UB)";
        case MoleculeMode::selected_tag: return "Selected tag";
    }
    return "Unknown";
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

}  // namespace

MainWindow::MainWindow() {
    setWindowTitle("BAM Seek — targeted BAM evidence");
    setAcceptDrops(true);
    resize(1240, 820);

    auto* root = new QWidget(this);
    auto* layout = new QVBoxLayout(root);
    auto* input_form = new QFormLayout();
    auto* bam_row = new QHBoxLayout();
    bam_path_ = new QLineEdit(root);
    bam_path_->setPlaceholderText("Drop a BAM/CRAM, enter a local path, or use an HTTPS URL");
    auto* browse_bam = new QPushButton("Browse…", root);
    bam_row->addWidget(bam_path_);
    bam_row->addWidget(browse_bam);
    input_form->addRow("BAM / CRAM", bam_row);
    auto* index_row = new QHBoxLayout();
    index_path_ = new QLineEdit(root);
    index_path_->setPlaceholderText("Optional explicit .bai/.csi/.crai path or HTTPS URL");
    auto* browse_index = new QPushButton("Browse…", root);
    index_row->addWidget(index_path_);
    index_row->addWidget(browse_index);
    input_form->addRow("Index", index_row);
    auto* reference_row = new QHBoxLayout();
    reference_path_ = new QLineEdit(root);
    reference_path_->setPlaceholderText("Optional for BAM; required for CRAM (hg19 FASTA)");
    auto* browse_reference = new QPushButton("Browse…", root);
    reference_row->addWidget(reference_path_);
    reference_row->addWidget(browse_reference);
    input_form->addRow("Reference", reference_row);
    layout->addLayout(input_form);

    auto* settings = new QHBoxLayout();
    molecule_mode_ = new QComboBox(root);
    molecule_mode_->addItem("Auto-detect (MI, RX, UB)", static_cast<int>(MoleculeMode::auto_detect));
    molecule_mode_->addItem("Raw reads only", static_cast<int>(MoleculeMode::raw_reads));
    molecule_mode_->addItem("Selected BAM tag", static_cast<int>(MoleculeMode::selected_tag));
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
    query_text_->setPlainText("# Variants use one-based genomic coordinates\n# chr7:140453136 A>T");
    layout->addWidget(query_text_, 1);

    auto* buttons = new QHBoxLayout();
    run_button_ = new QPushButton("Query evidence", root);
    pileup_button_ = new QPushButton("View pileup", root);
    export_button_ = new QPushButton("Export audit JSON…", root);
    status_ = new QLabel("Select an indexed BAM or CRAM, then enter queries.", root);
    buttons->addWidget(run_button_); buttons->addWidget(pileup_button_); buttons->addWidget(export_button_); buttons->addWidget(status_, 1);
    layout->addLayout(buttons);

    tabs_ = new QTabWidget(root);
    auto* evidence_tab = new QWidget(tabs_);
    auto* evidence_layout = new QVBoxLayout(evidence_tab);
    evidence_layout->setContentsMargins(0, 0, 0, 0);
    auto* splitter = new QSplitter(Qt::Vertical, evidence_tab);
    results_ = new QTableWidget(splitter);
    results_->setColumnCount(12);
    results_->setHorizontalHeaderLabels({"Query", "Status", "Depth", "Alt reads", "VAF", "Alt Fwd", "Alt Rev", "Strand P", "Alt molecules", "Ref molecules", "Molecule tag", "Notes"});
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
    auto* pileup_scroll = new QScrollArea(pileup_tab);
    pileup_view_ = new PileupView(pileup_scroll);
    pileup_scroll->setWidget(pileup_view_);
    pileup_scroll->setWidgetResizable(false);
    pileup_layout->addWidget(pileup_scroll);
    tabs_->addTab(pileup_tab, "Pileup");
    layout->addWidget(tabs_, 2);
    setCentralWidget(root);

    connect(browse_bam, &QPushButton::clicked, this, [this] { choose_bam(); });
    connect(browse_index, &QPushButton::clicked, this, [this] { choose_index(); });
    connect(browse_reference, &QPushButton::clicked, this, [this] { choose_reference(); });
    connect(run_button_, &QPushButton::clicked, this, [this] { run_queries(); });
    connect(pileup_button_, &QPushButton::clicked, this, [this] { show_pileup(); });
    connect(export_button_, &QPushButton::clicked, this, [this] { export_audit(); });
    connect(molecule_mode_, &QComboBox::currentIndexChanged, this, [this] {
        molecule_tag_->setEnabled(molecule_mode_->currentData().toInt() == static_cast<int>(MoleculeMode::selected_tag));
    });
    connect(results_, &QTableWidget::cellClicked, this, [this](const int row, const int column) { show_read_details(row, column); });
    connect(&watcher_, &QFutureWatcher<BatchEvidence>::finished, this, [this] { show_results(); });
    connect(&pileup_watcher_, &QFutureWatcher<PileupLoad>::finished, this, [this] { pileup_loaded(); });
    connect(&audit_watcher_, &QFutureWatcher<AuditSave>::finished, this, [this] { audit_saved(); });
    connect(group_pairs_, &QCheckBox::toggled, this, [this](const bool enabled) { pileup_view_->set_group_pairs(enabled); });
}

void MainWindow::dragEnterEvent(QDragEnterEvent* event) {
    if (event->mimeData()->hasUrls()) event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent* event) {
    const auto urls = event->mimeData()->urls();
    if (urls.isEmpty() || !urls.front().isLocalFile()) return;
    bam_path_->setText(urls.front().toLocalFile());
    status_->setText("BAM/CRAM selected. An accompanying index is required.");
    event->acceptProposedAction();
}

void MainWindow::choose_bam() {
    const auto file = QFileDialog::getOpenFileName(this, "Choose indexed BAM or CRAM", {}, "Alignment files (*.bam *.cram *.sam);;All files (*)");
    if (!file.isEmpty()) bam_path_->setText(file);
}

void MainWindow::choose_index() {
    const auto file = QFileDialog::getOpenFileName(this, "Choose alignment index", {}, "Alignment indexes (*.bai *.csi *.crai);;All files (*)");
    if (!file.isEmpty()) index_path_->setText(file);
}

void MainWindow::choose_reference() {
    const auto file = QFileDialog::getOpenFileName(this, "Choose hg19 reference FASTA", {}, "FASTA files (*.fa *.fasta *.fna);;All files (*)");
    if (!file.isEmpty()) reference_path_->setText(file);
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
    if (bam_path_->text().trimmed().isEmpty()) {
        QMessageBox::warning(this, "No alignment file", "Choose or drop an indexed BAM or CRAM first.");
        return;
    }
    if (bam_path_->text().startsWith("http://") || index_path_->text().startsWith("http://") || reference_path_->text().startsWith("http://")) {
        QMessageBox::warning(this, "Secure remote access", "BAM Seek accepts local paths or HTTPS resources only. HTTP URLs are not permitted.");
        return;
    }
    if (!vaf_->hasAcceptableInput() || !minimum_alt_reads_->hasAcceptableInput() || !minimum_alt_molecules_->hasAcceptableInput()
        || !minimum_mapq_->hasAcceptableInput() || !minimum_baseq_->hasAcceptableInput()
        || (molecule_mode_->currentData().toInt() == static_cast<int>(MoleculeMode::selected_tag) && !molecule_tag_->hasAcceptableInput())) {
        QMessageBox::warning(this, "Invalid filters", "Correct the highlighted numeric filters and enter a valid two-character BAM tag.");
        return;
    }
    const auto parsed = parse_queries(query_text_->toPlainText().toStdString());
    if (parsed.queries.empty()) {
        QMessageBox::warning(this, "No valid queries", QString::fromStdString(parsed.errors.empty() ? "Enter at least one query." : parsed.errors.front()));
        return;
    }
    last_filters_ = filters();
    last_bam_path_ = bam_path_->text().trimmed();
    last_index_path_ = index_path_->text().trimmed();
    last_reference_path_ = reference_path_->text().trimmed();
    last_query_text_ = query_text_->toPlainText();
    const auto bam = last_bam_path_.toStdString();
    const auto index = last_index_path_.toStdString();
    const auto reference = last_reference_path_.toStdString();
    run_button_->setEnabled(false);
    status_->setText("Querying indexed alignment evidence…");
    watcher_.setFuture(QtConcurrent::run([queries = parsed.queries, errors = parsed.errors, filter_values = last_filters_, bam, index, reference] {
        BatchEvidence batch;
        batch.errors = errors;
        try {
            igv::Resource resource{.uri = bam};
            if (!index.empty()) resource.index_uri = index;
            if (!reference.empty()) resource.reference_uri = reference;
            EvidenceEngine engine(std::move(resource));
            auto evaluated = engine.evaluate(queries, filter_values);
            batch.results = std::move(evaluated.results);
            batch.errors.insert(batch.errors.end(), evaluated.errors.begin(), evaluated.errors.end());
        } catch (const std::exception& error) {
            batch.errors.push_back(sanitized_error(error.what(), {bam, index, reference}));
        }
        return batch;
    }));
}

void MainWindow::show_results() {
    last_batch_ = watcher_.result();
    results_->setRowCount(0);
    for (int index = 0; index < static_cast<int>(last_batch_.results.size()); ++index) {
        results_->insertRow(index);
        const auto& result = last_batch_.results[static_cast<std::size_t>(index)];
        if (const auto* evidence = std::get_if<VariantEvidence>(&result)) {
            const auto& count = evidence->counts;
            QString notes = QString("%1 callable reads; threshold %2")
                .arg(evidence->reads.size()).arg(evidence->passes_thresholds ? "met" : "not met");
            if (evidence->reads_missing_molecule_tag > 0) notes += QString("; %1 missing molecule tag").arg(evidence->reads_missing_molecule_tag);
            const auto strand_p = count.strand_bias_p_value();
            const QList<QString> values{display_query(evidence->query), evidence->passes_thresholds ? "PRESENT" : "not detected", QString::number(count.depth()),
                QString::number(count.alternate_reads), QString::number(count.allele_fraction() * 100.0, 'f', 4) + '%', QString::number(count.alternate_forward_reads),
                QString::number(count.alternate_reverse_reads), strand_p ? QString::number(*strand_p, 'g', 3) : "N/A",
                evidence->molecule_counts_available ? QString::number(count.alternate_molecules) : "N/A",
                evidence->molecule_counts_available ? QString::number(count.reference_molecules) : "N/A",
                QString::fromStdString(evidence->molecule_counts_available ? evidence->molecule_tag_used : "Unavailable"), notes};
            for (int column = 0; column < values.size(); ++column) results_->setItem(index, column, new QTableWidgetItem(values[column]));
        } else {
            const auto& region = std::get<RegionEvidence>(result);
            const QList<QString> values{QString::fromStdString(region.query.source_text), region.candidates.empty() ? "NO CANDIDATES" : "CANDIDATES", "", QString::number(region.candidates.size()), "", "", "", "", "", "", "", QString::fromStdString(region.note)};
            for (int column = 0; column < values.size(); ++column) results_->setItem(index, column, new QTableWidgetItem(values[column]));
        }
    }
    run_button_->setEnabled(true);
    status_->setText(QString("%1 result(s), %2 issue(s). Select a row for read evidence.").arg(last_batch_.results.size()).arg(last_batch_.errors.size()));
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
        QString text = QString::fromStdString(region.note) + "\n\n";
        for (const auto& candidate : region.candidates) {
            text += display_query(candidate.query) + "  depth=" + QString::number(candidate.counts.depth())
                + "  alt=" + QString::number(candidate.counts.alternate_reads)
                + "  VAF=" + QString::number(candidate.counts.allele_fraction() * 100.0, 'f', 4) + "%\n";
        }
        read_details_->setPlainText(text);
        return;
    }
    QString text = "Query: " + display_query(evidence->query) + "\nMolecule grouping: "
        + QString::fromStdString(evidence->molecule_counts_available ? evidence->molecule_tag_used : "unavailable") + "\n\n";
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
    const auto bam = last_bam_path_.toStdString();
    const auto index = last_index_path_.toStdString();
    const auto reference = last_reference_path_.toStdString();
    const auto query = evidence->query;
    const auto filter_values = last_filters_;
    pileup_button_->setEnabled(false);
    status_->setText("Loading local alignment pileup…");
    pileup_watcher_.setFuture(QtConcurrent::run([bam, index, reference, query, filter_values] {
        PileupLoad loaded;
        try {
            igv::Resource resource{.uri = bam};
            if (!index.empty()) resource.index_uri = index;
            if (!reference.empty()) resource.reference_uri = reference;
            EvidenceEngine engine(std::move(resource));
            loaded.data = engine.pileup(query, filter_values);
        } catch (const std::exception& error) {
            loaded.error = sanitized_error(error.what(), {bam, index, reference});
        }
        return loaded;
    }));
}

void MainWindow::pileup_loaded() {
    const auto loaded = pileup_watcher_.result();
    pileup_button_->setEnabled(true);
    if (!loaded.error.empty()) {
        QMessageBox::warning(this, "Pileup unavailable", QString::fromStdString(loaded.error));
        status_->setText("Could not load pileup.");
        return;
    }
    pileup_view_->set_data(loaded.data);
    tabs_->setCurrentIndex(1);
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
        QJsonObject row{{"query", display_query(evidence.query)}, {"present", evidence.passes_thresholds}, {"depth", evidence.counts.depth()},
            {"ref_reads", evidence.counts.reference_reads}, {"other_reads", evidence.counts.other_reads},
            {"ref_forward_reads", evidence.counts.reference_forward_reads}, {"ref_reverse_reads", evidence.counts.reference_reverse_reads},
            {"alt_reads", evidence.counts.alternate_reads}, {"alt_forward_reads", evidence.counts.alternate_forward_reads},
            {"alt_reverse_reads", evidence.counts.alternate_reverse_reads}, {"vaf", evidence.counts.allele_fraction()},
            {"molecule_counts_available", evidence.molecule_counts_available}, {"reads_missing_molecule_tag", evidence.reads_missing_molecule_tag},
            {"molecule_tag", QString::fromStdString(evidence.molecule_tag_used)}};
        row["alt_molecules"] = evidence.molecule_counts_available ? QJsonValue(evidence.counts.alternate_molecules) : QJsonValue::Null;
        row["ref_molecules"] = evidence.molecule_counts_available ? QJsonValue(evidence.counts.reference_molecules) : QJsonValue::Null;
        const auto strand_p = evidence.counts.strand_bias_p_value();
        row["strand_bias_fisher_p"] = strand_p ? QJsonValue(*strand_p) : QJsonValue::Null;
        QJsonArray reads;
        for (const auto& read : evidence.reads) reads.append(QJsonObject{{"name", QString::fromStdString(read.read_name)}, {"allele", QString::fromStdString(allele_name(read.allele))}, {"reverse", read.reverse_strand}, {"mapq", read.mapping_quality}, {"minimum_baseq", read.minimum_base_quality}, {"molecule", QString::fromStdString(read.molecule_id)}, {"summary", QString::fromStdString(read.summary)}});
        row["reads"] = reads;
        return row;
    };
    for (const auto& item : last_batch_.results) {
        if (const auto* evidence = std::get_if<VariantEvidence>(&item)) {
            result_array.append(variant_json(*evidence));
        } else {
            const auto& region = std::get<RegionEvidence>(item);
            QJsonObject row{{"region", QString::fromStdString(region.query.source_text)}, {"note", QString::fromStdString(region.note)}};
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
    const auto alignment = last_bam_path_;
    auto index = last_index_path_;
    if (index.isEmpty() && !alignment.startsWith("https://")) index = index_path_for(alignment);
    const auto reference = last_reference_path_;
    audit_watcher_.setFuture(QtConcurrent::run([root = std::move(root), alignment, index, reference, path]() mutable {
        AuditSave result{path, {}};
        root["alignment"] = file_identity(alignment);
        root["index"] = index.isEmpty() ? QJsonObject{{"auto_detected", true}} : file_identity(index);
        root["reference"] = file_identity(reference);
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
    if (!result.error.empty()) {
        QMessageBox::critical(this, "Export failed", QString::fromStdString(result.error));
        status_->setText("Audit export failed.");
        return;
    }
    status_->setText("Audit export saved: " + result.path);
}

}  // namespace bamseek
