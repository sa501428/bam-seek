#include <bamseek/evidence.hpp>
#include <bamseek/pileup_view.hpp>

#include <QApplication>
#include <QImage>
#include <QPainter>
#include <QScrollArea>
#include <QScrollBar>

#include <cstdlib>
#include <iostream>

int main(int argc, char* argv[]) {
    QApplication application(argc, argv);
    bamseek::FilterSettings filters;
    filters.minimum_mapping_quality = 20;
    filters.minimum_base_quality = 20;
    filters.molecule_mode = bamseek::MoleculeMode::auto_detect;
    const bamseek::VariantQuery query{"JAK2 V617F", "chr9", 5073769, "G", "T"};
    bamseek::EvidenceEngine engine({.uri = BAM_SEEK_JAK2_PILEUP_BAM});
    const auto data = engine.pileup(query, filters, 40);
    if (data.alignments.size() != 150 || data.total_alignments != 150) return EXIT_FAILURE;

    QScrollArea scroll;
    scroll.resize(800, 600);
    scroll.setWidgetResizable(false);
    auto* view = new bamseek::PileupView(&scroll);
    scroll.setWidget(view);
    view->set_data(data);
    scroll.show();
    application.processEvents();
    if (view->width() < 900 || view->height() < 1000 || view->variant_x() <= 0 || view->variant_x() >= view->width()
        || scroll.horizontalScrollBar()->maximum() <= 0 || scroll.verticalScrollBar()->maximum() <= 0) {
        std::cerr << "Pileup widget did not acquire a scrollable content size\n";
        return EXIT_FAILURE;
    }

    QImage rendered(view->size(), QImage::Format_RGB32);
    rendered.fill(Qt::white);
    QPainter painter(&rendered);
    view->render(&painter);
    painter.end();
    std::size_t non_white = 0;
    for (int y = 0; y < rendered.height(); y += 5) {
        for (int x = 0; x < rendered.width(); x += 5) {
            if (rendered.pixelColor(x, y) != QColor(Qt::white)) ++non_white;
        }
    }
    if (non_white < 1000) {
        std::cerr << "Pileup rendering was unexpectedly blank\n";
        return EXIT_FAILURE;
    }
    const auto screenshot = qEnvironmentVariable("BAM_SEEK_PILEUP_SCREENSHOT");
    if (!screenshot.isEmpty() && !rendered.save(screenshot)) {
        std::cerr << "Could not save requested pileup screenshot\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
