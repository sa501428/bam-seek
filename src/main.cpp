#include <bamseek/main_window.hpp>

#include <QApplication>
#include <QIcon>
#include <QPainter>
#include <QPixmap>
#include <QSvgRenderer>

namespace {

QIcon applicationIcon() {
    QSvgRenderer renderer(QStringLiteral(":/logo.svg"));
    if (!renderer.isValid()) {
        return {};
    }

    QPixmap pixmap(512, 512);
    pixmap.fill(Qt::transparent);
    QPainter painter(&pixmap);
    renderer.render(&painter);
    return QIcon(pixmap);
}

}  // namespace

int main(int argc, char* argv[]) {
    QApplication application(argc, argv);
    application.setApplicationName("BAM Seek");
    application.setOrganizationName("BAM Seek");
    application.setWindowIcon(applicationIcon());
    bamseek::MainWindow window;
    window.show();
    return application.exec();
}
