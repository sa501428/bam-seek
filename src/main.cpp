#include <bamseek/main_window.hpp>

#include <QApplication>

int main(int argc, char* argv[]) {
    QApplication application(argc, argv);
    application.setApplicationName("BAM Seek");
    application.setOrganizationName("BAM Seek");
    bamseek::MainWindow window;
    window.show();
    return application.exec();
}
