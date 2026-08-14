#include <bamseek/igv_command_receiver.hpp>

#include <QCoreApplication>

#include <cstdlib>
#include <iostream>

namespace {

bool contains_all(const QString& text, const QStringList& expected) {
    for (const auto& value : expected) {
        if (!text.contains(value)) {
            std::cerr << "Missing expected receiver text: " << value.toStdString() << "\nActual:\n"
                      << text.toStdString() << '\n';
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {
    QCoreApplication application(argc, argv);

    const auto http = bamseek::IgvCommandReceiver::describe_request(
        "GET /load?file=https%253A%252F%252Fexample.org%252Fa.bam%2Chttps%253A%252F%252Fexample.org%252Fb.bam"
        "&index=https%3A%2F%2Fexample.org%2Fa.bam.bai&locus=chr7%3A140453136 HTTP/1.1");
    if (!contains_all(http, {"HTTP GET /load", "file[1]: https://example.org/a.bam", "file[2]: https://example.org/b.bam",
            "index: https://example.org/a.bam.bai", "locus: chr7:140453136"})) {
        return EXIT_FAILURE;
    }

    const auto raw = bamseek::IgvCommandReceiver::describe_request(
        "load \"/data/sample one.bam\" index=\"/data/sample one.bam.bai\" locus=chr1:10-20");
    if (!contains_all(raw, {"Port command: load", "file: /data/sample one.bam",
            "index: /data/sample one.bam.bai", "locus: chr1:10-20"})) {
        return EXIT_FAILURE;
    }

    const auto intercepted = bamseek::IgvCommandReceiver::bam_paths_from_request(
        "GET /load?file=/path/to/XYZ_consensus_filtered.bam&genome=hg19"
        "&locus=chr13:28608000-28608400&merge=false HTTP/1.1");
    if (intercepted != QStringList{"/path/to/XYZ_consensus_filtered.bam"}) {
        std::cerr << "Receiver did not isolate the /load file value\n";
        return EXIT_FAILURE;
    }
    if (!bamseek::IgvCommandReceiver::bam_paths_from_request(
            "GET /goto?file=%2Fpath%2Fto%2Fignored.bam&locus=chr13%3A1-2 HTTP/1.1").isEmpty()
        || !bamseek::IgvCommandReceiver::bam_paths_from_request(
            "HEAD /load?file=%2Fpath%2Fto%2Fignored.bam HTTP/1.1").isEmpty()) {
        std::cerr << "Only actionable GET /load requests may emit BAM paths\n";
        return EXIT_FAILURE;
    }
    const auto file_uri = bamseek::IgvCommandReceiver::bam_paths_from_request(
        "load file=file:///data/sample%20one.bam genome=hg19 locus=chr1:10-20");
    if (file_uri != QStringList{"/data/sample one.bam"}) {
        std::cerr << "Receiver did not normalize a local file URI\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
