#include <bamseek/igv_command_receiver.hpp>

#include <QCoreApplication>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QHostAddress>
#include <QTcpServer>
#include <QTcpSocket>

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

template <typename Predicate>
bool wait_until(Predicate&& predicate, const int timeout_ms = 2000) {
    QElapsedTimer timer;
    timer.start();
    while (!predicate() && timer.elapsed() < timeout_ms) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 20);
    }
    return predicate();
}

QTcpSocket* accept_connection(QTcpServer& server) {
    if (!wait_until([&server] { return server.hasPendingConnections(); })) return nullptr;
    return server.nextPendingConnection();
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

    QTcpServer mock_igv;
    if (!mock_igv.listen(QHostAddress::LocalHost, 0)) {
        std::cerr << "Skipping loopback relay integration checks: " << mock_igv.errorString().toStdString() << '\n';
        return EXIT_SUCCESS;
    }
    bamseek::IgvCommandReceiver relay;
    if (!relay.listen(0, mock_igv.serverPort()) || relay.upstream_port() != mock_igv.serverPort()) {
        std::cerr << "Could not start BAM Seek relay\n";
        return EXIT_FAILURE;
    }
    QStringList relayed_bams;
    QStringList descriptions;
    QObject::connect(&relay, &bamseek::IgvCommandReceiver::bam_load_requested,
        [&relayed_bams](const QStringList& paths) { relayed_bams.append(paths); });
    QObject::connect(&relay, &bamseek::IgvCommandReceiver::request_received,
        [&descriptions](const QString& description) { descriptions.append(description); });

    QTcpSocket http_client;
    http_client.connectToHost(QHostAddress::LocalHost, relay.port());
    if (!http_client.waitForConnected(2000)) {
        std::cerr << "HTTP client could not connect to relay\n";
        return EXIT_FAILURE;
    }
    const QByteArray http_request =
        "GET /load?file=%2Fdata%2Frelay.bam&genome=hg19&locus=chr1%3A10-20 HTTP/1.1\r\n"
        "Host: localhost:60151\r\nX-Verbatim: a  b\r\nConnection: close\r\n\r\n";
    http_client.write(http_request);
    auto* http_upstream = accept_connection(mock_igv);
    if (http_upstream == nullptr || !wait_until([&] { return http_upstream->bytesAvailable() >= http_request.size(); })) {
        std::cerr << "Mock IGV did not receive HTTP request\n";
        return EXIT_FAILURE;
    }
    const auto received_http = http_upstream->readAll();
    if (received_http != http_request) {
        std::cerr << "HTTP request changed while relaying\n";
        return EXIT_FAILURE;
    }
    const QByteArray http_response =
        "HTTP/1.1 200 OK\r\nContent-Type: text/plain\r\nContent-Length: 7\r\nConnection: close\r\n\r\nOK IGV\n";
    http_upstream->write(http_response);
    http_upstream->disconnectFromHost();
    if (!wait_until([&] { return http_client.bytesAvailable() >= http_response.size(); })
        || http_client.readAll() != http_response) {
        std::cerr << "IGV HTTP response changed while relaying\n";
        return EXIT_FAILURE;
    }
    if (relayed_bams != QStringList{"/data/relay.bam"}
        || descriptions.isEmpty() || !descriptions.front().contains("HTTP GET /load")) {
        std::cerr << "Relayed HTTP request was not inspected correctly\n";
        return EXIT_FAILURE;
    }

    QTcpSocket command_client;
    command_client.connectToHost(QHostAddress::LocalHost, relay.port());
    if (!command_client.waitForConnected(2000)) return EXIT_FAILURE;
    const QByteArray commands = "goto chr2:20-40\nload file=\"/data/raw relay.bam\" genome=hg19\n";
    command_client.write(commands);
    auto* command_upstream = accept_connection(mock_igv);
    if (command_upstream == nullptr || !wait_until([&] { return command_upstream->bytesAvailable() >= commands.size(); })
        || command_upstream->readAll() != commands) {
        std::cerr << "Persistent port commands were not forwarded verbatim\n";
        return EXIT_FAILURE;
    }
    const QByteArray command_response = "OK\nRECEIVED\n";
    command_upstream->write(command_response);
    if (!wait_until([&] { return command_client.bytesAvailable() >= command_response.size(); })
        || command_client.readAll() != command_response) {
        std::cerr << "Persistent IGV responses were not relayed verbatim\n";
        return EXIT_FAILURE;
    }
    if (!relayed_bams.contains("/data/raw relay.bam") || descriptions.size() != 3
        || !descriptions[1].contains("Port command: goto") || !descriptions[2].contains("Port command: load")) {
        std::cerr << "Non-BAM and BAM port commands were not both inspected\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
