#include <bamseek/igv_command_receiver.hpp>

#include <QHostAddress>
#include <QProcess>
#include <QTcpServer>
#include <QTcpSocket>
#include <QUrl>
#include <QUrlQuery>

#include <utility>

namespace bamseek {
namespace {

QString decoded(QString value) {
    // QUrlQuery performs normal URL decoding. A second layer is common in IGV
    // links because its command executor conditionally decodes file arguments.
    for (int pass = 0; pass < 2 && value.contains('%'); ++pass) {
        const auto next = QUrl::fromPercentEncoding(value.toUtf8());
        if (next == value) break;
        value = next;
    }
    return value;
}

void append_value(QStringList& output, const QString& label, const QString& value) {
    if (value.isEmpty()) return;
    output.append(label + ": " + decoded(value));
}

QString last_query_value(const QUrlQuery& query, const QString& key) {
    const auto values = query.allQueryItemValues(key, QUrl::FullyDecoded);
    return values.isEmpty() ? QString{} : values.back();
}

QString normalized_file(QString value) {
    value = decoded(std::move(value)).trimmed();
    const QUrl url(value);
    if (url.isLocalFile()) return url.toLocalFile();
    return value;
}

QStringList split_file_paths(const QString& value) {
    QStringList paths;
    for (const auto& part : decoded(value).split(',', Qt::SkipEmptyParts)) {
        const auto path = normalized_file(part);
        if (!path.isEmpty() && !paths.contains(path)) paths.append(path);
    }
    return paths;
}

QStringList load_files_from_http(const QString& request_line) {
    const auto fields = request_line.simplified().split(' ');
    if (fields.size() < 2 || fields.front().compare("GET", Qt::CaseInsensitive) != 0) return {};
    const QUrl target(fields[1]);
    if (target.path().compare("/load", Qt::CaseInsensitive) != 0) return {};
    return split_file_paths(last_query_value(QUrlQuery(target), "file"));
}

QStringList load_files_from_port_command(const QString& request_line) {
    const auto arguments = QProcess::splitCommand(request_line);
    if (arguments.isEmpty() || arguments.front().compare("load", Qt::CaseInsensitive) != 0) return {};
    QStringList paths;
    for (int index = 1; index < arguments.size(); ++index) {
        const auto& argument = arguments[index];
        const auto separator = argument.indexOf('=');
        if (separator > 0 && argument.left(separator).compare("file", Qt::CaseInsensitive) != 0) continue;
        const auto value = separator > 0 ? argument.mid(separator + 1) : argument;
        for (const auto& path : split_file_paths(value)) if (!paths.contains(path)) paths.append(path);
    }
    return paths;
}

void append_list_value(QStringList& output, const QString& label, const QString& value) {
    const auto values = decoded(value).split(',', Qt::SkipEmptyParts);
    if (values.size() <= 1) {
        append_value(output, label, value);
        return;
    }
    for (qsizetype index = 0; index < values.size(); ++index) {
        output.append(QString("%1[%2]: %3").arg(label).arg(index + 1).arg(values[index].trimmed()));
    }
}

QString describe_http(const QString& request_line) {
    const auto fields = request_line.simplified().split(' ');
    if (fields.size() < 2) return "HTTP request\nRaw: " + request_line;

    const QUrl target(fields[1]);
    const QUrlQuery query(target);
    const auto path = target.path().isEmpty() ? QStringLiteral("/") : target.path();
    QStringList output{"HTTP " + fields.front() + ' ' + path};

    const QStringList file_keys{"file", "bigDataURL", "sessionURL", "dataURL", "hubURL"};
    for (const auto& key : file_keys) {
        const auto value = last_query_value(query, key);
        if (!value.isEmpty()) append_list_value(output, key, value);
    }
    append_list_value(output, "index", last_query_value(query, "index"));
    append_value(output, "locus", last_query_value(query, "locus"));
    append_value(output, "genome", last_query_value(query, "genome"));
    append_list_value(output, "name", last_query_value(query, "name"));
    append_value(output, "command", last_query_value(query, "command"));
    output.append("Raw: " + request_line);
    return output.join('\n');
}

QString describe_port_command(const QString& request_line) {
    const auto arguments = QProcess::splitCommand(request_line);
    if (arguments.isEmpty()) return "Empty command";

    QStringList output{"Port command: " + arguments.front()};
    const auto command = arguments.front().toLower();
    if (command == "goto" || command == "reload") {
        append_value(output, "locus", arguments.mid(1).join(' '));
    } else if (command == "load") {
        for (int index = 1; index < arguments.size(); ++index) {
            const auto& argument = arguments[index];
            const auto separator = argument.indexOf('=');
            if (separator > 0) {
                const auto key = argument.left(separator);
                if (key.compare("file", Qt::CaseInsensitive) == 0 || key.compare("index", Qt::CaseInsensitive) == 0) {
                    append_list_value(output, key, argument.mid(separator + 1));
                } else {
                    append_value(output, key, argument.mid(separator + 1));
                }
            } else {
                append_list_value(output, "file", argument);
            }
        }
    } else if (command == "execute") {
        append_value(output, "command", arguments.mid(1).join(' '));
    } else {
        append_value(output, "arguments", arguments.mid(1).join(' '));
    }
    output.append("Raw: " + request_line);
    return output.join('\n');
}

QByteArray http_response(const bool options, const bool head) {
    if (options) {
        return "HTTP/1.1 204 No Content\r\n"
               "Access-Control-Allow-Origin: *\r\n"
               "Access-Control-Allow-Headers: access-control-allow-origin\r\n"
               "Access-Control-Allow-Private-Network: true\r\n"
               "Access-Control-Allow-Methods: HEAD, GET, OPTIONS\r\n"
               "Connection: close\r\n\r\n";
    }
    static const QByteArray body = "RECEIVED\n";
    auto response = "HTTP/1.1 202 Accepted\r\n"
           "Access-Control-Allow-Origin: *\r\n"
           "Access-Control-Allow-Private-Network: true\r\n"
           "Cache-Control: no-cache, no-store\r\n"
           "Content-Type: text/plain; charset=utf-8\r\n"
           "Content-Length: " + QByteArray::number(body.size()) + "\r\n"
           "Connection: close\r\n\r\n";
    if (!head) response += body;
    return response;
}

}  // namespace

IgvCommandReceiver::IgvCommandReceiver(QObject* parent) : QObject(parent), server_(new QTcpServer(this)) {
    connect(server_, &QTcpServer::newConnection, this, &IgvCommandReceiver::accept_connections);
}

bool IgvCommandReceiver::listen(const quint16 port) {
    close();
    // Match IGV's local command service while preventing access from other hosts.
    const auto started = server_->listen(QHostAddress::LocalHost, port);
    emit listening_changed(started, started ? server_->serverPort() : port, started ? QString{} : server_->errorString());
    return started;
}

void IgvCommandReceiver::close() {
    if (!server_->isListening()) return;
    const auto old_port = server_->serverPort();
    server_->close();
    emit listening_changed(false, old_port, {});
}

bool IgvCommandReceiver::is_listening() const { return server_->isListening(); }

quint16 IgvCommandReceiver::port() const { return server_->serverPort(); }

QString IgvCommandReceiver::error_string() const { return server_->errorString(); }

QString IgvCommandReceiver::describe_request(const QString& request_line) {
    const auto trimmed = request_line.trimmed();
    if (trimmed.startsWith("GET ", Qt::CaseInsensitive) || trimmed.startsWith("HEAD ", Qt::CaseInsensitive)
        || trimmed.startsWith("OPTIONS ", Qt::CaseInsensitive)) {
        return describe_http(trimmed);
    }
    return describe_port_command(trimmed);
}

QStringList IgvCommandReceiver::bam_paths_from_request(const QString& request_line) {
    const auto trimmed = request_line.trimmed();
    if (trimmed.startsWith("GET ", Qt::CaseInsensitive)) return load_files_from_http(trimmed);
    return load_files_from_port_command(trimmed);
}

void IgvCommandReceiver::accept_connections() {
    while (server_->hasPendingConnections()) {
        auto* client = server_->nextPendingConnection();
        connect(client, &QTcpSocket::readyRead, this, [this, client] { read_client(client); });
        connect(client, &QTcpSocket::disconnected, client, &QObject::deleteLater);
    }
}

void IgvCommandReceiver::read_client(QTcpSocket* client) {
    auto buffer = client->property("request_buffer").toByteArray();
    buffer += client->readAll();
    if (buffer.size() > 64 * 1024) {
        client->write("HTTP/1.1 413 Content Too Large\r\nConnection: close\r\n\r\n");
        client->disconnectFromHost();
        return;
    }
    const auto first_newline = buffer.indexOf('\n');
    if (first_newline < 0) {
        client->setProperty("request_buffer", buffer.left(64 * 1024));
        return;
    }

    const auto first_line = QString::fromUtf8(buffer.left(first_newline)).trimmed();
    const auto is_http = first_line.startsWith("GET ", Qt::CaseInsensitive)
        || first_line.startsWith("HEAD ", Qt::CaseInsensitive)
        || first_line.startsWith("OPTIONS ", Qt::CaseInsensitive);
    if (is_http) {
        if (!buffer.contains("\r\n\r\n") && !buffer.contains("\n\n")) {
            client->setProperty("request_buffer", buffer.left(64 * 1024));
            return;
        }
        if (!first_line.startsWith("OPTIONS ", Qt::CaseInsensitive)) {
            emit request_received(describe_request(first_line));
            const auto bams = bam_paths_from_request(first_line);
            if (!bams.isEmpty()) emit bam_load_requested(bams);
        }
        client->write(http_response(first_line.startsWith("OPTIONS ", Qt::CaseInsensitive),
            first_line.startsWith("HEAD ", Qt::CaseInsensitive)));
        client->disconnectFromHost();
        return;
    }

    while (true) {
        const auto newline = buffer.indexOf('\n');
        if (newline < 0) break;
        const auto line = QString::fromUtf8(buffer.left(newline)).trimmed();
        buffer.remove(0, newline + 1);
        if (line.isEmpty()) continue;
        emit request_received(describe_request(line));
        const auto bams = bam_paths_from_request(line);
        if (!bams.isEmpty()) emit bam_load_requested(bams);
        client->write("RECEIVED\n");
    }
    client->setProperty("request_buffer", buffer.left(64 * 1024));
}

}  // namespace bamseek
