#pragma once

#include <QObject>
#include <QString>
#include <QStringList>

class QTcpServer;
class QTcpSocket;

namespace bamseek {

// A localhost-only implementation of IGV's command protocol. Requests are
// decoded for review, and /load file values are emitted separately so the
// application can prepare the requested BAM without acting on other fields.
class IgvCommandReceiver final : public QObject {
    Q_OBJECT

public:
    explicit IgvCommandReceiver(QObject* parent = nullptr);

    [[nodiscard]] bool listen(quint16 port);
    void close();
    [[nodiscard]] bool is_listening() const;
    [[nodiscard]] quint16 port() const;
    [[nodiscard]] QString error_string() const;

    [[nodiscard]] static QString describe_request(const QString& request_line);
    [[nodiscard]] static QStringList bam_paths_from_request(const QString& request_line);

signals:
    void request_received(const QString& description);
    void bam_load_requested(const QStringList& bam_paths);
    void listening_changed(bool listening, quint16 port, const QString& error);

private:
    void accept_connections();
    void read_client(QTcpSocket* client);

    QTcpServer* server_{};
};

}  // namespace bamseek
