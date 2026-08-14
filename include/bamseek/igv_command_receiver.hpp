#pragma once

#include <QByteArray>
#include <QObject>
#include <QString>
#include <QStringList>

class QTcpServer;
class QTcpSocket;

namespace bamseek {

// A localhost-only transparent relay for IGV's command protocol. Inbound bytes
// are inspected for review and BAM loading, then forwarded unchanged to IGV.
class IgvCommandReceiver final : public QObject {
    Q_OBJECT

public:
    explicit IgvCommandReceiver(QObject* parent = nullptr);

    [[nodiscard]] bool listen(quint16 port, quint16 upstream_port = 60152);
    void close();
    [[nodiscard]] bool is_listening() const;
    [[nodiscard]] quint16 port() const;
    [[nodiscard]] quint16 upstream_port() const noexcept;
    [[nodiscard]] QString error_string() const;

    [[nodiscard]] static QString describe_request(const QString& request_line);
    [[nodiscard]] static QStringList bam_paths_from_request(const QString& request_line);

signals:
    void request_received(const QString& description);
    void bam_load_requested(const QStringList& bam_paths);
    void listening_changed(bool listening, quint16 port, const QString& error);
    void forwarding_changed(bool connected, quint16 port, const QString& error);

private:
    void accept_connections();
    void inspect_client_bytes(QTcpSocket* client, const QByteArray& bytes);
    void fail_client_relay(QTcpSocket* client);

    QTcpServer* server_{};
    quint16 upstream_port_ = 60152;
};

}  // namespace bamseek
