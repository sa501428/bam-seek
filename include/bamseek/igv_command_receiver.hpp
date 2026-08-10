#pragma once

#include <QObject>
#include <QString>

class QTcpServer;
class QTcpSocket;

namespace bamseek {

// A deliberately passive implementation of IGV's localhost command protocol.
// Requests are decoded and reported, but never dispatched to application actions.
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

signals:
    void request_received(const QString& description);
    void listening_changed(bool listening, quint16 port, const QString& error);

private:
    void accept_connections();
    void read_client(QTcpSocket* client);

    QTcpServer* server_{};
};

}  // namespace bamseek
