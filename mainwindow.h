#pragma once
#include <QMainWindow>
#include <QWidget>
#include <QPointF>
#include <QString>
#include <QPainter>
#include <cmath>
#include <vector>
#include <algorithm>

class QWheelEvent;
class QMouseEvent;
class QPaintEvent;
class QDoubleSpinBox;
class QTableWidget;
class QPushButton;
class QGroupBox;
class QComboBox;

struct Circle
{
    double x = 0, y = 0, r = 1;
};

struct TangentCircle
{
    Circle c;
    int s1 = 1, s2 = 1, s3 = 1;
};

enum class DrawMode
{
    TTT, TTR
};

inline bool tangent(const Circle& sol, const Circle& known, int s, double eps = 1e-6)
{
    double d = std::hypot(sol.x - known.x, sol.y - known.y);
    double expected = (s > 0) ? (sol.r + known.r) : std::abs(sol.r - known.r);
    double scale = std::max(1.0, std::max({sol.r, known.r, d}));
    return std::abs(d - expected) < eps * scale;
}

inline bool solution(const Circle& c1, const Circle& c2, const Circle& c3, const TangentCircle& sol, double eps = 1e-6)
{
    return tangent(sol.c, c1, sol.s1, eps)&& tangent(sol.c, c2, sol.s2, eps)&& tangent(sol.c, c3, sol.s3, eps);
}

class CanvasWidget : public QWidget
{
    Q_OBJECT

public:
    explicit CanvasWidget(QWidget* parent = nullptr);
    void setScene(const Circle& c1, const Circle& c2, const Circle& c3, const std::vector<TangentCircle>& sols, DrawMode mode);
    void renderScene(QPainter* p) const;

    void resetView();
    double zoom() const { return m_zoom; }
    QPointF pan() const { return m_pan; }
    void setView(double zoom, const QPointF& pan);

protected:
    void paintEvent(QPaintEvent*) override;
    void wheelEvent(QWheelEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;
    void mouseDoubleClickEvent(QMouseEvent*) override;

private:
    Circle m_c1{0,0,1}, m_c2{4,0,1}, m_c3{2,3,1};
    std::vector<TangentCircle> m_sols;
    QRectF worldBounds() const;
    DrawMode m_mode = DrawMode::TTT;
    double  m_zoom = 1.0;
    QPointF m_pan{0, 0};
    bool    m_dragging = false;
    QPoint  m_lastMousePos;
};

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);

private slots:
    void onSolve();
    void onExportSvg();
    void onExportCsv();
    void onVerifyExternal();
    void onModeChanged(int index);

private:
    void buildUi();
    Circle readCircle(int index) const;
    void fillTable(const std::vector<TangentCircle>& sols, bool isTTR);

    QDoubleSpinBox* inputs[3][3];
    QDoubleSpinBox* extX;
    QDoubleSpinBox* extY;
    QDoubleSpinBox* extR;

    QTableWidget* table = nullptr;
    CanvasWidget* canvas = nullptr;
    std::vector<TangentCircle> lastSolutions;

    void updateInputsUi();
    QGroupBox* groupBoxes[3];
    QComboBox* modeCombo = nullptr;
    QPushButton* verifyBtn = nullptr;
};