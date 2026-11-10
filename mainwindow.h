#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#pragma once
#include <QMainWindow>
#include <QWidget>
#include <cmath>
#include <vector>

struct Circle
{
    double x = 0, y = 0, r = 1;
};
struct TangentCircle
{
    Circle c;
    int s1 = 1, s2 = 1, s3 = 1;
};

inline bool is_tangent(const Circle& sol, const Circle& known, int s, double eps = 1e-6)
{
    double d = std::hypot(sol.x - known.x, sol.y - known.y);
    double expected = (s > 0) ? (sol.r + known.r) : std::abs(sol.r - known.r);
    double scale = std::max(1.0, std::max({sol.r, known.r, d}));
    return std::abs(d - expected) < eps * scale;
}

inline bool verify_solution(const Circle& c1, const Circle& c2, const Circle& c3,const TangentCircle& sol, double eps = 1e-6)
{
    return is_tangent(sol.c, c1, sol.s1, eps) && is_tangent(sol.c, c2, sol.s2, eps) && is_tangent(sol.c, c3, sol.s3, eps);
}


class CanvasWidget : public QWidget
{
    Q_OBJECT

public:
    explicit CanvasWidget(QWidget* parent = nullptr);
    void setScene(const Circle& c1, const Circle& c2, const Circle& c3,const std::vector<TangentCircle>& solutions);
    void renderScene(QPainter* p) const;

protected:
    void paintEvent(QPaintEvent*) override;

private:
    Circle m_c1{0,0,1}, m_c2{4,0,1}, m_c3{2,3,1};
    std::vector<TangentCircle> m_sols;
    QRectF worldBounds() const;
};


class QDoubleSpinBox;
class QTableWidget;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);

private slots:
    void onSolve();
    void onExportSvg();

private:
    void buildUi();
    Circle readCircle(int index) const;
    void fillTable(const std::vector<TangentCircle>& sols);

    QDoubleSpinBox* m_inputs[3][3];
    QTableWidget* m_table = nullptr;
    CanvasWidget* m_canvas = nullptr;
    std::vector<TangentCircle> m_lastSolutions;
};
#endif // MAINWINDOW_H
