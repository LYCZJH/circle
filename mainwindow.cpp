#include "mainwindow.h"
#include <QtWidgets>
#include <QtSvg/QSvgGenerator>
#include <QFile>
#include <QTextStream>
#include <QShortcut>
#include <algorithm>
#include <cmath>
#include <limits>


static std::vector<TangentCircle> solve_one(const Circle& c1, const Circle& c2, const Circle& c3,int s1, int s2, int s3)
{
    auto K = [](const Circle& c)
    {
        return c.x*c.x + c.y*c.y - c.r*c.r;
    };
    double K1 = K(c1), K2 = K(c2), K3 = K(c3);

    double A2 = 2*(c2.x-c1.x), B2 = 2*(c2.y-c1.y), C2 = 2.0*(s2*c2.r - s1*c1.r), D2 = K2-K1;
    double A3 = 2*(c3.x-c1.x), B3 = 2*(c3.y-c1.y), C3 = 2.0*(s3*c3.r - s1*c1.r), D3 = K3-K1;

    double M = A2*B3 - A3*B2;
    if (std::abs(M) < 1e-12)
        return {};

    double x0 = (D2*B3 - D3*B2) / M, xR = (C3*B2 - C2*B3) / M;
    double y0 = (A2*D3 - A3*D2) / M, yR = (A2*C3 - A3*C2) / M;

    double ax = x0 - c1.x, ay = y0 - c1.y;
    double alpha = xR*xR + yR*yR - 1.0;
    double beta  = 2.0*(ax*xR + ay*yR - s1*c1.r);
    double gamma = ax*ax + ay*ay - c1.r*c1.r;

    std::vector<TangentCircle> out;
    auto add = [&](double R)
    {
        if (R <= 1e-9) return;   //半径 0的点解不当作圆
        out.push_back({{x0 + xR*R, y0 + yR*R, R}, s1, s2, s3});
    };

    if (std::abs(alpha) < 1e-12)
    {
        if (std::abs(beta) > 1e-12) add(-gamma / beta);
    }
    else
    {
        double disc = beta*beta - 4*alpha*gamma;
        if (disc >= -1e-9)
        {
            if (disc < 0) disc = 0;
            double sq = std::sqrt(disc);
            add((-beta + sq) / (2*alpha));
            add((-beta - sq) / (2*alpha));
        }
    }
    return out;
}
static std::vector<TangentCircle> solve_ttr(const Circle& c1, const Circle& c2, double targetR)
{
    std::vector<TangentCircle> results;
    if (targetR <= 1e-9) return results;   //半径0不算圆

    for (int s1 : {1, -1})
    {
        for (int s2 : {1, -1})
        {
            double r1 = (s1 == 1) ? (targetR + c1.r) : std::abs(targetR - c1.r);
            double r2 = (s2 == 1) ? (targetR + c2.r) : std::abs(targetR - c2.r);

            double dx = c2.x - c1.x, dy = c2.y - c1.y;
            double d  = std::hypot(dx, dy);

            const double eps = 1e-9;
            if (d > r1 + r2 + eps) continue;             // 两轨迹圆外离
            if (d < std::abs(r1 - r2) - eps) continue;   // 一圆包含另一圆
            if (d < eps) continue;                       // 同心无解或无穷多解跳过

            double a  = (r1*r1 - r2*r2 + d*d) / (2.0*d);
            double h2 = r1*r1 - a*a;
            if (h2 < 0) h2 = 0;
            double h  = std::sqrt(h2);

            double xm = c1.x + a * dx / d;
            double ym = c1.y + a * dy / d;
            double rx = -dy * (h / d);
            double ry =  dx * (h / d);

            results.push_back({{xm + rx, ym + ry, targetR}, s1, s2, 0});
            if (h > eps)
                results.push_back({{xm - rx, ym - ry, targetR}, s1, s2, 0});
        }
    }
    return results;
}

static std::vector<TangentCircle> solve_ttt(const Circle& c1, const Circle& c2, const Circle& c3)
{
    std::vector<TangentCircle> all;
    for (int s1 : {1, -1})
        for (int s2 : {1, -1})
            for (int s3 : {1, -1})
            {
                auto part = solve_one(c1, c2, c3, s1, s2, s3);
                all.insert(all.end(), part.begin(), part.end());
            }
    return all;
}


CanvasWidget::CanvasWidget(QWidget* parent) : QWidget(parent)
{
    setMinimumSize(420, 420);
    setFocusPolicy(Qt::StrongFocus);   //允许获得焦点方便后续键盘操作
    setMouseTracking(true);
}

void CanvasWidget::setScene(const Circle& c1, const Circle& c2, const Circle& c3,const std::vector<TangentCircle>& sols, DrawMode mode)
{
    m_c1 = c1; m_c2 = c2; m_c3 = c3;
    m_sols = sols;
    m_mode = mode;
    m_selectedSolution = -1;
    resetView();
    update();
}

void CanvasWidget::setSelectedSolution(int index)
{
    m_selectedSolution = (index >= 0 && index < static_cast<int>(m_sols.size())) ? index : -1;
    update();
}

bool CanvasWidget::removeSolutionAt(int index)
{
    if (index < 0 || index >= static_cast<int>(m_sols.size()))
        return false;

    m_sols.erase(m_sols.begin() + index);
    if (m_selectedSolution == index)
        m_selectedSolution = -1;
    else if (m_selectedSolution > index)
        --m_selectedSolution;
    update();
    return true;
}

void CanvasWidget::resetView()
{
    m_zoom = 1.0;
    m_pan  = QPointF(0, 0);
    update();
}

void CanvasWidget::setView(double zoom, const QPointF& pan)
{
    m_zoom = std::clamp(zoom, 0.05, 200.0);
    m_pan  = pan;
    update();
}

QRectF CanvasWidget::worldBounds() const
{
    auto expand = [](QRectF& rc, const Circle& c)
    {
        rc |= QRectF(c.x - c.r, c.y - c.r, 2*c.r, 2*c.r);
    };
    QRectF rc(m_c1.x, m_c1.y, 0, 0);
    expand(rc, m_c1);
    expand(rc, m_c2);

    if (m_mode == DrawMode::TTT)
        expand(rc, m_c3);

    for (const auto& s : m_sols)
        expand(rc, s.c);

    double pad = std::max(rc.width(), rc.height()) * 0.15 + 1.0;
    return rc.adjusted(-pad, -pad, pad, pad);
}

QRectF CanvasWidget::viewRect() const
{
    return QRectF(rect()).adjusted(10, 10, -10, -10);
}

double CanvasWidget::viewScale() const
{
    QRectF world = worldBounds();
    QRectF view = viewRect();
    return std::min(view.width() / world.width(), view.height() / world.height()) * m_zoom;
}

QPointF CanvasWidget::worldToScreen(const QPointF& point) const
{
    QRectF world = worldBounds();
    QRectF view = viewRect();
    double scale = viewScale();
    return QPointF(view.center().x() + (point.x() - world.center().x()) * scale + m_pan.x(),
                   view.center().y() - (point.y() - world.center().y()) * scale + m_pan.y());
}

QPointF CanvasWidget::screenToWorld(const QPointF& point) const
{
    QRectF world = worldBounds();
    QRectF view = viewRect();
    double scale = viewScale();
    if (scale <= 0.0)
        return world.center();
    return QPointF(world.center().x() + (point.x() - view.center().x() - m_pan.x()) / scale,
                   world.center().y() - (point.y() - view.center().y() - m_pan.y()) / scale);
}

void CanvasWidget::renderScene(QPainter* p) const
{
    p->setRenderHint(QPainter::Antialiasing);
    QRectF view = viewRect();
    double scale = viewScale();
    auto toScreen = [&](double x, double y) { return worldToScreen(QPointF(x, y)); };
    auto drawCircle = [&](const Circle& c, QColor color, double width, bool dashed)
    {
        QPen pen(color, width);
        if (dashed) pen.setStyle(Qt::DashLine);
        p->setPen(pen);
        p->setBrush(Qt::NoBrush);
        p->drawEllipse(toScreen(c.x, c.y), c.r * scale, c.r * scale);
    };
    auto drawCenter = [&](const Circle& c, QColor color)
    {
        QPointF center = toScreen(c.x, c.y);
        p->setPen(QPen(color, 1));
        p->drawLine(center + QPointF(-4, 0), center + QPointF(4, 0));
        p->drawLine(center + QPointF(0, -4), center + QPointF(0, 4));
    };

    p->save();
    p->setClipRect(view);
    p->fillRect(view, Qt::white);

    const QPointF topLeft = screenToWorld(view.topLeft());
    const QPointF bottomRight = screenToWorld(view.bottomRight());
    const double minX = std::min(topLeft.x(), bottomRight.x());
    const double maxX = std::max(topLeft.x(), bottomRight.x());
    const double minY = std::min(topLeft.y(), bottomRight.y());
    const double maxY = std::max(topLeft.y(), bottomRight.y());
    const double rawStep = 60.0 / scale;
    const double magnitude = std::pow(10.0, std::floor(std::log10(rawStep)));
    const double fraction = rawStep / magnitude;
    const double step = (fraction <= 1.0 ? 1.0 : fraction <= 2.0 ? 2.0 : fraction <= 5.0 ? 5.0 : 10.0) * magnitude;
    const int decimals = std::clamp(static_cast<int>(std::ceil(-std::log10(step))), 0, 6);

    QFont tickFont = p->font();
    tickFont.setPointSize(8);
    p->setFont(tickFont);
    p->setPen(QPen(QColor(232, 236, 240), 1));
    for (double x = std::ceil(minX / step) * step; x <= maxX; x += step)
    {
        const double sx = toScreen(x, 0).x();
        p->drawLine(QPointF(sx, view.top()), QPointF(sx, view.bottom()));
    }
    for (double y = std::ceil(minY / step) * step; y <= maxY; y += step)
    {
        const double sy = toScreen(0, y).y();
        p->drawLine(QPointF(view.left(), sy), QPointF(view.right(), sy));
    }

    p->setPen(QPen(QColor(90, 100, 110), 1.5));
    if (minX <= 0 && maxX >= 0)
    {
        const double x0 = toScreen(0, 0).x();
        p->drawLine(QPointF(x0, view.top()), QPointF(x0, view.bottom()));
    }
    if (minY <= 0 && maxY >= 0)
    {
        const double y0 = toScreen(0, 0).y();
        p->drawLine(QPointF(view.left(), y0), QPointF(view.right(), y0));
    }

    p->setPen(QColor(80, 90, 100));
    for (double x = std::ceil(minX / step) * step; x <= maxX; x += step)
    {
        const double sx = toScreen(x, 0).x();
        p->drawText(QPointF(sx + 2, view.bottom() - 3), QString::number(x, 'f', decimals));
    }
    for (double y = std::ceil(minY / step) * step; y <= maxY; y += step)
    {
        const double sy = toScreen(0, y).y();
        p->drawText(QPointF(view.left() + 3, sy - 2), QString::number(y, 'f', decimals));
    }
    p->restore();

    if (m_mode == DrawMode::TTT)
    {
        drawCircle(m_c1, Qt::black, 2, false);
        drawCircle(m_c2, Qt::black, 2, false);
        drawCircle(m_c3, Qt::black, 2, false);
        for (size_t i = 0; i < m_sols.size(); ++i)
        {
            const auto& sol = m_sols[i];
            bool ok = solution(m_c1, m_c2, m_c3, sol);
            bool selected = static_cast<int>(i) == m_selectedSolution;
            QColor color = selected ? QColor(45, 115, 230) : (ok ? QColor(0,160,0) : QColor(220,0,0));
            drawCircle(sol.c, color, selected ? 3 : (ok ? 1.5 : 2), !selected && !ok);
        }
    }
    else
    {
        drawCircle(m_c1, Qt::black, 2, false);
        drawCircle(m_c2, Qt::black, 2, false);
        for (size_t i = 0; i < m_sols.size(); ++i)
        {
            const auto& sol = m_sols[i];
            bool ok = tangent(sol.c, m_c1, sol.s1) && tangent(sol.c, m_c2, sol.s2);
            bool selected = static_cast<int>(i) == m_selectedSolution;
            QColor color = selected ? QColor(45, 115, 230) : (ok ? QColor(0,160,0) : QColor(220,0,0));
            drawCircle(sol.c, color, selected ? 3 : 2, !selected && !ok);
        }
    }

    drawCenter(m_c1, QColor(40, 40, 40));
    drawCenter(m_c2, QColor(40, 40, 40));
    if (m_mode == DrawMode::TTT)
        drawCenter(m_c3, QColor(40, 40, 40));
    for (size_t i = 0; i < m_sols.size(); ++i)
        drawCenter(m_sols[i].c, static_cast<int>(i) == m_selectedSolution ? QColor(45, 115, 230) : QColor(0, 130, 0));
}

void CanvasWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.fillRect(rect(), Qt::white);
    renderScene(&p);
}

void CanvasWidget::wheelEvent(QWheelEvent* e)
{
    if (e->angleDelta().y() == 0)
    {
        e->ignore();
        return;
    }

    const double factor  = (e->angleDelta().y() > 0) ? 1.15 : (1.0 / 1.15);
    const double newZoom = std::clamp(m_zoom * factor, 0.05, 200.0);
    if (std::abs(newZoom - m_zoom) < 1e-9)
    {
        e->accept();
        return;
    }

    QRectF view = QRectF(rect()).adjusted(10, 10, -10, -10);

#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    const QPointF mouse = e->position();
#else
    const QPointF mouse = e->posF();
#endif

    const double ratio = newZoom / m_zoom;
    const double A = mouse.x() - view.center().x();
    const double B = mouse.y() - view.center().y();
    m_pan.setX(A * (1.0 - ratio) + m_pan.x() * ratio);
    m_pan.setY(B * (1.0 - ratio) + m_pan.y() * ratio);

    m_zoom = newZoom;
    update();
    e->accept();
}

void CanvasWidget::mousePressEvent(QMouseEvent* e)
{
    if (e->button() == Qt::LeftButton || e->button() == Qt::MiddleButton)
    {
        m_dragging = true;//按下左键或中键拖拽
        m_lastMousePos = e->pos();
        setFocus();
        setCursor(Qt::ClosedHandCursor);
        e->accept();
        return;
    }
    QWidget::mousePressEvent(e);
}

void CanvasWidget::mouseMoveEvent(QMouseEvent* e)
{
    if (m_dragging)//拖拽中平移
    {
        QPoint delta = e->pos() - m_lastMousePos;
        m_pan += QPointF(delta);
        m_lastMousePos = e->pos();
        update();
        emit hoverInfoChanged(hoverInfo(e->pos()));
        e->accept();
        return;
    }
    emit hoverInfoChanged(hoverInfo(e->pos()));
    QWidget::mouseMoveEvent(e);
}

void CanvasWidget::mouseReleaseEvent(QMouseEvent* e)
{
    if (m_dragging && (e->button() == Qt::LeftButton || e->button() == Qt::MiddleButton))
    {
        m_dragging = false;//释放鼠标结束拖拽
        unsetCursor();
        e->accept();
        return;
    }
    QWidget::mouseReleaseEvent(e);
}

void CanvasWidget::mouseDoubleClickEvent(QMouseEvent* e)
{
    if (e->button() == Qt::LeftButton)//左键重置视图
    {
        resetView();
        e->accept();
        return;
    }
    QWidget::mouseDoubleClickEvent(e);
}

void CanvasWidget::leaveEvent(QEvent* e)
{
    emit hoverInfoChanged(QString());
    QWidget::leaveEvent(e);
}

QString CanvasWidget::hoverInfo(const QPointF& screenPoint) const
{
    const QPointF point = screenToWorld(screenPoint);
    double nearestDistance = std::numeric_limits<double>::infinity();
    QString nearestName;

    auto consider = [&](const Circle& circle, const QString& name)
    {
        const double centerDistance = std::hypot(point.x() - circle.x, point.y() - circle.y);
        const double circumferenceDistance = std::abs(centerDistance - circle.r);
        if (circumferenceDistance < nearestDistance)
        {
            nearestDistance = circumferenceDistance;
            nearestName = name;
        }
    };

    consider(m_c1, "已知圆 1");
    consider(m_c2, "已知圆 2");
    if (m_mode == DrawMode::TTT)
        consider(m_c3, "已知圆 3");
    for (size_t i = 0; i < m_sols.size(); ++i)
        consider(m_sols[i].c, QString("候选解 %1").arg(i + 1));

    return QString("世界坐标: (%1, %2) | 最近圆: %3 | 到圆周距离: %4")
        .arg(point.x(), 0, 'f', 4)
        .arg(point.y(), 0, 'f', 4)
        .arg(nearestName)
        .arg(nearestDistance, 0, 'f', 4);
}


MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent)
{
    buildUi();
    setWindowTitle("相切圆算法测试器");
    resize(1100, 750); //窗口
    onModeChanged(0); // 初始化模式界面
}

void MainWindow::buildUi()
{
    auto* central = new QWidget(this);
    auto* rootLayout = new QHBoxLayout(central);

    auto* left = new QVBoxLayout();
    auto* modeBox = new QGroupBox("画圆模式");
    auto* modeLayout = new QHBoxLayout(modeBox);
    modeCombo = new QComboBox();
    modeCombo->addItems({"相切、相切、相切", "相切、相切、半径"});
    modeLayout->addWidget(modeCombo);
    left->addWidget(modeBox);
    const char* names[3] = {"已知圆 1", "已知圆 2", "已知圆 3"};
    for (int i = 0; i < 3; ++i)
    {
        groupBoxes[i] = new QGroupBox(names[i]);
        auto* form = new QHBoxLayout(groupBoxes[i]);
        const char* labels[3] = {"x", "y", "r"};
        for (int j = 0; j < 3; ++j)
        {
            inputs[i][j] = new QDoubleSpinBox();
            inputs[i][j]->setRange(-1000, 1000);
            inputs[i][j]->setDecimals(4);
            inputs[i][j]->setValue(j == 2 ? 1.0 : (i == 1 && j == 0 ? 4.0: (i == 2 && j == 1 ? 3.0 : 0.0)));
            form->addWidget(new QLabel(labels[j]));
            form->addWidget(inputs[i][j]);
        }
        left->addWidget(groupBoxes[i]);
    }

    auto* extBox = new QGroupBox("验证外部解");
    auto* extLayout = new QHBoxLayout(extBox);
    extX = new QDoubleSpinBox(); extX->setRange(-1000, 1000); extX->setDecimals(4);
    extY = new QDoubleSpinBox(); extY->setRange(-1000, 1000); extY->setDecimals(4);
    extR = new QDoubleSpinBox(); extR->setRange(0, 1000); extR->setDecimals(4); extR->setValue(1.0);

    extLayout->addWidget(new QLabel("X:")); extLayout->addWidget(extX);
    extLayout->addWidget(new QLabel("Y:")); extLayout->addWidget(extY);
    extLayout->addWidget(new QLabel("R:")); extLayout->addWidget(extR);
    verifyBtn = new QPushButton("验证");
    extLayout->addWidget(verifyBtn);
    left->addWidget(extBox);

    auto* btnRow  = new QHBoxLayout();
    auto* solveBtn = new QPushButton("求解全部 8 种组合");
    auto* svgBtn   = new QPushButton("导出 SVG");
    auto* csvBtn   = new QPushButton("导出 CSV");
    auto* resetBtn = new QPushButton("重置视图");
    btnRow->addWidget(solveBtn);
    btnRow->addWidget(svgBtn);
    btnRow->addWidget(csvBtn);
    btnRow->addWidget(resetBtn);
    left->addLayout(btnRow);

    table = new QTableWidget(0, 6);
    table->setHorizontalHeaderLabels({"s1,s2,s3", "x", "y", "R", "验证", "说明"});
    table->horizontalHeader()->setStretchLastSection(true);
    left->addWidget(table, 1);

    auto* leftWidget = new QWidget();
    leftWidget->setLayout(left);
    leftWidget->setMaximumWidth(500); //加宽

    canvas = new CanvasWidget();
    rootLayout->addWidget(leftWidget);
    rootLayout->addWidget(canvas, 1);
    setCentralWidget(central);

    connect(solveBtn, &QPushButton::clicked, this, &MainWindow::onSolve);
    connect(svgBtn,   &QPushButton::clicked, this, &MainWindow::onExportSvg);
    connect(csvBtn,   &QPushButton::clicked, this, &MainWindow::onExportCsv);
    connect(verifyBtn, &QPushButton::clicked, this, &MainWindow::onVerifyExternal);
    connect(resetBtn, &QPushButton::clicked, canvas, &CanvasWidget::resetView);
    connect(modeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::onModeChanged);
    connect(table, &QTableWidget::currentCellChanged, this, [this](int currentRow, int, int, int)
    {
        canvas->setSelectedSolution(currentRow);
    });
    connect(canvas, &CanvasWidget::hoverInfoChanged, this, [this](const QString& text)
    {
        if (text.isEmpty())
            statusBar()->clearMessage();
        else
            statusBar()->showMessage(text);
    });

    auto addDeleteSolutionShortcut = [this](QWidget* owner)
    {
        auto* shortcut = new QShortcut(QKeySequence(Qt::Key_Delete), owner);
        shortcut->setContext(Qt::WidgetShortcut);
        connect(shortcut, &QShortcut::activated, this, &MainWindow::onDeleteSelectedSolution);
    };
    addDeleteSolutionShortcut(canvas);
    addDeleteSolutionShortcut(table);
}

void MainWindow::onModeChanged(int index)
{
    updateInputsUi();
    table->setRowCount(0);
    lastSolutions.clear();
    canvas->setScene(readCircle(0), readCircle(1), readCircle(2), {},index == 0 ? DrawMode::TTT : DrawMode::TTR);
    canvas->setSelectedSolution(-1);
}

void MainWindow::onDeleteSelectedSolution()
{
    const int row = table->currentRow();
    if (row < 0 || row >= static_cast<int>(lastSolutions.size()))
        return;
    if (!canvas->removeSolutionAt(row))
        return;

    lastSolutions.erase(lastSolutions.begin() + row);
    table->removeRow(row);

    if (lastSolutions.empty())
    {
        canvas->setSelectedSolution(-1);
        return;
    }

    const int nextRow = std::min(row, static_cast<int>(lastSolutions.size()) - 1);
    table->setCurrentCell(nextRow, 0);
    canvas->setSelectedSolution(nextRow);
}

void MainWindow::updateInputsUi()
{
    int index = modeCombo->currentIndex();

    for (int i = 0; i < 3; ++i)
    {
        groupBoxes[i]->setEnabled(true);
        for (int j = 0; j < 3; ++j)
        {
            inputs[i][j]->setEnabled(true);
        }
    }
    if (index == 0)
    {
        groupBoxes[2]->setTitle("已知圆 3");
    }
    else
    {
        groupBoxes[2]->setTitle("目标半径 (R)");
        inputs[2][0]->setEnabled(false); // 禁用 x
        inputs[2][1]->setEnabled(false); // 禁用 y
    }

    verifyBtn->setEnabled(true);
}

Circle MainWindow::readCircle(int i) const
{
    return { inputs[i][0]->value(), inputs[i][1]->value(), inputs[i][2]->value() };
}

void MainWindow::onSolve()
{
    Circle c1 = readCircle(0), c2 = readCircle(1), c3 = readCircle(2);
    DrawMode mode = (modeCombo->currentIndex() == 0) ? DrawMode::TTT : DrawMode::TTR;

    table->setRowCount(0);

    if (mode == DrawMode::TTR)
    {
        lastSolutions = solve_ttr(c1, c2, c3.r);   //第三个框的r当目标半径
    }
    else
    {
        lastSolutions = solve_ttt(c1, c2, c3);
    }
    if (lastSolutions.empty())
    {
        table->setRowCount(1);
        table->setItem(0, 0, new QTableWidgetItem("无解"));
        table->setItem(0, 5, new QTableWidgetItem("当前参数下找不到相切圆"));
        table->item(0, 5)->setForeground(QColor(200,0,0));
    }
    else
    {
        fillTable(lastSolutions, mode == DrawMode::TTR);
    }

    canvas->setScene(c1, c2, c3, lastSolutions, mode);
    canvas->setSelectedSolution(lastSolutions.empty() ? -1 : 0);
    if (!lastSolutions.empty())
        table->setCurrentCell(0, 0);
}

void MainWindow::fillTable(const std::vector<TangentCircle>& sols,bool isTTR)
{
    table->setRowCount(0);
    Circle c1 = readCircle(0), c2 = readCircle(1), c3 = readCircle(2);

    for (const auto& sol : sols)
    {
        bool ok = false;
        if (isTTR)
        {
            ok = tangent(sol.c, c1, sol.s1) && tangent(sol.c, c2, sol.s2);
        }
        else
        {
            ok = solution(c1, c2, c3, sol);
        }

        int row = table->rowCount();
        table->insertRow(row);
        auto set = [&](int col, const QString& text)
        {
            table->setItem(row, col, new QTableWidgetItem(text));
        };

        if (isTTR)
        {
            set(0, QString("(%1, %2, -)").arg(sol.s1).arg(sol.s2));
            set(5, ok ? QString("圆1%1, 圆2%2").arg(sol.s1 == 1 ? "外切" : "内切") .arg(sol.s2 == 1 ? "外切" : "内切") : "相切关系不成立");
        }
        else
        {
            set(0, QString("(%1,%2,%3)").arg(sol.s1).arg(sol.s2).arg(sol.s3));
            set(5, ok ? "" : "相切关系不成立");
        }

        set(1, QString::number(sol.c.x, 'f', 4));
        set(2, QString::number(sol.c.y, 'f', 4));
        set(3, QString::number(sol.c.r, 'f', 4));
        set(4, ok ? "通过" : "失败");

        table->item(row, 4)->setForeground(ok ? QColor(0,140,0) : QColor(200,0,0));
    }
}

void MainWindow::onVerifyExternal()
{
    Circle c1 = readCircle(0), c2 = readCircle(1), c3 = readCircle(2);
    Circle extSol{extX->value(), extY->value(), extR->value()};
    DrawMode mode = (modeCombo->currentIndex() == 0) ? DrawMode::TTT : DrawMode::TTR;

    auto guessSign = [](const Circle& sol, const Circle& known)
    {
        double d = std::hypot(sol.x - known.x, sol.y - known.y);
        double d_out = std::abs(d - (sol.r + known.r));
        double d_in  = std::abs(d - std::abs(sol.r - known.r));
        return (d_out < d_in) ? 1 : -1;
    };

    int s1 = guessSign(extSol, c1);
    int s2 = guessSign(extSol, c2);
    int s3 = guessSign(extSol, c3);

    TangentCircle sol = {extSol, s1, s2, s3};

    bool ok = (mode == DrawMode::TTR) ? (tangent(sol.c, c1, sol.s1) && tangent(sol.c, c2, sol.s2)): solution(c1, c2, c3, sol);

    table->setRowCount(1);
    auto set = [&](int col, const QString& text) { table->setItem(0, col, new QTableWidgetItem(text)); };

    set(0, mode == DrawMode::TTR ? QString("(%1,%2,-)").arg(s1).arg(s2) : QString("(%1,%2,%3)").arg(s1).arg(s2).arg(s3));
    set(1, QString::number(extSol.x, 'f', 4));
    set(2, QString::number(extSol.y, 'f', 4));
    set(3, QString::number(extSol.r, 'f', 4));
    set(4, ok ? "通过" : "失败");

    if (mode == DrawMode::TTR)
    {
        set(5, ok ? "外部解与圆1、圆2相切" : "外部解与圆1、圆2不相切");
    }
    else
    {
        set(5, ok ? "外部解正确" : "外部解存在误差或完全不相切");
    }

    table->item(0, 4)->setForeground(ok ? QColor(0,140,0) : QColor(200,0,0));

    lastSolutions = {sol};
    canvas->setScene(c1, c2, c3, {sol}, mode);
    canvas->setSelectedSolution(0);
}

static QString csvEscape(const QString& s)
{
    if (s.contains(',') || s.contains('"') || s.contains('\n') || s.contains('\r'))
    {
        QString t = s;
        t.replace('"', "\"\"");
        return '"' + t + '"';
    }
    return s;
}

void MainWindow::onExportCsv()
{
    QString path = QFileDialog::getSaveFileName
        (this, "导出 CSV", "tri_tangent_results.csv", "CSV 文件 (*.csv)");//导出表格内容为 CSV
    if (path.isEmpty())
        return;

    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        statusBar()->showMessage("无法写入: " + path, 3000);
        return;
    }

    QTextStream out(&file);
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
    out.setEncoding(QStringConverter::Utf8);
#else
    out.setCodec("UTF-8");
#endif
    out << QChar(0xFEFF);

    out << "# Known circles\n";
    out << "circle,x,y,r\n";
    for (int i = 0; i < 3; ++i)
    {
        Circle c = readCircle(i);
        out << "C" << (i + 1) << ","
            << QString::number(c.x, 'g', 12) << ","
            << QString::number(c.y, 'g', 12) << ","
            << QString::number(c.r, 'g', 12) << "\n";
    }
    out << "\n";

    int cols = table->columnCount();// 表格内容
    for (int c = 0; c < cols; ++c)
    {
        if (c > 0) out << ",";
        auto* h = table->horizontalHeaderItem(c);
        out << csvEscape(h ? h->text() : QString());
    }
    out << "\n";

    for (int r = 0; r < table->rowCount(); ++r)
    {
        for (int c = 0; c < cols; ++c)
        {
            if (c > 0) out << ",";
            auto* item = table->item(r, c);
            out << csvEscape(item ? item->text() : QString());
        }
        out << "\n";
    }

    file.close();
    statusBar()->showMessage("已导出 CSV: " + path, 3000);
}

void MainWindow::onExportSvg()
{
    QString path = QFileDialog::getSaveFileName(this, "导出 SVG", "tri_tangent.svg", "SVG 文件 (*.svg)");
    if (path.isEmpty())
        return;

    QSvgGenerator gen;
    gen.setFileName(path);
    gen.setSize(QSize(800, 800));
    gen.setViewBox(QRect(0, 0, 800, 800));

    QPainter p(&gen);
    p.fillRect(QRect(0, 0, 800, 800), Qt::white);

    p.setPen(QPen(QColor(230, 230, 230), 1));
    for (int i = 0; i < 800; i += 50)
    {
        p.drawLine(i, 0, i, 800);
        p.drawLine(0, i, 800, i);
    }

    QSize   oldSize = canvas->size();
    double  oldZoom = canvas->zoom();
    QPointF oldPan  = canvas->pan();

    canvas->resetView();
    canvas->resize(800, 800);
    canvas->renderScene(&p);
    canvas->resize(oldSize);
    canvas->setView(oldZoom, oldPan);

    p.end();

    statusBar()->showMessage("已导出: " + path, 3000);
}