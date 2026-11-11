#include "mainwindow.h"
#include <QtWidgets>
#include <QtSvg/QSvgGenerator>


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
    if (std::abs(M) < 1e-12) return {};

    double x0 = (D2*B3 - D3*B2) / M, xR = (C3*B2 - C2*B3) / M;
    double y0 = (A2*D3 - A3*D2) / M, yR = (A2*C3 - A3*C2) / M;

    double ax = x0 - c1.x, ay = y0 - c1.y;
    double alpha = xR*xR + yR*yR - 1.0;
    double beta  = 2.0*(ax*xR + ay*yR - s1*c1.r);
    double gamma = ax*ax + ay*ay - c1.r*c1.r;

    std::vector<TangentCircle> out;
    auto add = [&](double R)
    {
        if (R < -1e-9) return;
        if (R < 0) R = 0;
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

static std::vector<TangentCircle> solve_all(const Circle& c1, const Circle& c2, const Circle& c3)
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
}

void CanvasWidget::setScene(const Circle& c1, const Circle& c2, const Circle& c3, const std::vector<TangentCircle>& solutions)
{
    m_c1 = c1; m_c2 = c2; m_c3 = c3; m_sols = solutions;
    update();
}

QRectF CanvasWidget::worldBounds() const
{
    auto expand = [](QRectF& rc, const Circle& c)
    {
        rc |= QRectF(c.x - c.r, c.y - c.r, 2*c.r, 2*c.r);
    };
    QRectF rc(m_c1.x, m_c1.y, 0, 0);
    expand(rc, m_c1); expand(rc, m_c2); expand(rc, m_c3);
    for (const auto& s : m_sols) expand(rc, s.c);
    double pad = std::max(rc.width(), rc.height()) * 0.15 + 1.0;
    return rc.adjusted(-pad, -pad, pad, pad);
}

void CanvasWidget::renderScene(QPainter* p) const
{
    p->setRenderHint(QPainter::Antialiasing);
    QRectF world = worldBounds();
    QRectF view = QRectF(rect()).adjusted(10, 10, -10, -10);
    double scale = std::min(view.width() / world.width(), view.height() / world.height());
    double ox = view.center().x() - world.center().x() * scale;
    double oy = view.center().y() - world.center().y() * scale;

    auto toScreen = [&](double x, double y)
    {
        return QPointF(ox + x * scale, oy + y * scale);
    };
    auto drawCircle = [&](const Circle& c, QColor color, double width, bool dashed)
    {
        QPen pen(color, width);
        if (dashed) pen.setStyle(Qt::DashLine);
        p->setPen(pen);
        p->setBrush(Qt::NoBrush);
        p->drawEllipse(toScreen(c.x, c.y), c.r * scale, c.r * scale);
    };

    drawCircle(m_c1, Qt::black, 2, false);
    drawCircle(m_c2, Qt::black, 2, false);
    drawCircle(m_c3, Qt::black, 2, false);

    for (const auto& sol : m_sols)
    {
        bool ok = verify_solution(m_c1, m_c2, m_c3, sol);
        drawCircle(sol.c, ok ? QColor(0,160,0) : QColor(220,0,0), ok ? 1.5 : 2, !ok);
    }
}

void CanvasWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    p.fillRect(rect(), Qt::white);
    renderScene(&p);
}


MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent)
{
    buildUi();
    setWindowTitle("三切圆算法测试器");
    resize(1100, 750); //窗口
}

void MainWindow::buildUi()
{
    auto* central = new QWidget(this);
    auto* rootLayout = new QHBoxLayout(central);

    auto* left = new QVBoxLayout();
    const char* names[3] = {"已知圆 1", "已知圆 2", "已知圆 3"};
    for (int i = 0; i < 3; ++i)
    {
        auto* box = new QGroupBox(names[i]);
        auto* form = new QHBoxLayout(box);
        const char* labels[3] = {"x", "y", "r"};
        for (int j = 0; j < 3; ++j)
        {
            m_inputs[i][j] = new QDoubleSpinBox();
            m_inputs[i][j]->setRange(-1000, 1000);
            m_inputs[i][j]->setDecimals(4);
            m_inputs[i][j]->setValue(j == 2 ? 1.0 : (i == 1 && j == 0 ? 4.0 : (i == 2 && j == 1 ? 3.0 : 0.0)));
            form->addWidget(new QLabel(labels[j]));
            form->addWidget(m_inputs[i][j]);
        }
        left->addWidget(box);
    }


    auto* extBox = new QGroupBox("验证外部解");
    auto* extLayout = new QHBoxLayout(extBox);
    extX = new QDoubleSpinBox(); extX->setRange(-1000, 1000); extX->setDecimals(4);
    extY = new QDoubleSpinBox(); extY->setRange(-1000, 1000); extY->setDecimals(4);
    extR = new QDoubleSpinBox(); extR->setRange(0, 1000); extR->setDecimals(4); extR->setValue(1.0);

    extLayout->addWidget(new QLabel("X:")); extLayout->addWidget(extX);
    extLayout->addWidget(new QLabel("Y:")); extLayout->addWidget(extY);
    extLayout->addWidget(new QLabel("R:")); extLayout->addWidget(extR);
    auto* verifyBtn = new QPushButton("验证");
    extLayout->addWidget(verifyBtn);
    left->addWidget(extBox);

    auto* btnRow = new QHBoxLayout();
    auto* solveBtn = new QPushButton("求解全部 8 种组合");
    auto* svgBtn = new QPushButton("导出 SVG");
    btnRow->addWidget(solveBtn);
    btnRow->addWidget(svgBtn);
    left->addLayout(btnRow);

    table = new QTableWidget(0, 6);
    table->setHorizontalHeaderLabels({"符号 (s1,s2,s3)", "x", "y", "R", "验证", "说明"});
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
    connect(svgBtn, &QPushButton::clicked, this, &MainWindow::onExportSvg);
    connect(verifyBtn, &QPushButton::clicked, this, &MainWindow::onVerifyExternal);
}

Circle MainWindow::readCircle(int i) const
{
    return { m_inputs[i][0]->value(), m_inputs[i][1]->value(), m_inputs[i][2]->value() };
}

void MainWindow::onSolve()
{
    Circle c1 = readCircle(0), c2 = readCircle(1), c3 = readCircle(2);
    lastSolutions = solve_all(c1, c2, c3);


    if (lastSolutions.empty())
    {
        table->setRowCount(1);
        table->setItem(0, 0, new QTableWidgetItem("无解"));
        table->setItem(0, 1, new QTableWidgetItem("-"));
        table->setItem(0, 2, new QTableWidgetItem("-"));
        table->setItem(0, 3, new QTableWidgetItem("-"));
        table->setItem(0, 4, new QTableWidgetItem("失败"));
        table->setItem(0, 5, new QTableWidgetItem("几何冲突，找不到合法的相切圆"));
        table->item(0, 4)->setForeground(QColor(200,0,0));
    }
    else
    {
        fillTable(lastSolutions);
    }

    canvas->setScene(c1, c2, c3, lastSolutions);
}

void MainWindow::fillTable(const std::vector<TangentCircle>& sols)
{
    table->setRowCount(0);
    for (const auto& sol : sols)
    {
        bool ok = verify_solution(readCircle(0), readCircle(1), readCircle(2), sol);
        int row = table->rowCount();
        table->insertRow(row);
        auto set = [&](int col, const QString& text)
        {
            table->setItem(row, col, new QTableWidgetItem(text));
        };
        set(0, QString("(%1,%2,%3)").arg(sol.s1).arg(sol.s2).arg(sol.s3));
        set(1, QString::number(sol.c.x, 'f', 4));
        set(2, QString::number(sol.c.y, 'f', 4));
        set(3, QString::number(sol.c.r, 'f', 4));
        set(4, ok ? "通过" : "失败");
        set(5, ok ? "" : "相切关系不成立");
        table->item(row, 4)->setForeground(ok ? QColor(0,140,0) : QColor(200,0,0));
    }
}

void MainWindow::onVerifyExternal()
{
    Circle c1 = readCircle(0), c2 = readCircle(1), c3 = readCircle(2);

    Circle extSol;
    extSol.x = extX->value();
    extSol.y = extY->value();
    extSol.r = extR->value();


    auto guessSign = [](const Circle& sol, const Circle& known)
    {
        double d = std::hypot(sol.x - known.x, sol.y - known.y);
        double d_out = std::abs(d - (sol.r + known.r));   // 外切误差
        double d_in  = std::abs(d - std::abs(sol.r - known.r)); // 内切误差
        return (d_out < d_in) ? 1 : -1;
    };

    int s1 = guessSign(extSol, c1);
    int s2 = guessSign(extSol, c2);
    int s3 = guessSign(extSol, c3);

    TangentCircle sol = {extSol, s1, s2, s3};


    bool ok = verify_solution(c1, c2, c3, sol, 1e-3);


    table->setRowCount(1);
    auto set = [&](int col, const QString& text)
    {
        table->setItem(0, col, new QTableWidgetItem(text));
    };
    set(0, QString("(%1,%2,%3)").arg(s1).arg(s2).arg(s3));
    set(1, QString::number(extSol.x, 'f', 4));
    set(2, QString::number(extSol.y, 'f', 4));
    set(3, QString::number(extSol.r, 'f', 4));
    set(4, ok ? "通过" : "失败");
    set(5, ok ? "外部解正确" : "外部解存在误差或完全不相切");
    table->item(0, 4)->setForeground(ok ? QColor(0,140,0) : QColor(200,0,0));

    canvas->setScene(c1, c2, c3, {sol});
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

    QSize oldSize = canvas->size();
    canvas->resize(800, 800);
    canvas->renderScene(&p);
    canvas->resize(oldSize);
    p.end();

    statusBar()->showMessage("已导出: " + path, 3000);
}