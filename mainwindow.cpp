#include "mainwindow.h"
#include <QtWidgets>
#include <QtSvg/QSvgGenerator>
#include <QFile>
#include <QTextStream>
#include <QShortcut>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

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
            if (d > r1 + r2 + eps) continue;// 两轨迹圆外离
            if (d < std::abs(r1 - r2) - eps) continue;// 一圆包含另一圆
            if (d < eps) continue;// 同心无解或无穷多解跳过

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

static std::vector<TangentCircle> dedupSolutions(std::vector<TangentCircle> sols)
{
    auto key = [](const TangentCircle& s)
    {
        return std::make_tuple(
            std::llround(s.c.x * 1e6),
            std::llround(s.c.y * 1e6),
            std::llround(s.c.r * 1e6));
    };

    std::sort(sols.begin(), sols.end(),[&](const TangentCircle& a, const TangentCircle& b)
    {return key(a) < key(b);});

    sols.erase(std::unique(sols.begin(), sols.end(),[&](const TangentCircle& a, const TangentCircle& b)
    {return key(a) == key(b);}),sols.end());

    return sols;
}

static bool modeUsesLine(DrawMode mode, int slot)
{
    switch (mode)
    {
    case DrawMode::LLL: return true;
    case DrawMode::TLL: return slot > 0;
    case DrawMode::TTL: return slot == 2;
    case DrawMode::LLP: return slot < 2;
    case DrawMode::TLP: return slot == 1;
    case DrawMode::LPP: return slot == 0;
    default: return false;
    }
}

static bool modeUsesPolygon(DrawMode mode, int slot)
{
    switch (mode)
    {
    case DrawMode::TTP: return slot == 2;
    case DrawMode::LLP: return slot == 2;
    case DrawMode::TLP: return slot == 2;
    case DrawMode::CPP: return slot > 0;
    case DrawMode::LPP: return slot > 0;
    case DrawMode::PPP: return true;
    default: return false;
    }
}

static bool modeHasPolygon(DrawMode mode)
{
    for (int i = 0; i < 3; ++i)
        if (modeUsesPolygon(mode, i))
            return true;
    return false;
}

static bool modeUsesCircle(DrawMode mode, int slot)
{
    if (mode == DrawMode::TTT)
        return true;
    if (mode == DrawMode::TTR)
        return slot < 2;
    if (mode == DrawMode::TTP)
        return slot < 2;
    if (mode == DrawMode::TLL || mode == DrawMode::TLP || mode == DrawMode::CPP)
        return slot == 0;
    if (mode == DrawMode::TTL)
        return slot < 2;
    return false;
}

static DrawMode modeFromIndex(int index)
{
    switch (index)
    {
    case 1: return DrawMode::TTR;
    case 2: return DrawMode::LLL;
    case 3: return DrawMode::TLL;
    case 4: return DrawMode::TTL;
    case 5: return DrawMode::TTP;
    case 6: return DrawMode::LLP;
    case 7: return DrawMode::TLP;
    case 8: return DrawMode::CPP;
    case 9: return DrawMode::LPP;
    case 10: return DrawMode::PPP;
    default: return DrawMode::TTT;
    }
}

static int polygonCountForSlot(DrawMode mode, int slot)
{
    int number = 0;
    for (int i = 0; i <= slot; ++i)
        if (modeUsesPolygon(mode, i))
            ++number;
    return number;
}

static Line polygonEdgeLine(const Polygon& polygon, int edgeIndex)
{
    const QPointF& start = polygon.vertices[static_cast<size_t>(edgeIndex)];
    const QPointF& end = polygon.vertices[static_cast<size_t>((edgeIndex + 1) % polygon.vertices.size())];
    const double dx = end.x() - start.x();
    const double dy = end.y() - start.y();
    const double length = std::hypot(dx, dy);
    const double a = -dy / length;
    const double b = dx / length;
    return {a, b, -(a * start.x() + b * start.y())};
}

static bool polygonTangent(const Circle& circle, const Polygon& polygon, int edgeIndex, int side, double eps = 1e-6)
{
    if (edgeIndex < 0 || edgeIndex >= static_cast<int>(polygon.vertices.size()))
        return false;
    const Line contactLine = polygonEdgeLine(polygon, edgeIndex);
    const double signedDistance = contactLine.a * circle.x + contactLine.b * circle.y + contactLine.c;
    const double scale = std::max(1.0, std::max(circle.r, std::abs(signedDistance)));
    const double tolerance = eps * scale;
    if (std::abs(signedDistance - side * circle.r) > tolerance)
        return false;

    const QPointF& start = polygon.vertices[static_cast<size_t>(edgeIndex)];
    const QPointF& end = polygon.vertices[static_cast<size_t>((edgeIndex + 1) % polygon.vertices.size())];
    const double dx = end.x() - start.x();
    const double dy = end.y() - start.y();
    const double lengthSquared = dx * dx + dy * dy;
    const double contactX = circle.x - side * circle.r * contactLine.a;
    const double contactY = circle.y - side * circle.r * contactLine.b;
    const double projection = ((contactX - start.x()) * dx + (contactY - start.y()) * dy) / lengthSquared;
    if (projection < -tolerance || projection > 1.0 + tolerance)
        return false;

    if (side > 0)
    {
        for (size_t i = 0; i < polygon.vertices.size(); ++i)
        {
            const Line edge = polygonEdgeLine(polygon, static_cast<int>(i));
            const double distance = edge.a * circle.x + edge.b * circle.y + edge.c;
            if (distance < circle.r - tolerance)
                return false;
        }
    }
    return true;
}

class VertexEditor : public QScrollArea
{
public:
    explicit VertexEditor(QWidget* parent = nullptr) : QScrollArea(parent)
    {
        setWidgetResizable(true);
        setFrameShape(QFrame::StyledPanel);
        setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        setMinimumHeight(120);
        setMaximumHeight(240);

        QWidget* container = new QWidget(this);
        m_layout = new QGridLayout(container);
        m_layout->setContentsMargins(8, 8, 8, 8);
        m_layout->setHorizontalSpacing(10);
        m_layout->setVerticalSpacing(6);
        m_layout->setColumnStretch(1, 1);
        m_layout->setColumnStretch(2, 1);
        setWidget(container);
    }

    void setVertexCount(int count)
    {
        while (m_x.count() > count)
        {
            const int row = m_x.count() - 1;
            QWidget* widgets[] = {m_labels[row], m_x[row], m_y[row]};
            for (QWidget* w : widgets)
            {
                m_layout->removeWidget(w);
                delete w;
            }
            m_labels.removeLast();
            m_x.removeLast();
            m_y.removeLast();
        }
        while (m_x.count() < count)
        {
            const int row = m_x.count();
            auto* label = new QLabel(QString("顶点 %1").arg(row + 1));
            auto* x = new QDoubleSpinBox;
            x->setRange(-1000.0, 1000.0);
            x->setDecimals(4);
            auto* y = new QDoubleSpinBox;
            y->setRange(-1000.0, 1000.0);
            y->setDecimals(4);
            m_layout->addWidget(label, row, 0);
            m_layout->addWidget(x, row, 1);
            m_layout->addWidget(y, row, 2);
            m_labels.append(label);
            m_x.append(x);
            m_y.append(y);
        }
    }

    int vertexCount() const { return m_x.count(); }

    QPointF vertex(int row) const { return {m_x[row]->value(), m_y[row]->value()}; }

    void setVertex(int row, const QPointF& value)
    {
        m_x[row]->setValue(value.x());
        m_y[row]->setValue(value.y());
    }

private:
    QGridLayout* m_layout = nullptr;
    QList<QLabel*> m_labels;
    QList<QDoubleSpinBox*> m_x;
    QList<QDoubleSpinBox*> m_y;
};

static QString sideName(double a, double b)
{
    const double deg = std::atan2(b, a) * 180.0 / M_PI;  // -180 ~ 180
    if (deg >= -22.5 && deg < 22.5)    return "右侧";
    if (deg >=  22.5 && deg < 67.5)    return "右上侧";
    if (deg >=  67.5 && deg < 112.5)   return "上侧";
    if (deg >= 112.5 && deg < 157.5)   return "左上侧";
    if (deg >= 157.5 || deg < -157.5)  return "左侧";
    if (deg >= -157.5 && deg < -112.5) return "左下侧";
    if (deg >= -112.5 && deg < -67.5)  return "下侧";
    return "右下侧";  // -67.5 ~ -22.5
}

static bool validLine(const Line& line)
{
    return std::hypot(line.a, line.b) > 0.5;
}

static int lineNumberForSlot(DrawMode mode, int slot)
{
    int number = 0;
    for (int i = 0; i <= slot; ++i)
        if (modeUsesLine(mode, i))
            ++number;
    return number;
}

static bool solutionMatches(const std::array<Circle, 3>& circles, const std::array<Line, 3>& lines, const TangentCircle& sol, DrawMode mode)
{
    if (mode == DrawMode::TTR)
        return tangent(sol.c, circles[0], sol.s1) && tangent(sol.c, circles[1], sol.s2);
    if (mode == DrawMode::TTT)
        return solution(circles[0], circles[1], circles[2], sol);

    const int signs[] = {sol.s1, sol.s2, sol.s3};
    for (int i = 0; i < 3; ++i)
    {
        if (modeUsesLine(mode, i))
        {
            if (!tangent(sol.c, lines[static_cast<size_t>(i)]))
                return false;
        }
        else if (!tangent(sol.c, circles[static_cast<size_t>(i)], signs[i]))
        {
            return false;
        }
    }
    return true;
}

static int polygonEdgeForSlot(const TangentCircle& sol, int slot)
{
    return slot == 0 ? sol.e1 : slot == 1 ? sol.e2 : sol.e3;
}

static bool solutionMatches(const std::array<Circle, 3>& circles, const std::array<Line, 3>& lines,const std::array<Polygon, 3>& polygons, const TangentCircle& sol, DrawMode mode)
{
    if (mode == DrawMode::TTR || mode == DrawMode::TTT)
        return solutionMatches(circles, lines, sol, mode);

    const int signs[] = {sol.s1, sol.s2, sol.s3};
    for (int i = 0; i < 3; ++i)
    {
        if (modeUsesLine(mode, i))
        {
            if (!tangent(sol.c, lines[static_cast<size_t>(i)]))
                return false;
        }
        else if (modeUsesPolygon(mode, i))
        {
            if (!polygonTangent(sol.c, polygons[static_cast<size_t>(i)], polygonEdgeForSlot(sol, i), signs[i]))
                return false;
        }
        else if (!tangent(sol.c, circles[static_cast<size_t>(i)], signs[i]))
        {
            return false;
        }
    }
    return true;
}

static std::vector<TangentCircle> solve_linear_circle_combo(const std::array<Circle, 3>& circles,const std::array<Line, 3>& lines,DrawMode mode, const std::array<Polygon, 3>& polygons,const std::array<int, 3>& polygonEdges)
{
    struct Constraint
    {
        bool isLine;
        bool isPolygon;
        Circle circle;
        Line line;
        int sign;
    };

    std::vector<TangentCircle> results;
    for (int mask = 0; mask < 8; ++mask)
    {
        std::array<Constraint, 3> constraints;
        int circleCount = 0;
        int anchor = -1;
        for (int i = 0; i < 3; ++i)
        {
            const int sign = (mask & (1 << i)) ? 1 : -1;
            const bool polygonConstraint = modeUsesPolygon(mode, i);
            const bool linearConstraint = modeUsesLine(mode, i) || polygonConstraint;
            constraints[i] = {linearConstraint, polygonConstraint, circles[i], lines[i], sign};
            if (!linearConstraint)
            {
                if (anchor < 0)
                    anchor = i;
                ++circleCount;
            }
        }

        std::vector<std::array<double, 4>> equations;
        for (const Constraint& constraint : constraints)
        {
            if (constraint.isLine)
            {
                equations.push_back({constraint.line.a, constraint.line.b,-static_cast<double>(constraint.sign), -constraint.line.c});
            }
        }

        if (anchor >= 0)
        {
            const Circle& base = constraints[anchor].circle;
            const double baseK = base.x * base.x + base.y * base.y - base.r * base.r;
            for (int i = 0; i < 3; ++i)
            {
                if (constraints[i].isLine || i == anchor)
                    continue;
                const Circle& circle = constraints[i].circle;
                const double k = circle.x * circle.x + circle.y * circle.y - circle.r * circle.r;
                equations.push_back({-2.0 * circle.x + 2.0 * base.x,-2.0 * circle.y + 2.0 * base.y, -2.0 * constraints[i].sign * circle.r + 2.0 * constraints[anchor].sign * base.r, -(k - baseK)});
            }
        }

        auto addCandidate = [&](double x, double y, double radius)
        {
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(radius) || radius <= 1e-9)
                return;
            TangentCircle candidate{{x, y, radius}, constraints[0].sign, constraints[1].sign, constraints[2].sign};
            for (int i = 0; i < 3; ++i)
            {
                bool valid;
                if (constraints[i].isPolygon)
                    valid = polygonTangent(candidate.c, polygons[static_cast<size_t>(i)], polygonEdges[static_cast<size_t>(i)], constraints[i].sign);
                else if (constraints[i].isLine)
                    valid = tangent(candidate.c, constraints[i].line);
                else
                    valid = tangent(candidate.c, constraints[i].circle, constraints[i].sign);
                if (!valid)
                    return;
            }
            for (const TangentCircle& existing : results)
            {
                if (std::hypot(existing.c.x - x, existing.c.y - y) < 1e-6 && std::abs(existing.c.r - radius) < 1e-6)
                    return;
            }
            candidate.e1 = polygonEdges[0];
            candidate.e2 = polygonEdges[1];
            candidate.e3 = polygonEdges[2];
            results.push_back(candidate);
        };

        if (circleCount == 0)
        {
            if (equations.size() != 3)
                continue;
            double matrix[3][4];
            for (int row = 0; row < 3; ++row)
                for (int col = 0; col < 4; ++col)
                    matrix[row][col] = equations[row][col];
            bool singular = false;
            for (int col = 0; col < 3; ++col)
            {
                int pivot = col;
                for (int row = col + 1; row < 3; ++row)
                    if (std::abs(matrix[row][col]) > std::abs(matrix[pivot][col]))
                        pivot = row;
                if (std::abs(matrix[pivot][col]) < 1e-12)
                {
                    singular = true;
                    break;
                }
                for (int j = col; j < 4; ++j)
                    std::swap(matrix[col][j], matrix[pivot][j]);
                const double divisor = matrix[col][col];
                for (int j = col; j < 4; ++j)
                    matrix[col][j] /= divisor;
                for (int row = 0; row < 3; ++row)
                {
                    if (row == col)
                        continue;
                    const double factor = matrix[row][col];
                    for (int j = col; j < 4; ++j)
                        matrix[row][j] -= factor * matrix[col][j];
                }
            }
            if (!singular)
                addCandidate(matrix[0][3], matrix[1][3], matrix[2][3]);
            continue;
        }

        if (equations.size() != 2)
            continue;
        const auto& first = equations[0];
        const auto& second = equations[1];
        const std::array<double, 3> direction
        {
            first[1] * second[2] - first[2] * second[1],
            first[2] * second[0] - first[0] * second[2],
            first[0] * second[1] - first[1] * second[0]
        };
        if (std::hypot(std::hypot(direction[0], direction[1]), direction[2]) < 1e-12)
            continue;

        double largestDeterminant = 0;
        int firstColumn = 0;
        int secondColumn = 1;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = i + 1; j < 3; ++j)
            {
                const double determinant = first[i] * second[j] - first[j] * second[i];
                if (std::abs(determinant) > largestDeterminant)
                {
                    largestDeterminant = std::abs(determinant);
                    firstColumn = i;
                    secondColumn = j;
                }
            }
        }
        if (largestDeterminant < 1e-12)
            continue;
        std::array<double, 3> point{0, 0, 0};
        point[firstColumn] = (first[3] * second[secondColumn] - first[secondColumn] * second[3]) / (first[firstColumn] * second[secondColumn] - first[secondColumn] * second[firstColumn]);
        point[secondColumn] = (first[firstColumn] * second[3] - first[3] * second[firstColumn]) / (first[firstColumn] * second[secondColumn] - first[secondColumn] * second[firstColumn]);

        const Circle& base = constraints[anchor].circle;
        const double signedRadius = constraints[anchor].sign * base.r;
        const double dx = point[0] - base.x;
        const double dy = point[1] - base.y;
        const double radiusOffset = point[2] + signedRadius;
        const double qa = direction[0] * direction[0] + direction[1] * direction[1] - direction[2] * direction[2];
        const double qb = 2.0 * (dx * direction[0] + dy * direction[1] - radiusOffset * direction[2]);
        const double qc = dx * dx + dy * dy - radiusOffset * radiusOffset;
        auto addAt = [&](double t)
        {
            addCandidate(point[0] + t * direction[0],point[1] + t * direction[1],point[2] + t * direction[2]);
        };

        if (std::abs(qa) < 1e-12)
        {
            if (std::abs(qb) > 1e-12)
                addAt(-qc / qb);
            continue;
        }
        double discriminant = qb * qb - 4.0 * qa * qc;
        if (discriminant < -1e-9)
            continue;
        discriminant = std::max(0.0, discriminant);
        const double root = std::sqrt(discriminant);
        addAt((-qb + root) / (2.0 * qa));
        if (root > 1e-12)
            addAt((-qb - root) / (2.0 * qa));
    }

    return results;
}

static std::vector<TangentCircle> solve_line_circle(const std::array<Circle, 3>& circles,const std::array<Line, 3>& lines,const std::array<Polygon, 3>& polygons, DrawMode mode)
{
    std::vector<TangentCircle> results;
    std::array<Line, 3> effectiveLines = lines;
    std::array<int, 3> polygonEdges{-1, -1, -1};
    auto visit = [&](auto&& self, int slot) -> void
    {
        if (slot == 3)
        {
            auto solutions = solve_linear_circle_combo(circles, effectiveLines, mode, polygons, polygonEdges);
            results.insert(results.end(), solutions.begin(), solutions.end());
            return;
        }
        if (!modeUsesPolygon(mode, slot))
        {
            self(self, slot + 1);
            return;
        }
        const Polygon& polygon = polygons[static_cast<size_t>(slot)];
        for (size_t edge = 0; edge < polygon.vertices.size(); ++edge)
        {
            polygonEdges[static_cast<size_t>(slot)] = static_cast<int>(edge);
            effectiveLines[static_cast<size_t>(slot)] = polygonEdgeLine(polygon, static_cast<int>(edge));
            self(self, slot + 1);
        }
        polygonEdges[static_cast<size_t>(slot)] = -1;
    };
    visit(visit, 0);
    return results;
}


CanvasWidget::CanvasWidget(QWidget* parent) : QWidget(parent)
{
    setMinimumSize(420, 420);
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);
}

void CanvasWidget::setScene(const Circle& c1, const Circle& c2, const Circle& c3,const std::vector<TangentCircle>& sols, DrawMode mode, const std::array<Line, 3>& lines, const std::array<Polygon, 3>& polygons)
{
    m_c1 = c1; m_c2 = c2; m_c3 = c3;
    m_lines = lines;
    m_polygons = polygons;
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
    QRectF rc;
    bool hasBounds = false;
    const Circle circles[] = {m_c1, m_c2, m_c3};
    for (int i = 0; i < 3; ++i)
    {
        if (!modeUsesCircle(m_mode, i))
            continue;
        const Circle& circle = circles[i];
        const QRectF circleBounds(circle.x - circle.r, circle.y - circle.r, 2 * circle.r, 2 * circle.r);
        if (!hasBounds)
        {
            rc = circleBounds;
            hasBounds = true;
        }
        else
        {
            rc |= circleBounds;
        }
    }

    if (!hasBounds)
    {
        rc = QRectF(-10, -10, 20, 20);
        hasBounds = true;
    }

    for (const auto& s : m_sols)
        expand(rc, s.c);
    if (modeHasPolygon(m_mode))
    {
        for (int i = 0; i < 3; ++i)
        {
            if (!modeUsesPolygon(m_mode, i))
                continue;
            for (const QPointF& vertex : m_polygons[static_cast<size_t>(i)].vertices)
            {
                const QRectF pointBounds(vertex.x(), vertex.y(), 0, 0);
                if (!hasBounds)
                {
                    rc = pointBounds;
                    hasBounds = true;
                }
                else
                {
                    rc |= pointBounds;
                }
            }
        }
    }

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
    auto toScreen = [&](double x, double y)
    { return worldToScreen(QPointF(x, y)); };
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

    if (modeHasPolygon(m_mode))
    {
        const QColor colors[] = {QColor(90, 120, 165), QColor(160, 105, 75), QColor(90, 155, 115)};
        p->save();
        p->setClipRect(view);
        for (int i = 0; i < 3; ++i)
        {
            const Polygon& polygon = m_polygons[static_cast<size_t>(i)];
            if (!modeUsesPolygon(m_mode, i) || polygon.vertices.size() < 3)
                continue;
            QPolygonF polygonPath;
            for (const QPointF& vertex : polygon.vertices)
                polygonPath << worldToScreen(vertex);
            p->setPen(QPen(colors[i], 2));
            p->setBrush(QColor(colors[i].red(), colors[i].green(), colors[i].blue(), 28));
            p->drawPolygon(polygonPath);
        }
        p->restore();
    }

    const QPointF bottomLeft = screenToWorld(QPointF(view.left(), view.bottom()));
    const QPointF topRight = screenToWorld(QPointF(view.right(), view.top()));
    const double worldMinX = std::min(bottomLeft.x(), topRight.x());
    const double worldMaxX = std::max(bottomLeft.x(), topRight.x());
    const double worldMinY = std::min(bottomLeft.y(), topRight.y());
    const double worldMaxY = std::max(bottomLeft.y(), topRight.y());
    p->save();
    p->setClipRect(view);
    for (int i = 0; i < 3; ++i)
    {
        if (!modeUsesLine(m_mode, i))
            continue;
        const Line& line = m_lines[static_cast<size_t>(i)];
        std::vector<QPointF> intersections;
        auto addIntersection = [&](double x, double y)
        {
            if (x < worldMinX - 1e-8 || x > worldMaxX + 1e-8|| y < worldMinY - 1e-8 || y > worldMaxY + 1e-8)
                return;
            for (const QPointF& point : intersections)
                if (std::hypot(point.x() - x, point.y() - y) < 1e-8)
                    return;
            intersections.emplace_back(x, y);
        };
        if (std::abs(line.b) > 1e-12)
        {
            addIntersection(worldMinX, -(line.a * worldMinX + line.c) / line.b);
            addIntersection(worldMaxX, -(line.a * worldMaxX + line.c) / line.b);
        }
        if (std::abs(line.a) > 1e-12)
        {
            addIntersection(-(line.b * worldMinY + line.c) / line.a, worldMinY);
            addIntersection(-(line.b * worldMaxY + line.c) / line.a, worldMaxY);
        }
        if (intersections.size() >= 2)
        {
            p->setPen(QPen(QColor(35, 40, 48), 2));
            p->drawLine(worldToScreen(intersections[0]), worldToScreen(intersections[1]));
        }
    }
    p->restore();

    const Circle circles[] = {m_c1, m_c2, m_c3};
    for (int i = 0; i < 3; ++i)
        if (modeUsesCircle(m_mode, i))
            drawCircle(circles[i], Qt::black, 2, false);

    for (size_t i = 0; i < m_sols.size(); ++i)
    {
        const auto& sol = m_sols[i];
        const bool ok = solutionMatches({m_c1, m_c2, m_c3}, m_lines, m_polygons, sol, m_mode);
        const bool selected = static_cast<int>(i) == m_selectedSolution;
        const QColor color = selected ? QColor(45, 115, 230) : (ok ? QColor(0,160,0) : QColor(220,0,0));
        const double width = selected ? 3.0 : (ok ? 1.5 : 2.0);
        drawCircle(sol.c, color, width, !selected && !ok);
    }

    for (int i = 0; i < 3; ++i)
        if (modeUsesCircle(m_mode, i))
            drawCenter(circles[i], QColor(40, 40, 40));
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
    QString distanceDescription;

    auto consider = [&](const Circle& circle, const QString& name)
    {
        const double centerDistance = std::hypot(point.x() - circle.x, point.y() - circle.y);
        const double circumferenceDistance = std::abs(centerDistance - circle.r);
        if (circumferenceDistance < nearestDistance)
        {
            nearestDistance = circumferenceDistance;
            nearestName = name;
            distanceDescription = "到圆周距离";
        }
    };

    const Circle circles[] = {m_c1, m_c2, m_c3};
    for (int i = 0; i < 3; ++i)
    {
        if (modeUsesCircle(m_mode, i))
            consider(circles[i], QString("已知圆 %1").arg(i + 1));
        if (modeUsesLine(m_mode, i) && validLine(m_lines[static_cast<size_t>(i)]))
        {
            const Line& line = m_lines[static_cast<size_t>(i)];
            const double distance = std::abs(line.a * point.x() + line.b * point.y() + line.c);
            if (distance < nearestDistance)
            {
                nearestDistance = distance;
                nearestName = QString("直线 %1").arg(lineNumberForSlot(m_mode, i));
                distanceDescription = "到直线距离";
            }
        }
    }
    for (size_t i = 0; i < m_sols.size(); ++i)
        consider(m_sols[i].c, QString("相切圆 %1").arg(i + 1));

    if (modeHasPolygon(m_mode))
    {
        for (int polygonIndex = 0; polygonIndex < 3; ++polygonIndex)
        {
            if (!modeUsesPolygon(m_mode, polygonIndex))
                continue;
            const Polygon& polygon = m_polygons[static_cast<size_t>(polygonIndex)];
            for (size_t i = 0; i < polygon.vertices.size(); ++i)
            {
                const QPointF& start = polygon.vertices[i];
                const QPointF& end = polygon.vertices[(i + 1) % polygon.vertices.size()];
                const double dx = end.x() - start.x();
                const double dy = end.y() - start.y();
                const double lengthSquared = dx * dx + dy * dy;
                const double projection = std::clamp(((point.x() - start.x()) * dx + (point.y() - start.y()) * dy) / lengthSquared, 0.0, 1.0);
                const double distance = std::hypot(point.x() - (start.x() + projection * dx),
                                                   point.y() - (start.y() + projection * dy));
                if (distance < nearestDistance)
                {
                    nearestDistance = distance;
                    nearestName = QString("多边形%1边%2").arg(polygonCountForSlot(m_mode, polygonIndex)).arg(i + 1);
                    distanceDescription = "到边界距离";
                }
            }
        }
    }

    return QString("世界坐标: (%1, %2) | 最近对象: %3 | %4: %5")
        .arg(point.x(), 0, 'f', 4)
        .arg(point.y(), 0, 'f', 4)
        .arg(nearestName)
        .arg(distanceDescription)
        .arg(nearestDistance, 0, 'f', 4);
}


MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent)
{
    buildUi();
    setWindowTitle("相切圆算法测试器");
    resize(1300, 850); //窗口
    onModeChanged(0); // 初始化模式界面
}

void MainWindow::buildUi()
{
    auto* central = new QWidget(this);
    auto* rootLayout = new QHBoxLayout(central);

    auto* left = new QVBoxLayout();
    auto* modeBox = new QGroupBox("相切测试模式");
    auto* modeLayout = new QHBoxLayout(modeBox);
    modeCombo = new QComboBox();
    modeCombo->addItems({"三圆相切", "两圆定半径", "三线相切", "两线一圆相切", "两圆一线相切","两圆一多边形", "两线一多边形", "一圆一线一多边形","一圆两多边形", "一线两多边形", "三多边形"});
    modeLayout->addWidget(modeCombo);
    left->addWidget(modeBox);
    const char* names[3] = {"已知圆 1", "已知圆 2", "已知圆 3"};
    for (int i = 0; i < 3; ++i)
    {
        groupBoxes[i] = new QGroupBox(names[i]);
        auto* groupLayout = new QVBoxLayout(groupBoxes[i]);
        entityEditors[i] = new QStackedWidget();
        auto* circlePage = new QWidget();
        auto* circleForm = new QHBoxLayout(circlePage);
        const char* labels[3] = {"x", "y", "r"};
        for (int j = 0; j < 3; ++j)
        {
            inputs[i][j] = new QDoubleSpinBox();
            inputs[i][j]->setRange(-1000, 1000);
            inputs[i][j]->setDecimals(4);
            inputs[i][j]->setValue(j == 2 ? 1.0 : (i == 1 && j == 0 ? 4.0: (i == 2 && j == 1 ? 3.0 : 0.0)));
            circleForm->addWidget(new QLabel(labels[j]));
            circleForm->addWidget(inputs[i][j]);
        }
        entityEditors[i]->addWidget(circlePage);

        auto* linePage = new QWidget();
        linePage->setToolTip("输入直线上两个不同点；正侧为从第一个点看向第二个点的左侧");
        auto* lineForm = new QGridLayout(linePage);
        const char* lineLabels[4] = {"x1", "y1", "x2", "y2"};
        const double defaults[3][4] = {
            {-3.0, -3.0, -3.0, 3.0},
            {-3.0, -3.0, 3.0, -3.0},
            {-3.0, 3.0, 3.0, -3.0}
        };
        for (int j = 0; j < 4; ++j)
        {
            lineInputs[i][j] = new QDoubleSpinBox();
            lineInputs[i][j]->setRange(-1000, 1000);
            lineInputs[i][j]->setDecimals(4);
            lineInputs[i][j]->setValue(defaults[i][j]);
            lineForm->addWidget(new QLabel(lineLabels[j]), j / 2, (j % 2) * 2);
            lineForm->addWidget(lineInputs[i][j], j / 2, (j % 2) * 2 + 1);
        }
        entityEditors[i]->addWidget(linePage);
        auto* polygonPage = new QWidget();
        auto* polygonPageLayout = new QVBoxLayout(polygonPage);
        auto* polygonCountLayout = new QHBoxLayout();
        polygonCountLayout->addWidget(new QLabel("边数："));
        polygonSideCounts[i] = new QSpinBox();
        polygonSideCounts[i]->setRange(3, 20);
        polygonSideCounts[i]->setValue(4);
        polygonCountLayout->addWidget(polygonSideCounts[i]);
        polygonCountLayout->addStretch();
        polygonPageLayout->addLayout(polygonCountLayout);

        polygonEditors[i] = new VertexEditor(polygonPage);
        polygonPageLayout->addWidget(polygonEditors[i]);

        const double centerX = i == 1 ? 5.0 : 0.0;
        const double centerY = i == 2 ? 5.0 : 0.0;
        const double xValues[] = {centerX - 2, centerX + 2, centerX + 2, centerX - 2};
        const double yValues[] = {centerY - 2, centerY - 2, centerY + 2, centerY + 2};
        polygonEditors[i]->setVertexCount(4);
        for (int v = 0; v < 4; ++v)
            polygonEditors[i]->setVertex(v, {xValues[v], yValues[v]});

        connect(polygonSideCounts[i], QOverload<int>::of(&QSpinBox::valueChanged), this,
                [this, i](int count)
                {
                    // 减少边数后再增加时保留已有顶点值
                    QList<QPointF> saved;
                    for (int row = 0; row < polygonEditors[i]->vertexCount(); ++row)
                        saved.append(polygonEditors[i]->vertex(row));
                    polygonEditors[i]->setVertexCount(count);
                    for (int row = 0; row < saved.count() && row < count; ++row)
                        polygonEditors[i]->setVertex(row, saved[row]);
                });
        entityEditors[i]->addWidget(polygonPage);
        groupLayout->addWidget(entityEditors[i]);
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
    auto* solveBtn = new QPushButton("求解全部相切组合");
    auto* svgBtn   = new QPushButton("导出 SVG");
    auto* csvBtn   = new QPushButton("导出 CSV");
    auto* resetBtn = new QPushButton("重置视图");
    btnRow->addWidget(solveBtn);
    btnRow->addWidget(svgBtn);
    btnRow->addWidget(csvBtn);
    btnRow->addWidget(resetBtn);
    left->addLayout(btnRow);

    table = new QTableWidget(0, 6);
    table->setHorizontalHeaderLabels({"类型/切向", "x", "y", "R", "验证", "说明"});
    table->horizontalHeader()->setStretchLastSection(true);
    left->addWidget(table, 1);

    auto* leftWidget = new QWidget();
    leftWidget->setLayout(left);
    leftWidget->setMaximumWidth(1000); //加宽
    table->setColumnWidth(0, 160);
    table->setColumnWidth(1, 60);    //x
    table->setColumnWidth(2, 60);    //y
    table->setColumnWidth(3, 60);    //r
    table->setColumnWidth(4, 40);
    table->setColumnWidth(5, 170);

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
    const DrawMode mode = modeFromIndex(index);
    const std::array<Line, 3> lines{readLine(0), readLine(1), readLine(2)};
    std::array<Polygon, 3> polygons;
    QString error;
    if (!readPolygons(polygons, mode, error))
        statusBar()->showMessage(error, 5000);
    canvas->setScene(readCircle(0), readCircle(1), readCircle(2), {}, mode, lines, polygons);
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
    const DrawMode mode = modeFromIndex(modeCombo->currentIndex());
    int circleNumber = 0;
    int lineNumber = 0;
    for (int i = 0; i < 3; ++i)
    {
        groupBoxes[i]->setEnabled(true);
        for (int j = 0; j < 3; ++j)
            inputs[i][j]->setEnabled(true);
        for (int j = 0; j < 4; ++j)
            lineInputs[i][j]->setEnabled(true);

        if (mode == DrawMode::TTR && i == 2)
        {
            groupBoxes[i]->setTitle("目标圆半径");
            entityEditors[i]->setCurrentIndex(0);
            inputs[i][0]->setEnabled(false);
            inputs[i][1]->setEnabled(false);
        }
        else if (modeUsesLine(mode, i))
        {
            groupBoxes[i]->setTitle(QString("直线 %1").arg(++lineNumber));
            entityEditors[i]->setCurrentIndex(1);
        }
        else if (modeUsesPolygon(mode, i))
        {
            groupBoxes[i]->setTitle(QString("多边形 %1").arg(polygonCountForSlot(mode, i)));
            entityEditors[i]->setCurrentIndex(2);
        }
        else
        {
            groupBoxes[i]->setTitle(QString("圆 %1").arg(++circleNumber));
            entityEditors[i]->setCurrentIndex(0);
        }
    }

    verifyBtn->setEnabled(true);
}

Circle MainWindow::readCircle(int i) const
{
    return { inputs[i][0]->value(), inputs[i][1]->value(), inputs[i][2]->value() };
}

Line MainWindow::readLine(int i) const
{
    const double x1 = lineInputs[i][0]->value();
    const double y1 = lineInputs[i][1]->value();
    const double x2 = lineInputs[i][2]->value();
    const double y2 = lineInputs[i][3]->value();
    const double dx = x2 - x1;
    const double dy = y2 - y1;
    const double length = std::hypot(dx, dy);
    if (length < 1e-9)
        return {};
    const double a = -dy / length;
    const double b =  dx / length;
    return {a, b, -(a * x1 + b * y1)};
}

bool MainWindow::readPolygon(int slot, Polygon& polygon, QString& error) const
{
    polygon.vertices.clear();
    if (slot < 0 || slot >= 3)
    {
        error = "多边形槽位无效";
        return false;
    }

    VertexEditor* editor = polygonEditors[slot];
    if (editor->vertexCount() != polygonSideCounts[slot]->value())
    {
        error = QString("多边形 %1 的边数与顶点数不一致").arg(polygonCountForSlot(modeFromIndex(modeCombo->currentIndex()), slot));
        return false;
    }
    for (int i = 0; i < editor->vertexCount(); ++i)
        polygon.vertices.emplace_back(editor->vertex(i));

    if (polygon.vertices.size() < 3)
    {
        error = QString("多边形 %1 至少需要 3 个顶点").arg(polygonCountForSlot(modeFromIndex(modeCombo->currentIndex()), slot));
        return false;
    }
    double twiceArea = 0;
    for (size_t i = 0; i < polygon.vertices.size(); ++i)
    {
        const QPointF& current = polygon.vertices[i];
        const QPointF& next = polygon.vertices[(i + 1) % polygon.vertices.size()];
        if (QLineF(current, next).length() < 1e-9)
        {
            error = QString("多边形 %1 不能有重合的相邻顶点").arg(polygonCountForSlot(modeFromIndex(modeCombo->currentIndex()), slot));
            return false;
        }
        twiceArea += current.x() * next.y() - next.x() * current.y();
    }
    if (std::abs(twiceArea) < 1e-9)
    {
        error = QString("多边形 %1 面积必须大于 0").arg(polygonCountForSlot(modeFromIndex(modeCombo->currentIndex()), slot));
        return false;
    }

    int orientation = 0;
    for (size_t i = 0; i < polygon.vertices.size(); ++i)
    {
        const QPointF& a = polygon.vertices[i];
        const QPointF& b = polygon.vertices[(i + 1) % polygon.vertices.size()];
        const QPointF& c = polygon.vertices[(i + 2) % polygon.vertices.size()];
        const double cross = (b.x() - a.x()) * (c.y() - b.y()) - (b.y() - a.y()) * (c.x() - b.x());
        if (std::abs(cross) < 1e-9)
        {
            error = QString("多边形 %1 不能有共线的连续边").arg(polygonCountForSlot(modeFromIndex(modeCombo->currentIndex()), slot));
            return false;
        }
        const int currentOrientation = cross > 0 ? 1 : -1;
        if (orientation != 0 && currentOrientation != orientation)
        {
            error = QString("多边形 %1 顶点必须按凸多边形边界顺序排列").arg(polygonCountForSlot(modeFromIndex(modeCombo->currentIndex()), slot));
            return false;
        }
        orientation = currentOrientation;
    }
    if (twiceArea < 0)
        std::reverse(polygon.vertices.begin(), polygon.vertices.end());

    for (size_t i = 0; i < polygon.vertices.size(); ++i)
    {
        const QPointF& a = polygon.vertices[i];
        const QPointF& b = polygon.vertices[(i + 1) % polygon.vertices.size()];
        const double dx = b.x() - a.x();
        const double dy = b.y() - a.y();
        for (size_t j = 0; j < polygon.vertices.size(); ++j)
        {
            if (j == i || j == (i + 1) % polygon.vertices.size())
                continue;
            const QPointF& point = polygon.vertices[j];
            const double cross = dx * (point.y() - a.y()) - dy * (point.x() - a.x());
            if (cross <= 1e-9)
            {
                error = QString("多边形 %1 顶点必须按凸多边形边界顺序排列，且不能有重复/共线顶点")
                    .arg(polygonCountForSlot(modeFromIndex(modeCombo->currentIndex()), slot));
                return false;
            }
        }
    }
    return true;
}

bool MainWindow::readPolygons(std::array<Polygon, 3>& polygons, DrawMode mode, QString& error) const
{
    for (int slot = 0; slot < 3; ++slot)
    {
        if (modeUsesPolygon(mode, slot) && !readPolygon(slot, polygons[static_cast<size_t>(slot)], error))
            return false;
    }
    return true;
}

void MainWindow::onSolve()
{
    const std::array<Circle, 3> circles{readCircle(0), readCircle(1), readCircle(2)};
    const std::array<Line, 3> lines{readLine(0), readLine(1), readLine(2)};
    const DrawMode mode = modeFromIndex(modeCombo->currentIndex());
    const Circle& c1 = circles[0];
    const Circle& c2 = circles[1];
    const Circle& c3 = circles[2];
    std::array<Polygon, 3> polygons;
    QString polygonError;

    table->setRowCount(0);
    lastSolutions.clear();
    if (!readPolygons(polygons, mode, polygonError))
    {
        table->setRowCount(1);
        table->setItem(0, 0, new QTableWidgetItem("多边形输入错误"));
        table->setItem(0, 5, new QTableWidgetItem(polygonError));
        table->item(0, 5)->setForeground(QColor(200, 0, 0));
        statusBar()->showMessage(polygonError, 5000);
        canvas->setScene(c1, c2, c3, {}, mode, lines, polygons);
        return;
    }
    canvas->setScene(c1, c2, c3, {}, mode, lines, polygons);

    for (int i = 0; i < 3; ++i)
    {
        if (modeUsesLine(mode, i) && !validLine(lines[static_cast<size_t>(i)]))
        {
            table->setRowCount(1);
            table->setItem(0, 0, new QTableWidgetItem("输入错误"));
            table->setItem(0, 5, new QTableWidgetItem(QString("直线 %1 的两点不能重合").arg(lineNumberForSlot(mode, i))));
            table->item(0, 5)->setForeground(QColor(200, 0, 0));
            statusBar()->showMessage(QString("直线 %1 的两点不能重合").arg(lineNumberForSlot(mode, i)), 4000);
            return;
        }
    }

    if (mode == DrawMode::TTR)
    {
        lastSolutions = solve_ttr(c1, c2, c3.r);   //第三个框的r当目标半径
    }
    else if (mode == DrawMode::TTT)
    {
        lastSolutions = solve_ttt(c1, c2, c3);
    }
    else
    {
        lastSolutions = solve_line_circle(circles, lines, polygons, mode);
    }
    lastSolutions = dedupSolutions(std::move(lastSolutions));
    if (lastSolutions.empty())
    {
        table->setRowCount(1);
        table->setItem(0, 0, new QTableWidgetItem("无解"));
        table->setItem(0, 5, new QTableWidgetItem("当前参数下找不到相切圆"));
        table->item(0, 5)->setForeground(QColor(200,0,0));
    }
    else
    {
        fillTable(lastSolutions, mode, polygons);
    }

    canvas->setScene(c1, c2, c3, lastSolutions, mode, lines, polygons);
    canvas->setSelectedSolution(lastSolutions.empty() ? -1 : 0);
    if (!lastSolutions.empty())
        table->setCurrentCell(0, 0);
}

void MainWindow::fillTable(const std::vector<TangentCircle>& sols, DrawMode mode, const std::array<Polygon, 3>& polygons)
{
    table->setRowCount(0);
    const std::array<Circle, 3> circles{readCircle(0), readCircle(1), readCircle(2)};
    const std::array<Line, 3> lines{readLine(0), readLine(1), readLine(2)};
    for (const auto& sol : sols)
    {
        const bool ok = solutionMatches(circles, lines, polygons, sol, mode);

        int row = table->rowCount();
        table->insertRow(row);
        auto set = [&](int col, const QString& text)
        {
            table->setItem(row, col, new QTableWidgetItem(text));
        };

        if (mode == DrawMode::TTR)
        {
            set(0, QString("(%1, %2, -)").arg(sol.s1).arg(sol.s2));
            set(5, ok ? QString("圆1%1，圆2%2").arg(sol.s1 == 1 ? "外切" : "内切").arg(sol.s2 == 1 ? "外切" : "内切") : "相切关系不成立");
        }
        else
        {
            const int signs[] = {sol.s1, sol.s2, sol.s3};
            QStringList col0Labels;
            QStringList col5Labels;
            for (int i = 0; i < 3; ++i)
            {
                if (modeUsesLine(mode, i))
                {
                    const Line& line = lines[static_cast<size_t>(i)];
                    const int number = lineNumberForSlot(mode, i);

                    //正侧用原法向量方向负侧把法向量取反，得到真实方向
                    const QString side = signs[i] > 0 ? sideName(line.a, line.b) : sideName(-line.a, -line.b);
                    const QString shortSide = side.left(side.size() - 1);

                    col0Labels << QString("线%1%2").arg(number).arg(shortSide);
                    col5Labels << QString("直线%1的%2").arg(number).arg(side);
                }
                else if (modeUsesPolygon(mode, i))
                {
                    const int edge = polygonEdgeForSlot(sol, i);
                    const QString label = QString("多边形%1%2(边%3)")
                        .arg(polygonCountForSlot(mode, i))
                        .arg(signs[i] > 0 ? "内切" : "外切")
                        .arg(edge + 1);
                    col0Labels << QString("P%1边%2").arg(polygonCountForSlot(mode, i)).arg(edge + 1);
                    col5Labels << label;
                }
                else
                {
                    const QString label = QString("圆%1%2").arg(i + 1).arg(signs[i] > 0 ? "外切" : "内切");
                    col0Labels << label;
                    col5Labels << label;
                }
            }
            set(0, col0Labels.join(" · "));
            set(5, ok ? QString("圆心位于：%1").arg(col5Labels.join("、")): "相切关系不成立");
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
    const std::array<Circle, 3> circles{readCircle(0), readCircle(1), readCircle(2)};
    const std::array<Line, 3> lines{readLine(0), readLine(1), readLine(2)};
    const Circle& c1 = circles[0];
    const Circle& c2 = circles[1];
    const Circle& c3 = circles[2];
    Circle extSol{extX->value(), extY->value(), extR->value()};
    const DrawMode mode = modeFromIndex(modeCombo->currentIndex());
    std::array<Polygon, 3> polygons;
    QString polygonError;

    table->setRowCount(0);
    lastSolutions.clear();
    if (!readPolygons(polygons, mode, polygonError))
    {
        table->setRowCount(1);
        table->setItem(0, 0, new QTableWidgetItem("多边形输入错误"));
        table->setItem(0, 5, new QTableWidgetItem(polygonError));
        table->item(0, 5)->setForeground(QColor(200, 0, 0));
        statusBar()->showMessage(polygonError, 5000);
        canvas->setScene(c1, c2, c3, {}, mode, lines, polygons);
        return;
    }
    canvas->setScene(c1, c2, c3, {}, mode, lines, polygons);
    for (int i = 0; i < 3; ++i)
    {
        if (modeUsesLine(mode, i) && !validLine(lines[static_cast<size_t>(i)]))
        {
            table->setRowCount(1);
            table->setItem(0, 0, new QTableWidgetItem("输入错误"));
            table->setItem(0, 5, new QTableWidgetItem(QString("直线 %1 的两点不能重合").arg(lineNumberForSlot(mode, i))));
            table->item(0, 5)->setForeground(QColor(200, 0, 0));
            statusBar()->showMessage(QString("直线 %1 的两点不能重合").arg(lineNumberForSlot(mode, i)), 4000);
            return;
        }
    }

    auto guessSign = [](const Circle& sol, const Circle& known)
    {
        double d = std::hypot(sol.x - known.x, sol.y - known.y);
        double d_out = std::abs(d - (sol.r + known.r));
        double d_in  = std::abs(d - std::abs(sol.r - known.r));
        return (d_out < d_in) ? 1 : -1;
    };

    int signs[3] = {1, 1, 1};
    int edges[3] = {-1, -1, -1};
    for (int i = 0; i < 3; ++i)
    {
        if (modeUsesLine(mode, i))
        {
            const Line& line = lines[static_cast<size_t>(i)];
            signs[i] = line.a * extSol.x + line.b * extSol.y + line.c >= 0 ? 1 : -1;
        }
        else if (modeUsesCircle(mode, i))
        {
            signs[i] = guessSign(extSol, circles[static_cast<size_t>(i)]);
        }
    }

    TangentCircle sol{extSol, signs[0], signs[1], signs[2]};
    bool ok = false;
    auto findPolygonContacts = [&](auto&& self, int slot) -> bool
    {
        if (slot == 3)
            return solutionMatches(circles, lines, polygons, sol, mode);
        if (!modeUsesPolygon(mode, slot))
            return self(self, slot + 1);
        const Polygon& polygon = polygons[static_cast<size_t>(slot)];
        for (size_t edge = 0; edge < polygon.vertices.size(); ++edge)
        {
            edges[slot] = static_cast<int>(edge);
            for (int side : {1, -1})
            {
                signs[slot] = side;
                sol.s1 = signs[0];
                sol.s2 = signs[1];
                sol.s3 = signs[2];
                sol.e1 = edges[0];
                sol.e2 = edges[1];
                sol.e3 = edges[2];
                if (self(self, slot + 1))
                    return true;
            }
        }
        edges[slot] = -1;
        return false;
    };
    ok = findPolygonContacts(findPolygonContacts, 0);
    if (!ok)
    {
        sol.s1 = signs[0];
        sol.s2 = signs[1];
        sol.s3 = signs[2];
        sol.e1 = edges[0];
        sol.e2 = edges[1];
        sol.e3 = edges[2];
    }

    table->setRowCount(1);
    auto set = [&](int col, const QString& text) { table->setItem(0, col, new QTableWidgetItem(text)); };

    if (mode == DrawMode::TTR)
        set(0, QString("(%1,%2,-)").arg(signs[0]).arg(signs[1]));
    else
    {
        QStringList col0Labels;
        QStringList col5Labels;
        for (int i = 0; i < 3; ++i)
        {
            if (modeUsesLine(mode, i))
            {
                const Line& line = lines[static_cast<size_t>(i)];
                const int number = lineNumberForSlot(mode, i);

                //正侧用原法向量方向负侧把法向量取反，得到真实方向
                const QString side = signs[i] > 0 ? sideName(line.a, line.b) : sideName(-line.a, -line.b);
                const QString shortSide = side.left(side.size() - 1);

                col0Labels << QString("线%1%2").arg(number).arg(shortSide);
                col5Labels << QString("直线%1的%2").arg(number).arg(side);
            }
            else if (modeUsesPolygon(mode, i))
            {
                const int edge = polygonEdgeForSlot(sol, i);
                const QString label = QString("多边形%1%2(边%3)")
                    .arg(polygonCountForSlot(mode, i))
                    .arg(signs[i] > 0 ? "内切" : "外切")
                    .arg(edge + 1);
                col0Labels << QString("P%1边%2").arg(polygonCountForSlot(mode, i)).arg(edge + 1);
                col5Labels << label;
            }
            else
            {
                const QString label = QString("圆%1%2").arg(i + 1).arg(signs[i] > 0 ? "外切" : "内切");
                col0Labels << label;
                col5Labels << label;
            }
        }
        set(0, col0Labels.join(" · "));
        set(5, ok ? QString("圆心位于：%1").arg(col5Labels.join("、")): "相切关系不成立");
    }
    set(1, QString::number(extSol.x, 'f', 4));
    set(2, QString::number(extSol.y, 'f', 4));
    set(3, QString::number(extSol.r, 'f', 4));
    set(4, ok ? "通过" : "失败");

    if (mode == DrawMode::TTR)
    {
        set(5, ok ? "外部解与圆1、圆2相切" : "外部解与圆1、圆2不相切");
    }
    else if (mode == DrawMode::TTT)
    {
        set(5, ok ? "外部解正确" : "外部解存在误差或完全不相切");
    }
    else
    {
        set(5, ok ? "外部解与直线/圆/多边形相切" : "外部解与直线/圆/多边形不相切");
    }

    table->item(0, 4)->setForeground(ok ? QColor(0,140,0) : QColor(200,0,0));

    lastSolutions = {sol};
    canvas->setScene(c1, c2, c3, {sol}, mode, lines, polygons);
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
        (this, "导出 CSV", "results.csv", "CSV 文件 (*.csv)");//导出表格内容为csv
    if (path.isEmpty())
        return;

    const DrawMode mode = modeFromIndex(modeCombo->currentIndex());
    std::array<Polygon, 3> polygons;
    QString polygonError;
    if (!readPolygons(polygons, mode, polygonError))
    {
        statusBar()->showMessage(polygonError, 5000);
        return;
    }

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

    out << "section,object,type,x,y,r,a,b,c\n";
    for (int i = 0; i < 3; ++i)
    {
        out << "known," << (i + 1) << ",";
        if (modeUsesLine(modeFromIndex(modeCombo->currentIndex()), i))
        {
            const Line line = readLine(i);
            out << "line,,," << ","
                << QString::number(line.a, 'g', 12) << ","
                << QString::number(line.b, 'g', 12) << ","
                << QString::number(line.c, 'g', 12) << "\n";
        }
        else if (modeUsesPolygon(mode, i))
        {
            out << "known," << (i + 1) << ",polygon,,,,,,\n";
        }
        else
        {
            const Circle circle = readCircle(i);
            const bool targetRadius = mode == DrawMode::TTR && i == 2;
            out << (targetRadius ? "target-radius," : "circle,")
                << QString::number(circle.x, 'g', 12) << ","
                << QString::number(circle.y, 'g', 12) << ","
                << QString::number(circle.r, 'g', 12) << ",,,\n";
        }
    }
    for (int polygonSlot = 0; polygonSlot < 3; ++polygonSlot)
    {
        if (!modeUsesPolygon(mode, polygonSlot))
            continue;
        const Polygon& polygon = polygons[static_cast<size_t>(polygonSlot)];
        for (size_t i = 0; i < polygon.vertices.size(); ++i)
        {
            out << "polygon," << polygonCountForSlot(mode, polygonSlot) << ",vertex " << (i + 1) << ","
                << QString::number(polygon.vertices[i].x(), 'g', 12) << ","
                << QString::number(polygon.vertices[i].y(), 'g', 12) << ",,,,\n";
        }
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
    QString path = QFileDialog::getSaveFileName(this, "导出 SVG", "tangent.svg", "SVG 文件 (*.svg)");
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