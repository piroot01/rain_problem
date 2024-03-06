#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <limits>
#include <iomanip>


class Vector
{
public:
    Vector() = default;

    Vector(const double x, const double y) :
        x(x), y(y)
    {
    }

    double getSize() const
    {
        return std::sqrt(x * x + y * y);
    }

    double dot(const Vector& vector) const
    {
        return this->x * vector.x + this->y * vector.y;
    }

    Vector operator+(const Vector vector) const
    {
        return {this->x + vector.x, this->y + vector.y};
    }

    Vector operator-(const Vector vector) const
    {
        return {this->x - vector.x, this->y - vector.y};
    }

    Vector operator*(const double scalar) const
    {
        return {this->x * scalar, this->y * scalar};
    }

    Vector operator/(const double scalar) const
    {
        return {this->x / scalar, this->y / scalar};
    }

public:
    double x = {0};
    double y = {0};

};


enum class Orientation
{
    CW,
    CCW,
};


class Thing
{
public:
    void loadVector(const Vector& vector)
    {
        m_vectors.push_back(vector);
    }

    using orientedProjection_t = std::pair<Vector, Orientation>;

    void calculateProjections()
    {
        m_projections.resize(m_vectors.size());
        std::size_t index = 0;
        for (const auto& vector : m_vectors)
        {
            Vector projection = project(vector, m_normal);
            Orientation orientation = getOrientation(projection, m_normal);
            m_projections[index] = {projection, orientation};
            ++index;
        }
    }

    void printProjections() const
    {
        for (const auto& vector : m_projections)
        {
            std::cout << "[ " << vector.first.x << " : " << vector.first.y << " ] - ";
            if (vector.second == Orientation::CCW)
                std::cout << "CCW";
            else
                std::cout << "CW";
            std::cout << '\n';
        }
    }

    double getProjectionSpan() const
    {
        std::pair<uint32_t, double> min = {0, std::numeric_limits<double>::max()};
        std::pair<uint32_t, double> max = {0, std::numeric_limits<double>::min()};

        std::size_t index = 0;
        for (const auto& projection : m_projections)
        {
            double orientedSize = (projection.second == Orientation::CCW) ? 
                projection.first.getSize() * (-1) : projection.first.getSize() * 1;
            if (orientedSize < min.second)
                min = {index, orientedSize};
            if (orientedSize > max.second)
                max = {index, orientedSize};
            index++;
        }
        //std::cout << "[ " << m_projections[min.first].first.x << " : " << m_projections[min.first].first.y << " ]\n";
        //std::cout << "[ " << m_projections[max.first].first.x << " : " << m_projections[max.first].first.y << " ]\n";
        return std::abs(max.second - min.second);
    }

    void updateNormal(const Vector& newNormal)
    {
        m_normal = newNormal;
    }

private:
    Vector project(const Vector& vector, const Vector& normal) const
    {
        return vector - (normal * (vector.dot(normal)) / (normal.getSize() * normal.getSize()));
    }

    Orientation getOrientation(const Vector& vector, const Vector& normal) const
    {
        if (normal.x * vector.y - normal.y * vector.x > 0)
            return Orientation::CCW;
        else
            return Orientation::CW;
    }


private:
    Vector m_normal;
    std::vector<Vector> m_vectors;
    std::vector<orientedProjection_t> m_projections;

};


class NumericMinimizer
{
public:
    NumericMinimizer(Thing& thing)
    {
        m_thing = &thing;
    }

    void loadInitialNormal(const Vector& normal)
    {
        m_normal = normal;
    }

    void compute()
    {
        std::size_t iter = 0;
        m_thing->updateNormal(m_normal);
        double initialVx = m_normal.x;
        std::cout << std::setprecision(20);

        for (; iter < m_iterCount; ++iter)
        {
            double spanLeft = iteration({m_normal.x - m_dx, m_normal.y});
            double spanRight = iteration({m_normal.x + m_dx, m_normal.y});
            if (spanLeft - spanRight > 0)
            {
                m_normal.x += m_dx;
            }
            else
            {
                m_normal.x -= m_dx;
            }
            if (iter % m_printEvery == 0)
                std::cout << "normal_x = " << m_normal.x << ", span: " << m_thing->getProjectionSpan()<< '\n';
        }
        m_result = {initialVx - m_normal.x, m_thing->getProjectionSpan()};
        m_normal = {initialVx, m_normal.y};
    }

    void computeOptimilized()
    {
        m_thing->updateNormal(m_normal);
        double initialX = m_normal.x;
        double minSpan = iteration(m_normal);
        
        double intervalStart = -100;
        double intervalEnd = 100;

        double dx = 0.1;
        double x = intervalStart;
        double minX = 0;
        std::uint32_t iter = 0;

        std::cout << std::setprecision(15);
        // converges pretty fast so the number of iterations should be small
        while (true)
        {
            while (x < intervalEnd)
            {
                double newSpan = iteration({x, m_normal.y});
                if (newSpan < minSpan)
                {
                    minSpan = newSpan;
                    minX = x;
                }
                x += dx;
            }
            x = minX - dx;
            intervalEnd = minX + dx;
            dx = dx / 2;
            iter++;
            if (iter % m_printEvery == 0)
                std::cout << "["<< iter << "]: " << "normal_x = " << x << ", span: " << m_thing->getProjectionSpan()<< '\n';
            if (iter > m_iterCount)
                break;
        }
        m_result = {initialX - x, m_thing->getProjectionSpan()};
    }

    auto getResult()
    {
        return m_result;
    }

private:
    double iteration(const Vector& normal)
    {
        m_thing->updateNormal(normal);
        m_thing->calculateProjections();
        return m_thing->getProjectionSpan();
    }

private:
    Thing* m_thing;
    Vector m_normal;

    uint32_t m_iterCount = {40};
    uint32_t m_printEvery = 5;
    double m_dx = {0.00000005};

    std::pair<double, double> m_result;

};


int main()
{
    // the thing vertecies
    const Vector v0 = {0, 0};
    const Vector v1 = {8, 5};
    const Vector v2 = {8, 8};
    const Vector v3 = {6, 7};
    const Vector v4 = {6, 6};

    const Vector normal = {1, -3.9};

    Thing thing;
    thing.loadVector(v0);
    thing.loadVector(v1);
    thing.loadVector(v2);
    thing.loadVector(v3);
    thing.loadVector(v4);

    NumericMinimizer numericMinimizer(thing);
    numericMinimizer.loadInitialNormal(normal);
    numericMinimizer.computeOptimilized();

    std::cout << "v_x = " << numericMinimizer.getResult().first << ", span = " << numericMinimizer.getResult().second << '\n';


    return 0;
}
