#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>


class Vector
{
public:
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
    Thing(const Vector& normal) :
        m_normal(normal)
    {
    }

    void loadVector(const Vector& vector)
    {
        m_vectors.push_back(vector);
    }

    using orientedProjection_t = std::pair<Vector, Orientation>;

    void calculateProjections()
    {
        for (const auto& vector : m_vectors)
        {
            Vector projection = project(vector, m_normal);
            Orientation orientation = getOrientation(projection, m_normal);
            m_projections.push_back({projection, orientation});
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

int main()
{
    // the thing vertecies
    const Vector v0 = {0, 0};
    const Vector v1 = {8, 5};
    const Vector v2 = {8, 8};
    const Vector v3 = {6, 7};
    const Vector v4 = {6, 6};

    const Vector normal = {1, -3.9};

    Thing thing(normal);
    thing.loadVector(v0);
    thing.loadVector(v1);
    thing.loadVector(v2);
    thing.loadVector(v3);
    thing.loadVector(v4);

    thing.calculateProjections();
    thing.printProjections();



    return 0;
}
