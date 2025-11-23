#include "geometry.hpp"



Vector3D Vector3D::operator+(const Vector3D &other) const
{
    return {x + other.x, y + other.y, z + other.z};
}

Vector3D Vector3D::operator-(const Vector3D &other) const
{
    return {x - other.x, y - other.y, z - other.z};
}

Vector3D Vector3D::operator-() const
{
    return {-x, -y, -z};
}

Vector3D Vector3D::operator*(double scalar) const
{
    return {x *scalar, y * scalar, z * scalar};
}

Vector3D Vector3D::operator/(double scalar) const 
{
    if (scalar == 0.0)
    {
        return *this;
    }
    return {x / scalar, y / scalar, z / scalar};
}

Vector3D &Vector3D::operator+=(const Vector3D &other)
{
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

Vector3D &Vector3D::operator-=(const Vector3D &other)
{
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}
Vector3D &Vector3D::operator*=(double scalar)
{
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}
Vector3D &Vector3D::operator/=(double scalar)
{
    if (scalar == 0.0)
    {
        return *this;
    }
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}
double Vector3D::len()
{
    return std::sqrt(x * x + y * y + z * z);
}


