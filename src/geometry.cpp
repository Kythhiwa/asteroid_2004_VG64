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
double Vector3D::len() const
{
    return std::sqrt(x * x + y * y + z * z);
}

Vector3D Vector3D::normalized() const {
        double le = len();
        if (le < 1e-18)
        { 
            return *this;
        }
        return Vector3D(x/le, y/le, z/le);
    }
    
void Vector3D::normalize() {
    double le = len();
    if (le > 1e-15) 
    {
        x /= le;
        y /= le; 
        z /= le;
    }
}
