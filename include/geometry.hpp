#pragma once
#include <cmath>

struct Vector3D
{
    double x;
    double y;
    double z;
    
    Vector3D()
        : x(0), y(0), z(0) {}

    Vector3D(double x, 
             double y, 
             double z)
        : x(x), y(y), z(z) {}

    Vector3D operator+(const Vector3D &other) const;
    Vector3D operator-(const Vector3D &other) const;
    Vector3D operator-() const;
    Vector3D operator*(double scalar) const;
    Vector3D operator/(double scalar) const;

    Vector3D &operator+=(const Vector3D &other);
    Vector3D &operator-=(const Vector3D &other);
    Vector3D &operator*=(double scalar);
    Vector3D &operator/=(double scalar);

    
    Vector3D normalized() const;
    void normalize();


    double len() const;

};

struct StateVector
{
    Vector3D x;
    Vector3D v;
};


struct BodyVector
{
    Vector3D x;
    Vector3D v;
    Vector3D a;
};
