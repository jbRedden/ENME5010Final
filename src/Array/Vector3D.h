#pragma once
#include "../Geometry/Point.h"
#include <cmath>

class Vector3D {
    public:
        Vector3D();
        Vector3D(const Vector3D &vec);
        Vector3D(double u, double v, double w);
        Vector3D(const Point &p1, const Point &p2);
        ~Vector3D();

       void Normalize(void);

       void FlipDirection(void);

       void PrintComponents(void);
       double GetU(void) const;
       double GetV(void) const;
       double GetW(void) const;
       void SetU(double);
       void SetV(double);
       void SetW(double);
       double GetMagnitude(void) const;

       double operator*(const Vector3D &vec) const;

       Vector3D operator*(const double &a) const;

       Vector3D operator&(const Vector3D &vec) const;

       Vector3D operator+(const Vector3D &vec) const;

       Vector3D operator/(const int &a) const;

       friend class Array2D;
       friend class VectorArray2D;

    private:
        double u;
        double v;
        double w;
        double mag;
};
