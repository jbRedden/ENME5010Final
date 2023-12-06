#pragma once
#include "../Geometry/Point.h"
#include <cmath>

class Vector2D {
    public:
        Vector2D();
        Vector2D(double vx, double vy);
        Vector2D(Vector2D &vec);
        Vector2D(const Point &p1, const Point &p2);
        ~Vector2D();

       void Normalize(void);

    private:
        double vx;
        double vy;
};
