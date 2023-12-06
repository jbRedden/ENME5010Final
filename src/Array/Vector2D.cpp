#include "Vector2D.h"

Vector2D::Vector2D(double vx, double vy) :
    vx(vx), vy(vy) {
    }

Vector2D::Vector2D(const Point &p1, const Point &p2) {
    vx = p2.GetX() - p1.GetX();
    vy = p2.GetY() - p1.GetY();
}

void Vector2D::Normalize(void) {
    double mag = std::sqrt(vx * vx + vy * vy);
    vx = vx / mag;
    vy = vy / mag;
}
