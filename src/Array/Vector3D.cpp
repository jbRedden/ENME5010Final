#include "Vector3D.h"

Vector3D::Vector3D() 
    : u(0), v(0), w(0), mag(0) {
        mag = 0;
    }

Vector3D::Vector3D(const Vector3D &vec) 
    : u(vec.u), v(vec.v), w(vec.w) {
        mag = sqrt(u*u + v*v + w*w);
    }

Vector3D::Vector3D(double u, double v, double w) :
    u(u), v(v), w(w) {
        mag = sqrt(u*u + v*v + w*w);
    }

Vector3D::Vector3D(const Point &p1, const Point &p2) {
    u = p2.GetX() - p1.GetX();
    v = p2.GetY() - p1.GetY();
    w = p2.GetZ() - p1.GetZ();
    mag = sqrt(u*u + v*v + w*w);
}

Vector3D::~Vector3D() { }

void Vector3D::PrintComponents() {
    std::cout << "\t" << u << "\t" << v << "\t" << w << std::endl;
}

double Vector3D::GetU() const {
    return u;
}

double Vector3D::GetV() const {
    return v;
}

double Vector3D::GetW() const {
    return w;
}

void Vector3D::SetU(double u) {
    this->u = u;
}

void Vector3D::SetV(double v) {
    this->v = v;
}

void Vector3D::SetW(double w) {
    this->w = w;
}

double Vector3D::GetMagnitude() const {
    return mag;
}

void Vector3D::Normalize(void) {
    u /= mag;
    v /= mag;
    w /= mag;
    mag = 1;
}

void Vector3D::FlipDirection() {
    u *= -1;
    v *= -1;
    w *= -1;
}
double Vector3D::operator*(const Vector3D &vec) const {
    return u * vec.u + v * vec.v + w * vec.w;
}

Vector3D Vector3D::operator*(const double &a) const {
    return Vector3D(a * u, a * v, a * w);
}

Vector3D Vector3D::operator&(const Vector3D &vec) const {
    double newU = v * vec.w - w * vec.v;
    double newV = w * vec.u - u * vec.w;
    double newW = u * vec.v - v * vec.u;
    return Vector3D(newU, newV, newW);
}

Vector3D Vector3D::operator+(const Vector3D &vec) const {
    return Vector3D(u + vec.u, v + vec.v, w + vec.w);
}

Vector3D Vector3D::operator/(const int &a) const {
    return Vector3D(u / a, v / a, w / a);
}
