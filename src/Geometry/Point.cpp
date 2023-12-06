#include "Point.h"
#include "../Array/Vector3D.h"
#include "Face.h"
#include "CellFace.h"

Point::Point() 
    : x(0), y(0), z(0) {
    }

Point::Point(const Point &p1) {
    x = p1.x;
    y = p1.y;
    z = p1.z;
    index = p1.index;
}

Point::Point(const Point &p1, double dx, double dy, double dz) {
    x = p1.x + dx;
    y = p1.y + dy;
    z = p1.z + dz;
    index = -1;
}

Point::Point(double x, double y, double z) 
    : x(x), y(y), z(z) {
    }

Point::Point(double x, double y, double z, int index) 
    : x(x), y(y), z(z), index(index) {
    }

Point::Point(const Point &p1, const Point &p2) 
{
    x = (p1.GetX() + p2.GetX()) / 2;
    y = (p1.GetY() + p2.GetY()) / 2;
    z = (p1.GetZ() + p2.GetZ()) / 2;
}

Point::Point(const std::vector<Point> &pointList) {
    size_t pointCount = pointList.size();
    double xAvg = 0;
    double yAvg = 0;
    double zAvg = 0;

    for (size_t i = 0; i < pointCount; i++) {
        const Point p = pointList[i];
        xAvg += p.GetX();
        yAvg += p.GetY();
        zAvg += p.GetZ();
    }

    x = xAvg / pointCount;
    y = yAvg / pointCount;
    z = zAvg / pointCount;
}

Point::Point(const std::vector<Point*> &pointVec) {
    std::vector<Point> newPointVec;
    for (Point *p : pointVec) {
        newPointVec.push_back(*p);
    }
    *this = Point(newPointVec);
}

Point::Point(const std::vector<const Face*> &faceVec) {
    std::vector<Point> pointVec;
    for (size_t i = 0; i < faceVec.size(); i++) {
        const Face *face = faceVec[i];
        pointVec.push_back(face->GetCenter());
    }

    *this = Point(pointVec);
}

Point::Point(const std::vector<CellFace*> &faceVec) {
    std::vector<const Face*> baseFaces;
    for (const CellFace *cFace : faceVec) {
        const Face *face = cFace->GetFace();
        baseFaces.push_back(face);
    }
    *this = Point(baseFaces);
}

Point::~Point() {}

void Point::PrintCoordinates(void) const {
    std::cout << "\t[ " << x << ", " << y << ", " << z << " ]\n";
}

double Point::GetX(void) const {
    return x;
}

double Point::GetY(void) const {
    return y;
}

double Point::GetZ(void) const {
    return z;
}

void Point::SetX(double x) {
    this->x = x;
}

void Point::SetY(double y) {
    this->y = y;
}

void Point::SetZ(double z) {
    this->z = z;
}

int Point::GetIndex(void) {
    return index;
}

void Point::SetIndex(int i) {
    index = i;
}

Vector3D Point::operator-(const Point &p) const {
    return Vector3D(p, *this);
}

bool Point::operator==(const Point &p) const {
    if (x == p.x && y == p.y && z == p.z) return true;
    else return false;
}

bool Point::operator!=(const Point &p) const {
    if (x == p.x && y == p.y && z == p.z) return false;
    else return true;
}
