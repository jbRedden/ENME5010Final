#include "Face.h"
#include "Point.h"

Face::Face() {}

Face::Face(Point *p1, Point *p2, Point *p3, Point *p4, int index) 
    : points({p1, p2, p3, p4}), fCtr(points), index(index) {
    CheckPointPointers();
    area = (*p1 - *p2).GetMagnitude() * (*p3 - *p2).GetMagnitude();
    CalculateNormal();
}

Face::Face(std::vector<Point*> &pointVec, int index) 
    : points(pointVec), fCtr(pointVec), index(index) {
    CheckPointPointers();
    assert(pointVec.size() == 4);
}

Face::~Face() {}

void Face::PrintPoints(void) const {
    for (size_t i = 0; i < points.size(); i++) {
        Point *point = points[i];
        std::cout << i << "\t";
        point->PrintCoordinates();
    }
}

void Face::AddCellIndex(int index) {
    cells.push_back(index);
}

Vector3D Face::CalculateNormal() {
    Vector3D vec1 = *points[0] - fCtr;
    Vector3D vec2 = *points[1] - fCtr;
    Vector3D norm = vec1 & vec2;
    norm.Normalize();
    return norm;
}

const Point Face::GetCenter() const {
    return fCtr;
}

int Face::GetIndex() const {
    return index;
}

const std::vector<int> &Face::GetCellIndices() const {
    return cells;
}

std::vector<Point*> Face::GetPoints() const {
    return points;
}

Vector3D Face::GetVelocity() const {
    return velocity;
}

void Face::SetVelocity(Vector3D u) {
    velocity = u;
}

double Face::GetArea(void) {
    return area;
}

void Face::SetArea(double area) {
    this->area = area;
}

void Face::CheckPointPointers() {
    for (Point *point : points) {
        if (point == nullptr) {
            throw std::runtime_error("One or more pointers are null");
        }
    }
}

bool Face::IsCoplanar(Face &face) {
    Vector3D norm = this->CalculateNormal();
    if ((norm & face.CalculateNormal()).GetMagnitude() != 0) return false;
    Vector3D vec = fCtr - face.fCtr;
    return (vec * norm == 0) ? true : false;
}

