#include "Cell.h"
#include <stdexcept>
#include <string>

Cell::Cell() {
}

Cell::Cell(std::vector<Point*> pointVec) 
    : velocity(0,0,0), cellCtr(pointVec) {
        for (Point* p : pointVec) {
            points.push_back(p);
        }
    }

Cell::Cell(std::vector<int> pointIndices, Point *pointList) 
    : velocity(0,0,0) {
        for (int i : pointIndices) {
            points.push_back(&pointList[i]);
        }
        DimensionCell();
    }

Cell::Cell(std::vector<int> pointIndices, Point *pointList, int index) 
    : velocity(0,0,0), index(index) {
        for (int i : pointIndices) {
            points.push_back(&pointList[i]);
        }
        DimensionCell();
    }

Cell::~Cell() {
}

std::vector<int> Cell::GetFaceIndices() {
    return faceIndices;
}

std::vector<CellFace>& Cell::GetCellFaces() {
    return cellFaces;
}

Point Cell::GetCenter() const {
    return cellCtr;
}

Vector3D Cell::GetVelocity() {
    return velocity;
}

void Cell::SetVelocity(Vector3D u) {
    velocity = u;
}

double Cell::GetPressure() {
    return pressure;
}

void Cell::SetPressure(double p) {
    pressure = p;
}

double Cell::GetDX() const {
    Point p1 = eastFace.GetFace()->GetCenter();
    Point p2 = westFace.GetFace()->GetCenter();
    return (p1 - p2).GetMagnitude();
}
double Cell::GetDY() const {
    Point p1 = northFace.GetFace()->GetCenter();
    Point p2 = southFace.GetFace()->GetCenter();
    return (p1 - p2).GetMagnitude();
}
double Cell::GetDZ() const {
    return 2 * cellCtr.GetZ();
}

std::string Cell::GetConnectivity() {
    std::string connectivity;
    for (Point *p : points) {
        int i = p->GetIndex();
        connectivity += std::to_string(i) + " ";
    }
    return connectivity;
}

/*
CellFace Cell::GetFace(int index) {
    for (CellFace f : faces) {
        if (f.GetFace()->index == index) return f;
    }
    throw std::runtime_error("Cell does not contain face with given index");
}
*/

void Cell::AddFace(int index) {
    faceIndices.push_back(index);
}

void Cell::AddCellFace(Face &face, int index) {
    CellFace cFace(face, cellCtr, index);
    Vector3D vec = cFace.GetCenter() - cellCtr;
    if (vec * cFace.normal < 0) cFace.FlipNormal();
    
    cellFaces.push_back(cFace);

    if (cFace.normal.GetU() == 1) eastFace = cFace;
    else if (cFace.normal.GetU() == -1) westFace = cFace;
    else if (cFace.normal.GetV() == 1) northFace = cFace;
    else if (cFace.normal.GetV() == -1) southFace = cFace;
}

void Cell::CheckFaceNormals() {
    for (CellFace &face : cellFaces) {
        Vector3D vec = face.GetCenter() - cellCtr;
       if (face.normal * vec < 0) face.FlipNormal();
    }
}

/*
void Cell::PrintFaceInfo() {
    for (CellFace face : faces) {
        face.PrintInfo();
    }
}
*/

void Cell::DimensionCell() {
    Point p(points);
    cellCtr = p;
}

int Cell::GetIndex() {
    return index;
}

void Cell::CalculateFluxes(const double &rho) {
    for (CellFace &face : cellFaces) {
        face.CalculatePhi(rho);
    }
}

