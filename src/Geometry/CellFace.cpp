#include "CellFace.h"

CellFace::CellFace() {
}

CellFace::CellFace(const CellFace &cellFace) 
    : face(cellFace.face), normal(cellFace.normal), index(cellFace.index){
    }

CellFace::CellFace(Face *face, int index)  
    : face(face), index(index) {
        CalculateNormal();
    }

CellFace::CellFace(Face &face, Point &cellCenter, int index)  
    : face(&face), index(index) {
        SetCellCenterDistance(&cellCenter);
        CalculateNormal();
        Vector3D centerToFace = face.GetCenter() - cellCenter;
        double dot = centerToFace * normal;
        if (dot < 0) FlipNormal();
    }

CellFace::~CellFace() { }

int CellFace::GetIndex() {
    return index;
}

Face* CellFace::GetFace(void) const {
    return face;
}

const Point CellFace::GetCenter(void) {
    return face->GetCenter();
}

Vector3D CellFace::GetNormal(void) {
    return normal;
}


void CellFace::PrintPoints(void) const {
    face->PrintPoints();
}

void CellFace::PrintNormal(void) const {
    std::cout << "Cell Face Normal\t[ " 
        << normal.GetU() << ", " << normal.GetV() 
        << ", " << normal.GetW() << " ]\n";
}

void CellFace::PrintInfo(void) {
    std::cout << "Index\t" << face->GetIndex() << "\nPoints\n";
    std::cout << "Working?\n";
//    for (Point *p : face->GetPoints()) p->PrintCoordinates();
    std::cout << "Face Center";
    face->GetCenter().PrintCoordinates();
    std::cout << "Normal\n";
    normal.PrintComponents();
}

    

void CellFace::SetCellCenterDistance(Point *cellCenter) {
    Vector3D p = *cellCenter - face->GetCenter();
    d = p.GetMagnitude();
}

void CellFace::CalculateNormal(void) {
    std::vector<Point*> points = face->GetPoints();
    Vector3D faceVec1 = *points[0] - *points[1];
    Vector3D faceVec2 = *points[0] - *points[2];
    Vector3D normal = faceVec1 & faceVec2;
    normal.Normalize();

    this->normal = normal;
}

void CellFace::FlipNormal(void) {
    normal = normal * -1;
    //normal.FlipDirection();
}

void CellFace::CalculatePhi(const double &rho) {
    const Vector3D &vel = face->GetVelocity();
    phi = rho * face->GetArea() * (vel * normal);
}

std::vector<int> CellFace::GetCellIndices() {
    return face->GetCellIndices();
}
