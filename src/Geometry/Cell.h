#pragma once
#include <vector>
#include "Face.h"
#include "CellFace.h"
#include "Point.h"
#include "../Array/Vector3D.h"

class Cell
{
    public:
        Cell();
        Cell(Face *f1, Face *f2, Face *f3, Face *f4, Face *f5, Face *f6);
        Cell(Point &p1, Point &p2, Point &p3, Point &p4, Point &p5, Point &p6, Point &p7, Point &p8);
        Cell(Point *points);
        Cell(std::vector<Point*>);
        Cell(std::vector<int> pointIndex, Point *pointList);
        Cell(std::vector<int> pointIndex, Point *pointList, int index);

        ~Cell();

        std::vector<int> GetFaceIndices(void);

        /*left, 1-bottom, 2-right, 3-top, 4-front, 5-back*/
        std::vector<CellFace>& GetCellFaces(void);

        Point GetCenter(void) const;

        double GetDX(void) const;

        double GetDY(void) const;

        double GetDZ(void) const;

        void GenerateCellFaces(void) const;

        void CheckFaceNormals(void);

        void DimensionCell(void);

        Vector3D GetVelocity(void);

        void SetVelocity(Vector3D u);

        double GetPressure(void);

        void SetPressure(double p);

        std::string GetConnectivity(void);

        CellFace GetFace(int index);

        void AddFace(int index);

        void AddCellFace(Face &f, int index);

        void PrintFaceInfo(void);

        int GetIndex(void);

        void CalculateFluxes(const double &rho);

        friend class Array2D;
        friend class VectorArray2D;
        friend class Mesh;

        CellFace eastFace;
        CellFace westFace;
        CellFace northFace;
        CellFace southFace;

    private:
        Vector3D velocity;
        int index = -1;
        Point cellCtr;
        std::vector<int> faceIndices;
        int type = 12;
        std::vector<Point*> points;
        std::vector<CellFace> cellFaces;
        double pressure = 0;

};
