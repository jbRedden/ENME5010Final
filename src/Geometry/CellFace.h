#pragma once
#include <iostream>
#include "Face.h"
#include "../Array/Vector3D.h"

//Getting points from face isn't working
//Bad face pointers, indices look good though. need to clean this up or just replace
class CellFace {

    public:
        CellFace();
        CellFace(const CellFace &cellFace);
        CellFace(Face *face, int index);
        CellFace(Face &face, Point &cellCenter, int index);
        ~CellFace();

        Face* GetFace(void) const;
        const Point GetCenter(void);
        int GetIndex(void);
        Vector3D GetNormal(void);
        void PrintPoints(void) const;
        void PrintNormal(void) const;
        void SetCellCenterDistance(Point *cellCenter);
        void PrintInfo(void);
        void FlipNormal(void);
        void CalculatePhi(const double &rho);
        std::vector<int> GetCellIndices(void);

        friend class Cell;

    private:
        Face *face;
        Vector3D normal;
        int index;
        double d = -1; //distance from cell center
        double phi;

        void CalculateNormal(void);
};
