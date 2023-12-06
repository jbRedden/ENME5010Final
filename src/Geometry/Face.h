#pragma once
#include <vector>
#include <assert.h>
#include <stdexcept>
#include "Point.h"
#include "../Array/Vector3D.h"

class Face {
    public:
        Face();
        Face(Point *p1, Point *p2, Point *p3, Point *p4, int index);
        Face(std::vector<Point*> &points, int index);
        ~Face();

        const Point GetCenter(void) const;
        std::vector<Point*> GetPoints(void) const;
        void PrintPoints(void) const;
        int GetIndex(void) const;
        Vector3D GetVelocity(void) const;
        void SetVelocity(Vector3D u);
        void AddCellIndex(int index);
        const std::vector<int>& GetCellIndices(void) const;
        bool IsCoplanar(Face &face);
        Vector3D CalculateNormal(void);
        void SetArea(double area);
        double GetArea(void);

        friend class Point;
        friend class Cell;

    private:
        const std::vector<Point*> points;
        const Point fCtr;
        int index;
        double area;
        std::vector<int> cells;
        Vector3D velocity;

        void CheckPointPointers();
};
