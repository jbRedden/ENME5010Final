#pragma once
#include <vector>
#include <iostream>

class Vector3D;
class Face;
class CellFace;

class Point
{
    public:
        Point();
        Point(const Point &p1);
        Point(const Point &p1, double dx, double dy, double dz);
        Point(double x, double y, double z);
        Point(double x, double y, double z, int index);
        Point(const Point &p1, const Point &p2);
        Point(const std::vector<Point> &pointVec);
        Point(const std::vector<Point*> &pointVec);
        Point(const std::vector<const Face*> &faceVec);
        Point(const std::vector<CellFace*> &faceVec);
        ~Point();

        Vector3D operator-(const Point& p) const;
        bool operator==(const Point& p) const;
        bool operator!=(const Point& p) const;
        double GetX(void) const;
        double GetY(void) const;
        double GetZ(void) const;
        void PrintCoordinates(void) const;
        void SetIndex(int index);
        int GetIndex(void);

        friend class Face;
        friend class Cell;

    private:
        void SetX(double x);
        void SetY(double y);
        void SetZ(double z);
        double x;
        double y;
        double z;
        int index = -1;
};
