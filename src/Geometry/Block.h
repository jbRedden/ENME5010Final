#pragma once

#include <new>
#include "Point.h"

class Block {

    public:
        Block();
        Block(double x1, double x2, double y1, 
                double y2, double z1, double z2, 
                size_t nx, size_t ny, size_t nz);
        ~Block();

        void PrintPoints(void);

    private:
        Point *points;
        double extents[6];
        size_t nx;
        size_t ny;
        size_t nz;
        size_t nPnts;

        friend class Mesh;


        void AllocateMemory(size_t PointCount);

        void MakeUniformBlock(void);


};
