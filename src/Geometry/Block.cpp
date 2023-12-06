#include "Block.h"

Block::Block() {
}

Block::Block(double x1, double x2, double y1, double y2, 
        double z1, double z2, size_t nx, size_t ny, size_t nz) 
    : nx(nx), ny(ny), nz(nz), extents {x1, x2, y1, y2, z1, z2} {
        nPnts = (nx + 1) * (ny + 1) * (nz + 1);
        AllocateMemory(nPnts);
        MakeUniformBlock();
}

Block::~Block() {
    delete[] points;
}

void Block::PrintPoints() {
    for (size_t i = 0; i < nPnts; i++) {
        points[i].PrintCoordinates();
    }
}

void Block::AllocateMemory(size_t PointCount) {
    points = new Point[PointCount];

    if (points == nullptr) {
        throw std::bad_alloc();
    }
}

void Block::MakeUniformBlock() {
    double dx = std::abs(extents[1] - extents[0]) / nx; 
    double dy = std::abs(extents[3] - extents[2]) / ny;
    double dz = std::abs(extents[5] - extents[4]) / nz; 
    
    /*make block*/
    int index = 0;
    for (size_t i = 0; i <= nz; i++) {
        for (size_t j = 0; j <= ny; j++) {
            for (size_t k = 0; k <= nx; k++) {
                double z = i * dz;
                double y = j * dy;
                double x = k * dx;
                Point p(x, y, z, index);
                points[index++] = p;
            }
        }
    }
}
