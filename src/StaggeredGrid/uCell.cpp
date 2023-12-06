#include "uCell.h"

uCell::uCell() 
    : StaggeredCell (Point(0,0,0), 0, 0, 0, 0) {}
    uCell::uCell(Point cellCenter, size_t i, size_t j)
    : StaggeredCell (cellCenter, 0, 0, i, j) {}
uCell::uCell(Point cellCenter, double dx, double dy, size_t i, size_t j)
    : StaggeredCell(cellCenter, dx, dy, i, j) {
    }
uCell::~uCell(){}
