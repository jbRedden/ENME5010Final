#include "StaggeredCell.h"
#include <stdexcept>

StaggeredCell::StaggeredCell(){}
StaggeredCell::StaggeredCell(Point cellCenter, double dx, double dy, size_t i, size_t j) 
    : ctr(cellCenter), dx(dx), dy(dy), i(i), j(j) {}
StaggeredCell::~StaggeredCell(){}

const Point& StaggeredCell::GetCellCenter() const {
    return ctr;
}

const double& StaggeredCell::GetDX() const {
    return dx;
}

const double& StaggeredCell::GetDY() const {
    return dy;
}



