#include "vCell.h"

vCell::vCell() 
    : StaggeredCell (Point(0,0,0), 0, 0, 0, 0) {
    }

vCell::vCell(Point cellCenter, size_t i, size_t j)
    : StaggeredCell (cellCenter, 0, 0, i, j) {
    }

vCell::vCell(Point cellCenter, double dx, double dy, size_t i, size_t j)
    : StaggeredCell(cellCenter, dx, dy, i, j) {
    }

vCell::~vCell(){}
