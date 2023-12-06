#include "pCell.h"

pCell::pCell(){}
pCell::pCell(Point cellCenter, double dx, double dy, size_t i, size_t j)
    : StaggeredCell(cellCenter, dx, dy, i, j) {
    }
pCell::~pCell(){}
