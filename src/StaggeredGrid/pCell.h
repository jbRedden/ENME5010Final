#pragma once

#include "StaggeredCell.h"

class pCell : public StaggeredCell {
    public:
        pCell();
        pCell(Point cellCenter, double dx, double dy, size_t i, size_t j);
        ~pCell();
    private:
};
