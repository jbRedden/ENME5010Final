#pragma once

#include "StaggeredCell.h"

class uCell : public StaggeredCell {
    public:
        uCell();
        uCell(Point cellCenter, size_t i, size_t j);
        uCell(Point cellCenter, double dx, double dy, size_t i, size_t j);
        ~uCell();

    private:
};
