#pragma once

#include "StaggeredCell.h"

class vCell : public StaggeredCell {
    public:
        vCell();
        vCell(Point cellCenter, size_t i, size_t j);
        vCell(Point cellCenter, double dx, double dy, size_t i, size_t j);
        ~vCell();
    private:
};
