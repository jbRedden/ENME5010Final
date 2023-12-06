#pragma once

#include "../Geometry/Geometry.h"
#include "pCell.h"
#include "uCell.h"
#include "vCell.h"
#include <fstream>

class StaggeredGrid {
    public:
        StaggeredGrid();
        StaggeredGrid(Mesh *mesh);
        ~StaggeredGrid();

        void PrintPGrid(void) const;
        void PrintUGrid(void) const;
        void PrintVGrid(void) const;
        void WriteCSV(void) const;
        void Testing(void);

        double GetUNX(void);
        double GetUNY(void);
        uCell** GetUCells(void);

        double GetVNX(void);
        double GetVNY(void);
        vCell** GetVCells(void);

        double GetPNX(void);
        double GetPNY(void);
        pCell** GetPCells(void);

    private:
        Mesh *mesh;
        size_t pnx;
        size_t pny;
        size_t unx;
        size_t uny;
        size_t vnx;
        size_t vny;
        pCell **pCells;
        uCell **uCells;
        vCell **vCells;

        void AllocateMemory(void);
        void PopulateGrids(void);
};
