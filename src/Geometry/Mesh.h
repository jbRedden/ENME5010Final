#pragma once

#include <new>
#include "Cell.h"
#include "Block.h"
#include "../Array/Vector3D.h"

class Mesh {

    public:

        Mesh();
        Mesh(size_t nx, size_t ny, Cell *cells);
        Mesh(Block *block);
        ~Mesh();

        
        Point* GetPoints(void);
        size_t GetPointCount(void);
        size_t GetCellCount(void);
        Cell* GetCellList(void);
        size_t GetNX(void);
        size_t GetNY(void);

        size_t GetRowCount(void);

        size_t GetColumnCount(void);

        Cell ***GetArray(void);

        Cell ***GetBlockArray(void);

        std::vector<Face>& GetFaceList(void);

        std::vector<Face*>& GetBoundaryFaces(void);
        std::vector<Face*>& GetInternalFaces(void);

        void PrintFaceInfo(void);

        void PrintCellFaceInfo(void);

        void CalculateFluxes(const double &rho);

        void StaggerToCell(double **u, double **v, double **p);


    private:
        size_t nx;
        size_t ny;
        Cell ***array;
        Cell *cells;
        std::vector<Face> faces;
        std::vector<Face*> boundaryFaces;
        std::vector<Face*> internalFaces;
        std::string type;
        Block *block;
        Cell ***blockCellArray;
        Cell *cellList;
        size_t cellCount;

        void AllocateMemory(void);
        void AllocateMemoryBlock(void);
        void GenerateCells(void);
        void PopulateCells(void);
        void GenerateFaces(void);
        void GenerateFacesBlock(void);
        void GenerateCellFaces(void);
        void BlockToCells(void);
        void PopulateCellsBlock(void);
        void CollectBoundaryFaces(void);
};
