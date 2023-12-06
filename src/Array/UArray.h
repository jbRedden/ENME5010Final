#pragma once

#include "../Geometry/Geometry.h"

class UArray {
    public:
        UArray(Mesh *mesh);
        ~UArray();

        double** GetU(void);
        size_t GetNX(void);
        size_t GetNY(void);

    private:
        Mesh *mesh;
        size_t nx;
        size_t ny;
        double **u;
        Face ***faces;
        std::vector<Face> psuedoFaces;

        void AllocateMemory(void);
        void PopulateFaceArray(void);

};
