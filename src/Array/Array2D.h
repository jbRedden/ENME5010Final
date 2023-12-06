#pragma once

#include "../Geometry/Geometry.h"
#include "Vector3D.h"
#include "../Boundary/Boundary.h"

class Array2D {

    public:
        enum class Value { p, u, v };

        Array2D();
        Array2D(Mesh *mesh);
        Array2D(Mesh *mesh, Value value);
        ~Array2D();

        double **GetArray(void);

        size_t GetNX(void);
        size_t GetNY(void);

        void PopulateUField(Boundary &left, Boundary &right, Boundary &top, Boundary &bottom);

    private:
        Mesh *mesh;
        Value value;
        size_t nx;
        size_t ny;
        double **array;

        void AllocateMemory(void);
        void AllocateMemory2(void);
        void PopulateArray(void);
        void PopulatePField(void);
        void PopulateVField(void);
};
