#pragma once

#include "../Geometry/Cell.h"
#include "Vector3D.h"
#include "../Geometry/Mesh.h"

class VectorArray2D {

    public:
        enum class Value { U };

        VectorArray2D();
        VectorArray2D(Mesh *mesh, Value value);
        ~VectorArray2D();

        Vector3D ***GetArray(void);

    private:
        Mesh *mesh;
        Value value;
        size_t nx;
        size_t ny;
        Vector3D ***array;

        void AllocateMemory(void);
        void PopulateArray(void);
        void PopulateUField(void);
};
