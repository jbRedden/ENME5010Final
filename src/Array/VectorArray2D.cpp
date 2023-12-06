#include "VectorArray2D.h"
#include "Vector3D.h"

VectorArray2D::VectorArray2D() {};
VectorArray2D::VectorArray2D(Mesh *mesh, Value value) 
: mesh(mesh), value(value), nx(mesh->GetRowCount()), ny(mesh->GetColumnCount()) {
    AllocateMemory();
    PopulateArray();
};

VectorArray2D::~VectorArray2D() {
    for (size_t i = 0; i < nx; i++) {
        delete[] array[i];
    }
    delete[] array;
}

Vector3D ***VectorArray2D::GetArray() {
    return array;
}

void VectorArray2D::AllocateMemory() {
    array = new Vector3D**[this->nx];
    for (size_t i = 0; i < nx; i++) {
        array[i] = new Vector3D*[ny];
    }
}

void VectorArray2D::PopulateArray() {
    switch (value) {
        case (Value::U) : PopulateUField();
    }
}

void VectorArray2D::PopulateUField() {
    Cell ***cArr = mesh->GetBlockArray();
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            Cell *cell = cArr[i][j];
            array[i][j] = &cell->velocity;
        } 
    }
}
