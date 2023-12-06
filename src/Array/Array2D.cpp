#include "Array2D.h"

Array2D::Array2D() {};
Array2D::Array2D(Mesh *mesh, Value value) 
: mesh(mesh), value(value), nx(mesh->GetRowCount()), ny(mesh->GetColumnCount()) {
    if (value == Value::u) {
        nx += 1;
        ny += 2;
    }
    AllocateMemory();
};

Array2D::Array2D(Mesh *mesh) 
: mesh(mesh), nx(mesh->GetRowCount()), ny(mesh->GetColumnCount()) {
    AllocateMemory();
};

Array2D::~Array2D() {
    for (size_t i = 0; i < nx; i++) {
        delete[] array[i];
    }
    delete[] array;
}

double **Array2D::GetArray() {
    return array;
}

void Array2D::AllocateMemory() {
    array = new double*[this->nx];
    for (size_t i = 0; i < nx; i++) {
        array[i] = new double[ny];
    }
}

void Array2D::AllocateMemory2() {
    array = new double*[this->nx + 1];
    for (size_t i = 0; i < nx + 1; i++) {
        array[i] = new double[ny + 2];
    }
}

void Array2D::PopulateArray() {
    switch (value) {
        case (Value::p) : PopulatePField(); break;
//        case (Value::u) : PopulateUField(); break;
//        case (Value::v) : PopulateVField(); break;
    }
}

void Array2D::PopulatePField() {
    Cell ***cArr = mesh->GetBlockArray();
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            Cell *cell = cArr[i][j];
            array[i][j] = cell->pressure;
        } 
    }
}

void Array2D::PopulateUField(Boundary &left, Boundary &right, Boundary &top, Boundary &bottom) {
    Cell ***cArr = mesh->GetBlockArray();
    size_t unx = nx - 1;
    size_t uny = ny - 2;
    for (size_t j = 0; j < ny; j++) {
        array[0][j] = left.GetVelocity().GetU();
        array[nx-1][j] = right.GetVelocity().GetU();
    }
    for (size_t i = 0; i < nx; i++) {
        array[i][0] = bottom.GetVelocity().GetU();
        array[i][ny-1] = top.GetVelocity().GetU();
    }

    for (size_t j = 1; j < ny-1; j++) {
        for (size_t i = 1; i < nx-1; i++) {
            Cell *cell = cArr[i-1][j-2];
            array[i][j] = cell->eastFace.GetFace()->GetVelocity().GetU();
        } 
    }
    
}

void Array2D::PopulateVField() {
}
