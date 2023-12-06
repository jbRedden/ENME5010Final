#include "UArray.h"

UArray::UArray(Mesh *mesh) 
    : mesh(mesh), nx(mesh->GetNX() + 1), ny(mesh->GetNY() + 2) {
        UArray::AllocateMemory();
}

UArray::~UArray() {
    for (size_t i = 0; i < nx; i++) {
        delete[] u[i];
        delete[] faces[i];
    }
    delete[] u;
    delete[] faces;
}

double** UArray::GetU() {
    return u;
}

size_t UArray::GetNX() {
    return nx;
}

size_t UArray::GetNY() {
    return ny;
}

void UArray::AllocateMemory() {
    u = new double*[nx];
    faces = new Face**[nx];
    for (size_t i = 0; i < nx; i++) {
        u[i] = new double[ny];
        faces[i] = new Face*[ny];
    }
}
