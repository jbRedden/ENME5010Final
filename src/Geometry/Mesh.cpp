#include "Mesh.h"
#include "../Array/Vector3D.h"

Mesh::Mesh() {};

Mesh::Mesh(size_t nx, size_t ny, Cell *cells) 
    : nx(nx), ny(ny), cells(cells) {
        AllocateMemory();
        PopulateCells();
//        GenerateFaces();
    }

Mesh::Mesh(Block *block) 
    : nx(block->nx), ny(block->ny), block(block) {
        cellCount = nx * ny;
        AllocateMemoryBlock();
        BlockToCells();
        GenerateFacesBlock();
        CollectBoundaryFaces();
        GenerateCellFaces();
    }

Mesh::~Mesh() {
    /*
    std::cout << "test1\n";
    if (array != nullptr) {
        for (size_t i = 0; i < nx; i++) {
            delete[] array[i];
        }
        delete[] array;
        array = nullptr;
    }
    */
    if (blockCellArray != nullptr) {
        for (size_t i = 0; i < nx; i++) {
            delete[] blockCellArray[i];
        }
        delete[] blockCellArray;
        blockCellArray = nullptr;
    }

    if (cellList != nullptr) {
        delete[] cellList;
        cellList = nullptr;
    }
}

void Mesh::AllocateMemory() {
    array = new Cell**[this->nx];
    for (size_t i = 0; i < nx; i++) {
        array[i] = new Cell*[ny];
    }
}

void Mesh::AllocateMemoryBlock() {
    cellList = new Cell[cellCount];
    blockCellArray = new Cell**[nx];
    for (size_t i = 0; i < nx; i++) {
        blockCellArray[i] = new Cell*[ny];
    }

    if (blockCellArray == nullptr) throw std::bad_alloc();
}

Point* Mesh::GetPoints() {
    return block->points;
}

size_t Mesh::GetPointCount() {
    return block->nPnts;
}

size_t Mesh::GetCellCount() {
    return cellCount;
}

Cell* Mesh::GetCellList() {
    return cellList;
}

size_t Mesh::GetNX() {
    return nx;
}

size_t Mesh::GetNY() {
    return ny;
}

void Mesh::BlockToCells() {
    size_t pIndex = 0;
    size_t cCount = 0;
    size_t nz = block->nz;
    size_t xyCount = (nx + 1) * (ny + 1);
    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                std::vector<int> indices;
                indices.push_back(pIndex);
                indices.push_back(pIndex + 1);
                indices.push_back(pIndex + nx + 2);
                indices.push_back(pIndex + nx + 1);
                indices.push_back(pIndex + xyCount);
                indices.push_back(pIndex + xyCount + 1);
                indices.push_back(pIndex + xyCount + nx + 2);
                indices.push_back(pIndex + xyCount + nx + 1);
                Cell cell(indices, &block->points[0], cCount);
                cellList[cCount] = cell;
                blockCellArray[i][j] = &cellList[cCount++];
                pIndex++;
            }
            pIndex++;
        }
        pIndex += nx + 1;
    }
}

void Mesh::PopulateCellsBlock() {
    size_t n = 0;
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            blockCellArray[i][j] = &cellList[n++];
        }
    }
}

void Mesh::PopulateCells() {
    size_t n = 0;
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            array[i][j] = &cells[n++];
        }
    }
}

void Mesh::PrintFaceInfo() {
    for (Face f : faces) {
        std::cout << "Index\t" << f.GetIndex() << 
            "\nCenter\n";
        f.GetCenter().PrintCoordinates();
        std::cout << "Points\n";

        for (Point *p : f.GetPoints()) {
            p->PrintCoordinates();
        }
        std::cout << "\n";
    }
}

/*
void Mesh::GenerateFaces() {
    int index = 0;
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            Cell *cell = array[i][j];
            std::vector<Point*> p = cell->points;
            Face f1(p[0], p[1], p[2], p[3], index++); //front face
            Face f2(p[4], p[5], p[6], p[7], index++); //back face
            Face f3(p[1], p[2], p[5], p[6], index++); //right face
            Face f4(p[2], p[3], p[6], p[7], index++); //top face

            //Add faces to face vector
            faces.push_back(f1);
            faces.push_back(f2);
            faces.push_back(f3);
            faces.push_back(f4);

//            std::vector<int> vec = faces[f1.GetIndex()].GetCellIndices();
 //           std::cout << vec.size() << std::endl;

            //Add cell faces to cell
            cell->AddFace(&faces[f1.GetIndex()], f1.GetIndex());
            cell->AddFace(&faces[f2.GetIndex()], f2.GetIndex());
            cell->AddFace(&faces[f3.GetIndex()], f3.GetIndex());
            cell->AddFace(&faces[f4.GetIndex()], f4.GetIndex());

            //Add left face and cell face if at beginning of row
            if (i == 0) {
                Face f5(p[0], p[3], p[4], p[7], index++); //left face
                faces.push_back(f5);
                cell->AddFace(&faces[f5.GetIndex()], f5.GetIndex());
            }

            //Add bottom face if at bottom of column
            if (j == 0) {
                Face f6(p[0], p[1], p[4], p[5], index++); //bottom face
                faces.push_back(f6);
                cell->AddFace(&faces[f6.GetIndex()], f6.GetIndex());
            }

            //Add left cell face to next cell in row if it exists
            if (i < nx - 1) {
                Cell *cellRight = array[i+1][j];
                cellRight->AddFace(&faces[f3.GetIndex()], f3.GetIndex());
            }

            //Add bottom cell face to next cell up in column if it exists
            if (j < ny - 1) {
                Cell *cellUp = array[i][j+1];
                cellUp->AddFace(&faces[f4.GetIndex()], f4.GetIndex());
            }
        }
    }
}
*/

void Mesh::GenerateFacesBlock() {
    int index = 0;
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            Cell *cell = blockCellArray[i][j];
            std::vector<Point*> &p = cell->points;
            Face f1(p[0], p[1], p[2], p[3], index++); //front face
            Face f2(p[5], p[4], p[7], p[6], index++); //back face
            Face f3(p[1], p[5], p[6], p[2], index++); //right face
            Face f4(p[3], p[2], p[6], p[7], index++); //top face

            //Add faces to face vector
            faces.push_back(f1);
            faces.push_back(f2);
            faces.push_back(f3);
            faces.push_back(f4);

            //Add cell indices to faces
            faces[f1.GetIndex()].AddCellIndex(cell->index);
            faces[f2.GetIndex()].AddCellIndex(cell->index);
            faces[f3.GetIndex()].AddCellIndex(cell->index);
            faces[f4.GetIndex()].AddCellIndex(cell->index);

            cell->AddFace(f3.GetIndex()); //right
            cell->AddFace(f4.GetIndex()); //top
            cell->AddFace(f1.GetIndex()); //front
            cell->AddFace(f2.GetIndex()); //back

            if (i == 0) {
                Face f5(p[4], p[0], p[3], p[7], index++); //left face
                faces.push_back(f5);
                faces[f5.GetIndex()].AddCellIndex(cell->index);
                cell->AddFace(f5.GetIndex());
            }

            //Add left cell face to next cell in row if it exists
            if (i < nx - 1) {
                Cell *cellRight = blockCellArray[i+1][j];
                faces[f3.GetIndex()].AddCellIndex(cellRight->index);
                cellRight->AddFace(f3.GetIndex());
            }

            //Add bottom face if at bottom of column
            if (j == 0) {
                Face f6(p[4], p[5], p[1], p[0], index++); //bottom face
                faces.push_back(f6);
                faces[f6.GetIndex()].AddCellIndex(cell->index);
                cell->AddFace(f6.GetIndex());
            }

            //Add bottom cell face to next cell up in column if it exists
            if (j < ny - 1) {
                Cell *cellUp = blockCellArray[i][j+1];
                faces[f4.GetIndex()].AddCellIndex(cellUp->index);
                cellUp->AddFace(f4.GetIndex());
            }

        }
    }
}

void Mesh::GenerateCellFaces() {
    for (size_t i = 0; i < cellCount; i++) {
        Cell &c = cellList[i];
        std::vector<int> faceIndices = c.GetFaceIndices();
        for (int index : faceIndices) {
            Face &f = faces[index];
            c.AddCellFace(f, index);
        }
    }
}

void Mesh::CollectBoundaryFaces() {
    for (Face &face : faces) {
        if (face.GetCellIndices().size() == 1) {
            boundaryFaces.push_back(&face);
        }
        else internalFaces.push_back(&face);
    }
}

std::vector<Face*>& Mesh::GetBoundaryFaces() {
    return boundaryFaces;
}

std::vector<Face*>& Mesh::GetInternalFaces() {
    return internalFaces;
}

size_t Mesh::GetRowCount() {
    return nx;
}

size_t Mesh::GetColumnCount() {
    return ny;
}

Cell ***Mesh::GetArray() {
    return array;
}

Cell ***Mesh::GetBlockArray() {
    return blockCellArray;
}

std::vector<Face>& Mesh::GetFaceList() {
    return faces;
}

void Mesh::CalculateFluxes(const double &rho) {
    for (size_t i = 0; i < cellCount; i++) {
        Cell &c = cellList[i];
    }
}

void Mesh::StaggerToCell(double **u, double **v, double **p) {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            double newU = (u[i][j] + u[i+1][j]) / 2;
            double newV = (v[i][j] + v[i][j+1]) / 2;
            double newP = p[i][j];
            blockCellArray[i][j]->SetVelocity(Vector3D(newU, newV, 0));
            blockCellArray[i][j]->SetPressure(newP);
        }
    }
}
