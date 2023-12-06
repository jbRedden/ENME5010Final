#include "StaggeredGrid.h"
#include "StaggeredCell.h"
#include "pCell.h"

StaggeredGrid::StaggeredGrid() {}

StaggeredGrid::StaggeredGrid(Mesh *mesh) 
    : mesh(mesh), pnx(mesh->GetNX()), pny(mesh->GetNY()), 
    unx(pnx + 1), uny(pny), vnx(pnx), vny(pny + 1) {
        AllocateMemory();
        PopulateGrids();
    }

StaggeredGrid::~StaggeredGrid() {
    for (size_t i = 0; i < pnx; i++) {
        delete[] pCells[i];
    }
    for (size_t i = 0; i < unx; i++) {
        delete[] uCells[i];
    }
    for (size_t i = 0; i < vnx; i++) {
        delete[] vCells[i];
    }

    delete[] pCells;
    delete[] uCells;
    delete[] vCells;
}

void StaggeredGrid::AllocateMemory() {
    pCells = new pCell*[pnx];
    uCells = new uCell*[unx];
    vCells = new vCell*[vnx];

    for (size_t i = 0; i < pnx; i++) {
        pCells[i] = new pCell[pny];
    }

    for (size_t i = 0; i < unx; i++) {
        uCells[i] = new uCell[uny];
    }

    for (size_t i = 0; i < vnx; i++) {
        vCells[i] = new vCell[vny];
    }
}

void StaggeredGrid::PopulateGrids() {
    Cell ***cArray = mesh->GetBlockArray();
    Cell c1 = *cArray[0][0];
    for (size_t j = 0; j < pny; j++) {
        for (size_t i = 0; i < pnx; i++) {
            Cell &c = *cArray[i][j];
            pCells[i][j] = pCell(c.GetCenter(), c.GetDX(), c.GetDY(), i, j);

            if (i < pnx - 1) {
                uCells[i+1][j] = uCell(c.eastFace.GetCenter(), cArray[i+1][j]->GetCenter().GetX() - c.GetCenter().GetX(), c.GetDY(), i+1, j);
            }
            else if (i == pnx - 1) uCells[i+1][j] = uCell(c.eastFace.GetCenter(), i+1, j);

            if (j < pny - 1) {
                vCells[i][j+1] = vCell(c.northFace.GetCenter(), c.GetDX(), cArray[i][j+1]->GetCenter().GetY() - c.GetCenter().GetY(), i, j+1);
            }
            else if (j == pny - 1) {
                vCells[i][j+1] = vCell(c.northFace.GetCenter(), i, j+1);
            }

            if (i == 0) {
                uCells[0][j] = uCell(c.westFace.GetCenter(), 0, j);
            }
            if (j == 0) {
                vCells[i][0] = vCell(c.southFace.GetCenter(), i, 0);
            }
        }
    }
}

/*
void StaggeredGrid::Testing() {
    for (size_t j = 1; j < uny - 1; j++) {
        for (size_t i = 1; i < unx - 1; i++) {
            uCell &u = uCells[i][j];
            u.GetCellCenter().PrintCoordinates();
            double d = u.LambdaN(reinterpret_cast<StaggeredCell **>(uCells), unx, uny);
            double d2 = u.LambdaS(reinterpret_cast<StaggeredCell **>(uCells), unx, uny);
            std::cout << d << "\t" << d2 << std::endl;
        }
    }
}
*/

void StaggeredGrid::PrintPGrid() const {
    std::cout << "Pressure Grid:\n";
    for (size_t j = 0; j < pny; j++) {
        for (size_t i = 0; i < pnx; i++) {
            pCells[i][j].GetCellCenter().PrintCoordinates();
        }
    }
}

void StaggeredGrid::PrintUGrid() const {
    std::cout << "U Grid:\n";
    for (size_t j = 0; j < uny; j++) {
        for (size_t i = 0; i < unx; i++) {
            uCells[i][j].GetCellCenter().PrintCoordinates();
        }
    }
}

void StaggeredGrid::PrintVGrid() const {
    std::cout << "V Grid:\n";
    for (size_t j = 0; j < vny; j++) {
        for (size_t i = 0; i < vnx; i++) {
            vCells[i][j].GetCellCenter().PrintCoordinates();
        }
    }
}

void StaggeredGrid::WriteCSV() const {
    std::ofstream file;
    file.open("p.csv");
    file << "x,y,z\n";
    for (size_t j = 0; j < pny; j++) {
        for (size_t i = 0; i < pnx; i++) {
            const Point &p = pCells[i][j].GetCellCenter();
            file << p.GetX() << ',' << p.GetY() << ',' << p.GetZ() << std::endl; 
        }
    }
    file.close();

    file.open("u.csv");
    file << "x,y,z\n";
    for (size_t j = 0; j < uny; j++) {
        for (size_t i = 0; i < unx; i++) {
            const Point &p = uCells[i][j].GetCellCenter();
            file << p.GetX() << ',' << p.GetY() << ',' << p.GetZ() << std::endl; 
        }
    }
    file.close();

    file.open("v.csv");
    file << "x,y,z\n";
    for (size_t j = 0; j < vny; j++) {
        for (size_t i = 0; i < vnx; i++) {
            const Point &p = vCells[i][j].GetCellCenter();
            file << p.GetX() << ',' << p.GetY() << ',' << p.GetZ() << std::endl; 
        }
    }
    file.close();
}

double StaggeredGrid::GetUNX() {
    return unx;
}

double StaggeredGrid::GetUNY() {
    return uny;
}

uCell** StaggeredGrid::GetUCells() {
    return uCells;
}

double StaggeredGrid::GetVNX() {
    return vnx;
}

double StaggeredGrid::GetVNY() {
    return vny;
}

vCell** StaggeredGrid::GetVCells() {
    return vCells;
}

double StaggeredGrid::GetPNX() {
    return pnx;
}

double StaggeredGrid::GetPNY() {
    return pny;
}

pCell** StaggeredGrid::GetPCells() {
    return pCells;
}
