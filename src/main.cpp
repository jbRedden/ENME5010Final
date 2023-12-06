#include "Geometry/Geometry.h"
#include "Utilities/Utilities.h"
#include "StaggeredGrid/Staggered.h"
#include <string>

int main(void) {

    size_t nx = 150;
    size_t ny = 150;
    size_t nz = 1;
    double xMin = 0;
    double xMax = 1;
    double yMin = 0;
    double yMax = 1;
    double zMin = 0;
    double zMax = 0.001;
    double dx = (xMax - xMin) / nx;
    double dy = (yMax - yMin) / ny;
    double dz = (zMax - zMin) / nz;
    Block block(xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz);
    Mesh mesh(&block);


    size_t cCount = mesh.GetCellCount();
    Cell *cells = mesh.GetCellList();

    StaggeredGrid grid(&mesh);

    //Currenly only works with density of 1 *********************************************************
    double rho = 1;
    double mu = 0.01;
    double dt = 0.001;

    Ueqn2 uEqn2(&grid, rho, mu, dt, 0, 0, 1, 0);
    Veqn2 vEqn(&grid, rho, mu, dt, 0, 0, 0, 0);
    Peqn2 pEqn(&grid, dt);

    pEqn.InitField();
    pEqn.EnforceBC();
    vEqn.InitField();
    vEqn.EnforceBC();
    uEqn2.InitField();
    uEqn2.EnforceBC();

    uEqn2.UpdateCoefficients(vEqn, pEqn);
    vEqn.UpdateCoefficients(uEqn2, pEqn);


    for (size_t outter = 0; outter < 1000000; outter++) {
        for (size_t iter = 0; iter < 50000; iter++) {
            uEqn2.Solve(pEqn);
            vEqn.Solve(pEqn);
            pEqn.UpdateCoefficients(uEqn2, vEqn);
            pEqn.Solve(uEqn2, vEqn);
            pEqn.EnforceBC();
            uEqn2.EnforceBC();
            vEqn.EnforceBC();

            if (pEqn.pCorrMax < 1e-6) break;
            std::cout << pEqn.pCorrMax << "\t" << pEqn.uCorrMax << "\t" << pEqn.vCorrMax << std::endl;
            if (pEqn.uCorrMax < 1e-4 && pEqn.vCorrMax < 1e-4) break;
        }
        if (pEqn.pCorrMax < 1e-6 || (pEqn.uCorrMax < 1e-10 && pEqn.vCorrMax < 1e-10)) break;
        pEqn.Overwrite();
        uEqn2.Overwrite();
        vEqn.Overwrite();
        uEqn2.UpdateCoefficients(vEqn, pEqn);
        vEqn.UpdateCoefficients(uEqn2, pEqn);
    }

    pEqn.EnforceBC();
    uEqn2.EnforceBC();
    vEqn.EnforceBC();


    mesh.StaggerToCell(uEqn2.GetU(), vEqn.GetV(), pEqn.GetP());
    IO::WriteVTU(mesh, "points");

    return 0;
}
