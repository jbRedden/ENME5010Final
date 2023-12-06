#include "Peqn.h"
#include "StaggeredGrid.h"

Peqn::Peqn() {}
Peqn::Peqn(StaggeredGrid *grid, double timestep) 
    : grid(grid), nx(grid->GetPNX()), ny(grid->GetPNY()), dt(timestep) {
        AllocateMemory();
    }

Peqn::~Peqn() {
    for (size_t i = 0; i < nx; i++) {
        delete[] p[i];
        delete[] pStar[i];
    }

    for (size_t i = 0; i < nx-2; i++) {
        delete[] A[i];
        delete[] B[i];
        delete[] C[i];
        delete[] D[i];
        delete[] E[i];
        delete[] F[i];
    }

    delete[] p;
    delete[] pStar;
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;
    delete[] E;
    delete[] F;
}

void Peqn::AllocateMemory() {
    p = new double*[nx];
    pStar = new double*[nx];
    A = new double*[nx-2];
    B = new double*[nx-2];
    C = new double*[nx-2];
    D = new double*[nx-2];
    E = new double*[nx-2];
    F = new double*[nx-2];

    for (size_t i = 0; i < nx; i++) {
        p[i] = new double[ny];
        pStar[i] = new double[ny];
    }

    for (size_t i = 0; i < nx-2; i++) {
        A[i] = new double[ny-2];
        B[i] = new double[ny-2];
        C[i] = new double[ny-2];
        D[i] = new double[ny-2];
        E[i] = new double[ny-2];
        F[i] = new double[ny-2];
    }
}

void Peqn::InitField(double p0) {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            p[i][j] = p0;
            pStar[i][j] = p0;
        }
    }
}

void Peqn::EnforceBC() {
    for (size_t i = 0; i < nx; i++) {
        pStar[i][0] = pStar[i][1];
        pStar[i][ny-1] = pStar[i][ny-1];
    }

    for (size_t j = 0; j < ny; j++) {
        pStar[0][j] = pStar[1][j];
        pStar[nx-1][j] = pStar[nx-2][j];
    }
}

void Peqn::Overwrite() {
    for (size_t j = ny-1; j > 0; j--) {
        for (size_t i = 0; i < nx; i++) {
            p[i][j] = pStar[i][j];
        }
    }
}

void Peqn::Print() {
    for (size_t j = ny-1; j > 0; j--) {
        for (size_t i = 0; i < nx; i++) {
            std::cout << p[i][j] << "\t";
        }
            std::cout << "\n";
    }
}

double** Peqn::GetP() {
    return p;
}

double** Peqn::GetPStar() {
    return pStar;
}

void Peqn::UpdateCoefficients(Ueqn &uEqn, Veqn &vEqn) {
    pCell **cells = grid->GetPCells();
    double **vStar = vEqn.GetVStar();
    double **uStar = uEqn.GetUStar();
    for (size_t j = 1; j < ny-1; j++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t i2 = i - 1;
            size_t j2 = j - 1;
            const pCell &cell = cells[i][j];
            double dx = cell.GetDX();
            double dy = cell.GetDY();
            //Note time steps are halved for ADI, no pressure terms yet. Need to add p* terms for after pressure is corrected
            A[i-1][j-1] = 0.5 * dy * dy * (1 / uEqn.A[i2][j2] + 1 / uEqn.A[i2+1][j2]) + 0.5 * dx * dx * (1 / vEqn.A[i2][j2] + 1 / vEqn.A[i2][j2+1]);
            B[i-1][j-1] = -0.5 * dy * dy / uEqn.A[i2+1][j2];
            C[i-1][j-1] = -0.5 * dy * dy / uEqn.A[i2][j2];
            D[i-1][j-1] = -0.5 * dx * dx / vEqn.A[i2][j2+1];
            E[i-1][j-1] = -0.5 * dx * dx / vEqn.A[i2][j2];
            F[i-1][j-1] = dy * (uStar[i][j] - uStar[i+1][j]) + dx * (vStar[i][j] - vStar[i][j+1]);
        }
    }
}

void Peqn::Solve(Ueqn &uEqn, Veqn &vEqn) {
    pCell **cells = grid->GetPCells();
    size_t nx2 = nx - 2;
    size_t ny2 = ny - 2;
    double pCorr[nx2][ny2];
    double A2[nx-2][ny-2];
    double B2[nx-2][ny-2];
    double C2[nx-2][ny-2];
    double D2[nx-2][ny-2];
    double E2[nx-2][ny-2];
    double F2[nx-2][ny-2];
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            double dy = cells[i+1][j+1].GetDY();
            A2[i][j] = A[i][j];
            B2[i][j] = B[i][j];
            C2[i][j] = C[i][j];
            D2[i][j] = D[i][j];
            E2[i][j] = E[i][j];
            F2[i][j] = F[i][j];
        }
    }

    for (size_t j = 0; j < ny-2; j++) {
        B2[nx2-1][j] = 0;
        A2[nx2-1][j] += B[nx2-1][j];
        //F2[nx2-1][j] -= B[nx2-1][j] * rightBC;//uStar[nx-1][j+1];
        C2[0][j] = 0;
        A2[0][j] += C[0][j];
        //F2[0][j] -= C[0][j] * leftBC;//uStar[0][j+1];
    }

    //first half step
    for (size_t j = 0; j < ny2; j++) {
        double a[nx2];
        double b[nx2];
        double c[nx2];
        double d[nx2];

        for (size_t i = 0; i < nx2; i++) {
            a[i] = C2[i][j];
            b[i] = A2[i][j];
            c[i] = B2[i][j];
            d[i] = F2[i][j];
        }

        double f[nx2];
        ThomasAlg(&a[0], &b[0], &c[0], &d[0], &f[0], nx2);
        for (size_t i = 0; i < nx2; i++) pCorr[i][j] = f[i];
    }

    //second half step
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            double dy = cells[i+1][j+1].GetDY();
            A2[i][j] = A[i][j];
            B2[i][j] = B[i][j];
            C2[i][j] = C[i][j];
            D2[i][j] = D[i][j];
            E2[i][j] = E[i][j];
            F2[i][j] = F[i][j]; 
        }
    }

    for (size_t i = 0; i < nx-2; i++) {
        D2[i][ny2-1] = 0;
        A2[i][ny2-1] += D[i][ny2-1];
        //F2[i][ny2-1] -= D[i][ny2-1] * topBC;
        E2[i][0] = 0;
        A2[i][0] += E[i][ny2-1];
        //F2[i][0] -= E[i][ny2-1] * botBC;
    }

    for (size_t i = 0; i < nx2; i++) {
        double a[nx2];
        double b[nx2];
        double c[nx2];
        double d[nx2];

        for (size_t j = 0; j < ny2; j++) {
            a[i] = E2[i][j];
            b[i] = A2[i][j];
            c[i] = D2[i][j];
            d[i] = F2[i][j];
        }

        double f[nx2];
        ThomasAlg(&a[0], &b[0], &c[0], &d[0], &f[0], nx2);
        for (size_t j = 0; j < ny2; j++) pCorr[i][j] = f[i];
    }

    /*
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            std::cout << pCorr[i][j] << "\t";
            pStar[i+1][j+1] += pCorr[i][j];
        }
        std::cout << "\n";
    }
    */

    uCell **uCells = grid->GetUCells();
    for (size_t j = 0; j < uEqn.ny-2; j++) {
        for (size_t i = 0; i < uEqn.nx-2; i++) {
            uEqn.uStar[i+1][j+1] -= 0.5 * uCells[i+1][j+1].GetDY() * (pCorr[i][j]-pCorr[i-1][j]) / uEqn.A[i][j];
        }
    }

    vCell **vCells = grid->GetVCells();
    for (size_t j = 0; j < vEqn.ny-2; j++) {
        for (size_t i = 0; i < vEqn.nx-2; i++) {
            vEqn.vStar[i+1][j+1] -= 0.5 * vCells[i+1][j+1].GetDX() * (pCorr[i][j]-pCorr[i][j-1]) / vEqn.A[i][j];
        }
    }
}

void Peqn::ThomasAlg(double *a, double *b, double *c, double *d, double *f, size_t n) {
    for (size_t i = 1; i < n; i++) {
        double w = a[i] / b[i-1];
        b[i] -= w * c[i-1];
        d[i] -= w * d[i-1];
    }

    f[n-1] = d[n-1] / b[n-1];

    for (size_t i = n-1; i > 0; i--) {
        f[i-1] = (d[i-1] - c[i-1] * f[i]) / b[i-1];
    }
}
