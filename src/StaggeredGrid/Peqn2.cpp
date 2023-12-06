#include "Peqn2.h"
#include "StaggeredGrid.h"

Peqn2::Peqn2() {}
Peqn2::Peqn2(StaggeredGrid *grid, double timestep) 
    : grid(grid), nx(grid->GetPNX()), ny(grid->GetPNY()), dt(timestep) {
        AllocateMemory();
    }

Peqn2::~Peqn2() {
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

void Peqn2::AllocateMemory() {
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

void Peqn2::InitField(double p0) {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            p[i][j] = p0;
            pStar[i][j] = p0;
        }
    }
}

void Peqn2::EnforceBC() {

    for (size_t j = 0; j < ny; j++) {
        pStar[0][j] = pStar[1][j];
        pStar[nx-1][j] = pStar[nx-2][j];
        p[0][j] = p[1][j];
        p[nx-1][j] = p[nx-2][j];
    }
    for (size_t i = 0; i < nx; i++) {
        pStar[i][0] = pStar[i][1];
        pStar[i][ny-1] = pStar[i][ny-2];
        p[i][0] = p[i][1];
        p[i][ny-1] = p[i][ny-2];
    }
}

void Peqn2::Overwrite() {
    for (size_t j = ny-1; j > 0; j--) {
        for (size_t i = 0; i < nx; i++) {
            p[i][j] = pStar[i][j];
        }
    }
}

void Peqn2::Print() {
    for (size_t j = ny-1; j > 0; j--) {
        for (size_t i = 0; i < nx; i++) {
            std::cout << p[i][j] << "\t";
        }
            std::cout << "\n";
    }
}

double** Peqn2::GetP() {
    return p;
}

double** Peqn2::GetPStar() {
    return pStar;
}

void Peqn2::UpdateCoefficients(Ueqn2 &uEqn, Veqn2 &vEqn) {
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

void Peqn2::Solve(Ueqn2 &uEqn, Veqn2 &vEqn) {
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
            d[i] = F2[i][j];// - D2[i][j]*pCorr[i][j] - E2[i][j]*pCorr[i][j];
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

    //I think there was an error here, was E[i][ny2-1]
    for (size_t i = 0; i < nx-2; i++) {
        D2[i][ny2-1] = 0;
        A2[i][ny2-1] += D[i][ny2-1];

        E2[i][0] = 0;
        A2[i][0] += E[i][0];
    }

    for (size_t i = 0; i < nx2; i++) {
        double a[ny2];
        double b[ny2];
        double c[ny2];
        double d[ny2];

        for (size_t j = 0; j < ny2; j++) {
            a[j] = E2[i][j];
            b[j] = A2[i][j];
            c[j] = D2[i][j];
            d[j] = F2[i][j];
//            if (j != 0) d[j] -= E[i][j] * pCorr[i][j-1];
//            if (j != ny2-1) d[j] -= D[i][j] * pCorr[i][j+1];
            if (i == 0) b[j] += C[i][j];
            else if (i != 0) d[j] -= C[i][j] * pCorr[i-1][j];
            if (i == nx - 1) b[j] += B[i][j];
            else if (i != nx - 1) d[j] -= B[i][j] * pCorr[i+1][j];
            
        }

        double f[ny2];
        ThomasAlg(&a[0], &b[0], &c[0], &d[0], &f[0], ny2);
        for (size_t j = 0; j < ny2; j++) pCorr[i][j] = f[j];
    }

        double relax = 0.3;
        pCorrMax = 0;
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            double pC = pCorr[i][j];
            if (std::abs(pC) > pCorrMax) pCorrMax = std::abs(pC);
            pC *= relax;
            pC = trunc(pC * 1e14) / 1e14;
            pStar[i+1][j+1] += pC;
        }
    }

    uCell **uCells = grid->GetUCells();
    uCorrMax = 0;
    for (size_t j = 0; j < uEqn.ny-2; j++) {
        for (size_t i = 1; i < uEqn.nx-3; i++) {
            double dx = cells[i+1][j+1].GetDX();
            double dy = cells[i+1][j+1].GetDY();
            double uCorr = -0.5 * uCells[i+1][j+1].GetDY() * (pCorr[i][j]-pCorr[i-1][j]) / (uEqn.A[i][j] - dx * dy / dt);
            if (std::abs(uCorr) > uCorrMax) uCorrMax = std::abs(uCorr);
            uCorr *= relax;
            uCorr = trunc(uCorr * 1e14) / 1e14;
            uEqn.uStar[i+1][j+1] += uCorr;
        }
    }

    vCell **vCells = grid->GetVCells();
    vCorrMax = 0;
    for (size_t j = 1; j < vEqn.ny-3; j++) {
        for (size_t i = 0; i < vEqn.nx-2; i++) {
            double dx = cells[i+1][j+1].GetDX();
            double dy = cells[i+1][j+1].GetDY();
            double vCorr = -0.5 * vCells[i+1][j+1].GetDX() * (pCorr[i][j]-pCorr[i][j-1]) / (vEqn.A[i][j] - dx * dy / dt);
            if (vCorr > vCorrMax) vCorrMax = std::abs(vCorr);
            vCorr *= relax;
            vCorr = trunc(vCorr * 1e14) / 1e14;
            vEqn.vStar[i+1][j+1] += vCorr;
        }
    }
}

void Peqn2::ThomasAlg(double *a, double *b, double *c, double *d, double *f, size_t n) {
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
