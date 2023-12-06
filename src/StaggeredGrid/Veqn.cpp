#include "Veqn.h"
#include "Ueqn.h"
#include "Peqn.h"
#include "StaggeredGrid.h"

Veqn::Veqn() {}
Veqn::Veqn(StaggeredGrid *grid, double rho, double mu, 
        double timestep, double leftBC, double rightBC, 
        double topBC, double botBC) 
    : grid(grid), nx(grid->GetVNX()), ny(grid->GetVNY()), 
    rho(rho), mu(mu), dt(timestep), leftBC(leftBC),
    rightBC(rightBC), topBC(topBC), botBC(botBC) {
        AllocateMemory();
    }

Veqn::~Veqn() {
    for (size_t i = 0; i < nx; i++) {
        delete[] v[i];
        delete[] vStar[i];
    }

    for (size_t i = 0; i < nx; i++) {
        v[i] = new double[ny];
        vStar[i] = new double[ny];
    }

    for (size_t i = 0; i < nx-2; i++) {
        delete[] A[i];
        delete[] F[i];
        delete[] B[i];
        delete[] C[i];
        delete[] D[i];
        delete[] E[i];
    }


    delete[] v;
    delete[] vStar;
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;
    delete[] E;
    delete[] F;
}

void Veqn::AllocateMemory() {
    v = new double*[nx];
    vStar = new double*[nx];
    A = new double*[nx-2];
    B = new double*[nx-2];
    C = new double*[nx-2];
    D = new double*[nx-2];
    E = new double*[nx-2];
    F = new double*[nx-2];

    for (size_t i = 0; i < nx; i++) {
        v[i] = new double[ny];
        vStar[i] = new double[ny];
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

void Veqn::InitField(double v0) {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            v[i][j] = v0;
            vStar[i][j] = v0;
        }
    }
}

void Veqn::EnforceBC() {
    for (size_t j = 0; j < ny; j++) {
        v[0][j] = leftBC;
        v[nx-1][j] = rightBC;
        vStar[0][j] = leftBC;
        vStar[nx-1][j] = rightBC;
    }

    for (size_t i = 0; i < nx; i++) {
        v[i][0] = botBC;
        v[i][ny-1] = topBC;
        vStar[i][0] = botBC;
        vStar[i][ny-1] = topBC;
    }
}

void Veqn::Overwrite() {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            v[i][j] = vStar[i][j];
        }
    }
}

void Veqn::Print() {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            std::cout << v[i][j] << "\t";
        }
            std::cout << "\n";
    }
}

void Veqn::PrintVStar() {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            std::cout << vStar[i][j] << "\t";
        }
            std::cout << "\n";
    }
}

double** Veqn::GetV() {
    return v;
}

double** Veqn::GetVStar() {
    return vStar;
}

void Veqn::UpdateCoefficients(Ueqn &uEqn, Peqn &pEqn) {
    vCell **cells = grid->GetVCells();
    double **u = uEqn.GetU();
    double **p = pEqn.GetP();
    for (size_t j = 1; j < ny-1; j++) {
        for (size_t i = 1; i < nx-1; i++) {
            const vCell &cell = cells[i][j];
            double le = cell.LambdaE(cells, nx, ny);
            double lw = cell.LambdaW(cells, nx, ny);
            double ln = cell.LambdaN(cells, nx, ny);
            double ls = cell.LambdaS(cells, nx, ny);
            double delE = cell.DelE(cells, nx, ny);
            double delW = cell.DelW(cells, nx, ny);
            double delN = cell.DelN(cells, nx, ny);
            double delS = cell.DelS(cells, nx, ny);
            double dx = cell.GetDX();
            double dy = cell.GetDY();
            A[i-1][j-1] = 2 * rho * dx * dy / dt + rho * dx * ((1-ln)*(1-ln) * v[i][j] + ln * (1-ln) * v[i][j+1]
                    - (1-ls)*(1-ls) * v[i][j] - ls * (1-ls) * v[i][j-1]) + 0.25 * rho * dy * ((u[i+1][j] + u[i+1][j-1]) * (1-le)
                    - (u[i][j] + u[i][j-1]) * (1-lw)) + 0.5 * mu * dy * (delE - delW) + 0.5 * mu * dx * (delN - delS);
            B[i-1][j-1] = 0.25 * rho * dy * le * (u[i+1][j] + u[i+1][j-1]) - 0.5 * mu * dy * delE;
            C[i-1][j-1] = -0.25 * rho * dy * lw * (u[i][j] + u[i][j-1]) + 0.5 * mu * dy * delW;

            D[i-1][j-1] = rho * dx * (ln * (1-ln) * v[i][j] + ln * ln * v[i][j+1]) - 0.5 * mu * dx * delN;
            E[i-1][j-1] = -rho * dx * (ls * (1-ls) * v[i][j] + ls * ls * v[i][j-1]) + 0.5 * mu * dx * delS;
            F[i-1][j-1] = 2 * rho * dx * dy / dt * v[i][j]
                - 0.5 * dx * (p[i][j] - p[i][j-1]) 
                + 0.25 * rho * dy * ((u[i][j] + u[i][j-1]) * ((1-lw) * v[i][j] + lw * v[i-1][j])
                        - (u[i+1][j] + u[i+1][j-1])*((1-le)*v[i][j] + le * v[i+1][j]))
                + 0.5 * mu * dy * (delE * (v[i+1][j] - v[i][j]) - delW * (v[i-1][j] - v[i][j]))
                + 0.5 * mu * dx * (delN * (v[i][j+1] - v[i][j]) - delS * (v[i][j-1] - v[i][j]));
        }
    }
/*
    std::ofstream f;
    f.open("vEqn.csv");
    f << "\nA\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << A[i][j] << ",";
        }
        f << "\n";
    }

    f << "\nB\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << B[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nC\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << C[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nD\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << D[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nE\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << E[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nF\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << F[i][j] << ",";
        }
        f << "\n";
    }

    f.close();
    */
}

void Veqn::Solve(Peqn &pEqn) {
    uCell **cells = grid->GetUCells();
    double **pStar = pEqn.GetPStar();
    size_t nx2 = nx - 2;
    size_t ny2 = ny - 2;
    double A2[nx-2][ny-2];
    double B2[nx-2][ny-2];
    double C2[nx-2][ny-2];
    double D2[nx-2][ny-2];
    double E2[nx-2][ny-2];
    double F2[nx-2][ny-2];

    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            double dx = cells[i+1][j+1].GetDX();
            A2[i][j] = A[i][j];
            B2[i][j] = B[i][j];
            C2[i][j] = C[i][j];
            D2[i][j] = D[i][j];
            E2[i][j] = E[i][j];
            F2[i][j] = F[i][j] - D[i][j] * vStar[i+1][j+2] - E[i][j] * vStar[i+1][j]
                - 0.5 * dx * (pStar[i+1][j+1] - pStar[i+1][j]);
        }
    }

    for (size_t j = 0; j < ny-2; j++) {
        B2[nx2-1][j] = 0;
        F2[nx2-1][j] -= B[nx2-1][j] * rightBC;//uStar[nx-1][j+1];
        C2[0][j] = 0;
        F2[0][j] -= C[0][j] * leftBC;//uStar[0][j+1];
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
        for (size_t i = 0; i < nx2; i++) vStar[i+1][j+1] = f[i];
    }

    //second half step
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            double dx = cells[i+1][j+1].GetDX();
            A2[i][j] = A[i][j];
            B2[i][j] = B[i][j];
            C2[i][j] = C[i][j];
            D2[i][j] = D[i][j];
            E2[i][j] = E[i][j];
            F2[i][j] = F[i][j] - B[i][j] * vStar[i+2][j+1] - C[i][j] * vStar[i][j+1]
                - 0.5 * dx * (pStar[i+1][j+1] - pStar[i+1][j]);
        }
    }

    for (size_t i = 0; i < nx-2; i++) {
        D2[i][ny2-1] = 0;
        F2[i][ny2-1] -= D[i][ny2-1] * topBC;
        E2[i][0] = 0;
        F2[i][0] -= E[i][ny2-1] * botBC;
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
        for (size_t j = 0; j < ny2; j++) vStar[i+1][j+1] = f[i];
    }

    std::ofstream f;
    f.open("testv.csv");

    f << "A\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << A[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nB\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << B[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nC\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << C[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nD\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << D[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nE\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << E[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nF\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << F[i][j] << ",";
        }
        f << "\n";
    }


    f << "A2\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << A2[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nB2\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << B2[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nC2\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << C2[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nD2\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << D2[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nE2\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << E2[i][j] << ",";
        }
        f << "\n";
    }
    f << "\nF2\n";
    for (size_t j = 0; j < ny-2; j++) {
        for (size_t i = 0; i < nx-2; i++) {
            f << F2[i][j] << ",";
        }
        f << "\n";
    }
    f.close();
}

void Veqn::ThomasAlg(double *a, double *b, double *c, double *d, double *f, size_t n) {
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
