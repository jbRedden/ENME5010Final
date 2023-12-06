#include "Ueqn.h"
#include "StaggeredGrid.h"
#include "Veqn.h"

Ueqn::Ueqn() {}
Ueqn::Ueqn(StaggeredGrid *grid, double rho, double mu, 
        double timestep, double leftBC, double rightBC, 
        double topBC, double botBC) 
    : grid(grid), nx(grid->GetUNX()), ny(grid->GetUNY()), 
    rho(rho), mu(mu), dt(timestep), leftBC(leftBC),
    rightBC(rightBC), topBC(topBC), botBC(botBC) {
        AllocateMemory();
    }

Ueqn::~Ueqn() {
    for (size_t i = 0; i < nx; i++) {
        delete[] u[i];
        delete[] uStar[i];
    }

    for (size_t i = 0; i < nx-2; i++) {
        delete[] A[i];
        delete[] B[i];
        delete[] C[i];
        delete[] D[i];
        delete[] E[i];
        delete[] F[i];
    }

    delete[] u;
    delete[] uStar;
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;
    delete[] E;
    delete[] F;
}

void Ueqn::AllocateMemory() {
    u = new double*[nx];
    uStar = new double*[nx];
    A = new double*[nx-2];
    B = new double*[nx-2];
    C = new double*[nx-2];
    D = new double*[nx-2];
    E = new double*[nx-2];
    F = new double*[nx-2];

    for (size_t i = 0; i < nx; i++) {
        u[i] = new double[ny];
        uStar[i] = new double[ny];
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

void Ueqn::InitField(double u0) {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            u[i][j] = u0;
            uStar[i][j] = u0;
        }
    }
}

void Ueqn::EnforceBC() {
    for (size_t j = 0; j < ny; j++) {
        u[0][j] = leftBC;
        u[nx-1][j] = rightBC;
        uStar[0][j] = leftBC;
        uStar[nx-1][j] = rightBC;
    }

    for (size_t i = 0; i < nx; i++) {
        u[i][0] = botBC;
        u[i][ny-1] = topBC;
        uStar[i][0] = botBC;
        uStar[i][ny-1] = topBC;
    }
}

void Ueqn::Print() {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            std::cout << u[i][j] << "\t";
        }
            std::cout << "\n";
    }
}

void Ueqn::Overwrite() {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            u[i][j] = uStar[i][j];
        }
    }
}

double** Ueqn::GetU() {
    return u;
}

double** Ueqn::GetUStar() {
    return uStar;
}

void Ueqn::UpdateForNextTimeStep() {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            u[i][j] = uStar[i][j];
        }
    }
}

//note - check delW and delS derivative definitions.
void Ueqn::UpdateCoefficients(Veqn &vEqn, Peqn &pEqn) {
    uCell **cells = grid->GetUCells();
    double **v = vEqn.GetV();
    double **p = pEqn.GetP();
    for (size_t j = 1; j < ny-1; j++) {
        for (size_t i = 1; i < nx-1; i++) {
            const uCell &cell = cells[i][j];
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
            //Note time steps are halved for ADI, no pressure terms yet. Need to add p* terms for after pressure is corrected
            A[i-1][j-1] = 2 * rho * dx * dy / dt 
                + rho * dy * ((1-le) * (1-le) * u[i][j] + le * (1-le) * u[i+1][j] 
                        - (1-lw) * (1-lw) * u[i][j] - lw * (1-lw) * u[i-1][j])
                + 0.25 * rho * dx * ((v[i-1][j+1] + v[i][j+1]) * (1 - ln)
                        - (v[i-1][j] + v[i][j]) * (1 - ls))
                + 0.5 * mu * (dy * (delE - delW) + dx * (delN - delS));
            B[i-1][j-1] = rho * dy * (le * (1-le) * u[i][j] + le * le * u[i+1][j]) - 0.5 * mu * dy * delE;
            C[i-1][j-1] = -rho * dy * (lw * (1-lw) * u[i][j] + lw * lw * u[i-1][j]) + 0.5 * mu * dy * delW;

            D[i-1][j-1] = 0.25 * rho * dx * ln * (v[i-1][j+1] + v[i][j+1]) - 0.5 * mu * dx * delN;
            E[i-1][j-1] = -0.25 * rho * dx * ls * (v[i-1][j] + v[i][j]) + 0.5 * mu * dx * delS;
            F[i-1][j-1] = 2 * rho * dx * dy / dt * u[i][j]
                + 0.25 * rho * dx *((v[i-1][j] + v[i][j]) * ((1 - ls) * u[i][j] + ls * u[i][j-1])
                        - (v[i-1][j+1] + v[i][j+1]) * ((1 - ln) * u[i][j] + ln * u[i][j+1]))
                + 0.5 * mu * dy * (delE * (u[i+1][j] - u[i][j]) - delW * (u[i-1][j] - u[i][j]))
                + 0.5 * mu * dx * (delN * (u[i][j+1] - u[i][j]) - delS * (u[i][j-1] - u[i][j]))
                - 0.5 * dy * (p[i][j] - p[i-1][j]);
        }
    }
}


void Ueqn::Solve(size_t j) {
    size_t row = nx - 2;
    double u1[nx];
    double uNew[row];
    double a1[row];
    double b1[row];
    double c1[row];
    double d1[row];

    for (size_t i = 0; i < nx; i++) {
        u1[i] = u[i][1];
    }

    for (size_t i = 0; i < nx - 2; i++) {
        b1[i] = A[i][0];
        d1[i] = F[i][0];
        if (i > 0) a1[i-1] = C[i-1][j-1];
        if (i < nx - 3) c1[i] = B[i][j-1];
    }
    a1[0] = 0;
    c1[row-1] = 0;
    d1[0] += a1[0] * u[0][j];
    d1[row-1] += c1[row-2] * u[nx-1][j];

    /*Forward Sweep*/
    for (size_t i = 1; i < row; i++) {
        if (b1[i-1] == 0) throw std::runtime_error("Error: division by zero.");
        double w = a1[i] / b1[i-1];
        b1[i] -= w * c1[i-1];
        d1[i] -= w * d1[i-1];
    }

    /*Backward Sweep*/
    u1[row-1] = d1[row-1] / b1[row-1];
    for (size_t i = row - 1; i > 0; i--) {
        u[i][j] = (d1[i] - c1[i] * u1[i+1]) / b1[i];
    }

    for (size_t i = 0; i < nx; i++) {
        std::cout << u[i][j] << "\t";
    }
    std::cout << "\n";
}

void Ueqn::Solv2(size_t j, Peqn &pEqn) {
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
            double dy = cells[i+1][j+1].GetDY();
            A2[i][j] = A[i][j];
            B2[i][j] = B[i][j];
            C2[i][j] = C[i][j];
            D2[i][j] = D[i][j];
            E2[i][j] = E[i][j];
            F2[i][j] = F[i][j] - D[i][j] * uStar[i+1][j+2] - E[i][j] * uStar[i+1][j]
                - 0.5 * dy * (pStar[i+1][j+1] - pStar[i][j+1]);
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
        for (size_t i = 0; i < nx2; i++) uStar[i+1][j+1] = f[i];
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
            F2[i][j] = F[i][j] - B[i][j] * uStar[i+2][j+1] - C[i][j] * uStar[i][j+1]
                - 0.5 * dy * (pStar[i+1][j+1] - pStar[i][j+1]);
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
        for (size_t j = 0; j < ny2; j++) uStar[i+1][j+1] = f[i];
    }

    std::ofstream f;
    f.open("test.csv");

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




    /*
    for (int row = 1; row < ny-1; row++) {
        size_t n = nx - 2;
        double a[n];
        double b[n];
        double c[n];
        double d[n];
        double f[n];

        for (size_t i = 1; i < nx - 1; i++) {
            b[i-1] = A[i-1][row-1];
            a[i-1] = C[i-1][row-1];
            c[i-1] = B[i-1][row-1];
            d[i-1] = F[i-1][row-1] - E[i-1][row-1] * u[i][row+1] - D[i-1][row-1] * u[i][row-1];

            if (i == nx -2) std::cout << u[i][row+1] << "\tCool!\n\n";
        }

        a[0] = 0;
        c[n-1] = 0;
        d[0] -= C[0][row-1] * u[0][row];
        d[n-1] -= B[n-1][row-1] * u[nx-1][row];

        ThomasAlg(&a[0], &b[0], &c[0], &d[0], &f[0], n);
        for (size_t i = 1; i < nx-1; i++) {
            std::cout << f[i-1] << "\t";
            uStar[i][row] = f[i-1];
        }
    }
    */
}
void Ueqn::ThomasAlg(double *a, double *b, double *c, double *d, double *f, size_t n) {
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


void Ueqn::WriteInfo() {
    std::ofstream f;
    f.open("output");

    f << "U Field\n";
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            f << u[i][j] << ",";
        }
        f << "\n";
    }

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
}

void Ueqn::PrintUStar() {
    for (size_t j = 0; j < ny; j++) {
        for (size_t i = 0; i < nx; i++) {
            std::cout << uStar[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}
