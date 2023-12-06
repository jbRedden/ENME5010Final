#pragma once

#include "StaggeredGrid.h"
#include "Veqn.h"
#include "Peqn.h"

class Ueqn {
    public:
        Ueqn();
        Ueqn(StaggeredGrid *grid, double rho, double mu, 
                double timestep, double leftBC, 
                double rightBC, double topBC, double botBC);
        ~Ueqn();

        double** GetU(void);
        double** GetUStar(void);

        void InitField(double u0 = 0);

        /*Updates A, B, C, D, E, F coefficients. Only call once per timestep*/
        void UpdateCoefficients(Veqn &vEqn, Peqn &pEqn);

        /*Copies uStar values to u*/
        void UpdateForNextTimeStep(void);

        void EnforceBC(void);

        void Print(void);
        void PrintUStar(void);

        void Solve(size_t j);
        void Solv2(size_t j, Peqn &pEqn);
        void WriteInfo(void);

        void Overwrite(void);

        friend class Peqn;

    private:
        StaggeredGrid *grid;
        size_t nx;
        size_t ny;
        double rho;
        double mu;
        double dt;
        double leftBC;
        double rightBC;
        double topBC;
        double botBC;
        double **u;
        double **uStar;
        double **A;
        double **B;
        double **C;
        double **D;
        double **E;
        double **F;

        void ThomasAlg(double *a, double *b, double *c, double *d, double *f, size_t n);

        void AllocateMemory(void);

};
