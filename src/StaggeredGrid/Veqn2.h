#pragma once

#include "StaggeredGrid.h"

class Ueqn2;
class Peqn2;

class Veqn2 {
    public:
        Veqn2();
        Veqn2(StaggeredGrid *grid, double rho, double mu, 
                double timestep, double leftBC, 
                double rightBC, double topBC, double botBC);
        ~Veqn2();

        void InitField(double v0 = 0);

        void UpdateCoefficients(Ueqn2 &uEqn, Peqn2 &pEqn);

        void EnforceBC(void);

        void Print(void);

        double** GetV(void);

        void PrintVStar(void);

        double** GetVStar(void);

        void Solve(Peqn2 &pEqn);

        void ThomasAlg(double *a, double *b, double *c, double *d, double *f, size_t n);

        void Overwrite(void);

        friend class Peqn2;

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
        double **v;
        double **vStar;
        double **A;
        double **B;
        double **C;
        double **D;
        double **E;
        double **F;

        void AllocateMemory(void);
};
