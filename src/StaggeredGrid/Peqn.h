#pragma once

#include "StaggeredGrid.h"
#include "Veqn.h"
#include "Ueqn.h"

class Peqn {
    public:
        Peqn();
        Peqn(StaggeredGrid *grid, double timestep);
        ~Peqn();

        void InitField(double p0 = 0);

        void UpdateCoefficients(Ueqn &uEqn, Veqn &vEqn);
        void Solve(Ueqn &uEqn, Veqn &vEqn);

        void Print(void);

        double** GetP(void);

        double** GetPStar(void);

        void EnforceBC(void);

        void Overwrite(void);

        

    private:
        StaggeredGrid *grid;
        size_t nx;
        size_t ny;
        double dt;
        double **p;
        double **pStar;
        double **A;
        double **B;
        double **C;
        double **D;
        double **E;
        double **F;

        void AllocateMemory(void);
        void ThomasAlg(double *a, double *b, double *c, double *d, double *f, size_t n);
        
};
