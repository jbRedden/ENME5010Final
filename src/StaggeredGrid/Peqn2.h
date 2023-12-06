#pragma once

#include "StaggeredGrid.h"
#include "Veqn2.h"
#include "Ueqn2.h"

class Peqn2 {
    public:
        Peqn2();
        Peqn2(StaggeredGrid *grid, double timestep);
        ~Peqn2();

        void InitField(double p0 = 0);

        void UpdateCoefficients(Ueqn2 &uEqn, Veqn2 &vEqn);
        void Solve(Ueqn2 &uEqn, Veqn2 &vEqn);

        void Print(void);

        double** GetP(void);

        double** GetPStar(void);

        void EnforceBC(void);

        void Overwrite(void);

        double uCorrMax = 0;
        double vCorrMax = 0;
        double pCorrMax = 0;

        

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
