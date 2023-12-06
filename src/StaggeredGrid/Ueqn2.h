#pragma once

#include "StaggeredGrid.h"
#include "Veqn2.h"
#include "Peqn2.h"

class Ueqn2 {
    public:
        Ueqn2();
        Ueqn2(StaggeredGrid *grid, double rho, double mu, 
                double timestep, double leftBC, 
                double rightBC, double topBC, double botBC);
        ~Ueqn2();

        double** GetU(void);
        double** GetUStar(void);

        void InitField(double u0 = 0);

        /*Updates A, B, C, D, E, F coefficients. Only call once per timestep*/
        void UpdateCoefficients(Veqn2 &vEqn, Peqn2 &pEqn);

        /*Copies uStar values to u*/
        void UpdateForNextTimeStep(void);

        void EnforceBC(void);

        void Print(void);
        void PrintUStar(void);

        void Solve(Peqn2 &pEqn);
        void WriteInfo(void);

        void Overwrite(void);

        void Testing(void);

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
