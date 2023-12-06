#pragma once
#include "../Geometry/Geometry.h"
#include <type_traits>

class StaggeredCell {
    public:
        StaggeredCell();
        StaggeredCell(Point cellCenter, double dx, double dy, size_t i, size_t j);
        ~StaggeredCell();

        const Point& GetCellCenter(void) const;
        const double& GetDX(void) const;
        const double& GetDY(void) const;
        const size_t& Geti(void) const;
        const size_t& Getj(void) const;

        void SetCellCenter(Point);
        void SetDX(double&);
        void SetDY(double&);

        template <typename T>
            double LambdaE(T **cArray, size_t nx, size_t ny) const {
                CheckIfDerived(cArray);
                CheckIfBoundary(nx, ny);
                double xij = ctr.GetX();
                double xe = xij + dx/2;
                double xip1j = cArray[i+1][j].ctr.GetX();
                return (xe - xij) / (xip1j - xij);
            }

        template <typename T>
            double LambdaW(T **cArray, size_t nx, size_t ny) const {
                CheckIfDerived(cArray);
                CheckIfBoundary(nx, ny);
                double xij = ctr.GetX();
                double xw = xij - dx/2;
                double xim1j = cArray[i-1][j].ctr.GetX();
                return (xw - xij) / (xim1j - xij);
            }

        template <typename T>
            double LambdaN(T **cArray, size_t nx, size_t ny) const {
                CheckIfDerived(cArray);
                CheckIfBoundary(nx, ny);
                double yij = ctr.GetY();
                double yn = yij + dy/2;
                double yijp1 = cArray[i][j+1].ctr.GetY();
                return (yn - yij) / (yijp1 - yij);
            }

        template <typename T>
            double LambdaS(T **cArray, size_t nx, size_t ny) const {
                CheckIfDerived(cArray);
                CheckIfBoundary(nx, ny);
                double yij = ctr.GetY();
                double ys = yij - dy/2;
                double yijm1 = cArray[i][j-1].ctr.GetY();
                return (ys - yij) / (yijm1 - yij);
            }

        template <typename T>
            double DelE(T **cArray, size_t nx, size_t ny) const {
                CheckIfDerived(cArray);
                CheckIfBoundary(nx, ny);
                double xij = ctr.GetX();
                double xip1j = cArray[i+1][j].ctr.GetX();
                return 1 / (xip1j - xij);
            }

        template <typename T>
            double DelW(T **cArray, size_t nx, size_t ny) const {
                CheckIfDerived(cArray);
                CheckIfBoundary(nx, ny);
                double xij = ctr.GetX();
                double xim1j = cArray[i-1][j].ctr.GetX();
                return 1 / (xim1j - xij);
            }

        template <typename T>
            double DelN(T **cArray, size_t nx, size_t ny) const {
                CheckIfDerived(cArray);
                CheckIfBoundary(nx, ny);
                double yij = ctr.GetY();
                double yijp1 = cArray[i][j+1].ctr.GetY();
                return 1 / (yijp1 - yij);
            }

        template <typename T>
            double DelS(T **cArray, size_t nx, size_t ny) const {
                CheckIfDerived(cArray);
                CheckIfBoundary(nx, ny);
                double yij = ctr.GetY();
                double yijm1 = cArray[i][j-1].ctr.GetY();
                return 1 / (yijm1 - yij);
            }

    private:
        Point ctr;
        double dx;
        double dy;
        size_t i;
        size_t j;

        void CheckIfBoundary(size_t nx, size_t ny) const {
            if (i == 0 || j == 0 || i == nx-1 || j == ny-1) 
                throw std::runtime_error("Attempted to interpolate at boundary cell.");
        }

        template <typename T>
        void CheckIfDerived(T **array) const {
            static_assert(std::is_base_of<StaggeredCell, T>::value, "T must be derived from 'StaggeredCell'");
        }

};
