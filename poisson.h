#ifndef _POISSON_H
#define _POISSON_H
#include <omp.h>
#include <cmath>
#include <stdio.h>
template<typename dtype>
class PoissonSolver {
  public:
    PoissonSolver() = delete;
    PoissonSolver(const int, const int, const double, const double, const int);
    ~PoissonSolver();
    virtual void operator()(dtype**, dtype**);
  private:
    const dtype tol = 1e-10;
    const int iter_max;
    const int nx;
    const int ny;
    const double dx;
    const double dy;
    dtype** buffer;
};

template<typename dtype>
inline PoissonSolver<dtype>::PoissonSolver(const int nx, const int ny, const double lx, const double ly, const int iter_max) :
    iter_max(iter_max), nx(nx), ny(ny), dx(lx/nx), dy(lx/ny)
{
    this->buffer = new dtype* [this->nx];
    for (int i=0; i<nx; i++) this->buffer[i] = new dtype[this->ny];
}



template<typename dtype>
inline PoissonSolver<dtype>::~PoissonSolver() {
    for (int i=0; i<nx; i++) delete [] this->buffer[i];
    delete [] this->buffer;
}


template<typename dtype>
inline void PoissonSolver<dtype>::operator() (dtype** sol, dtype** src) {
    const dtype dx2 = dx * dx;
    const dtype dy2 = dy * dy;
    const dtype err_norm = 1.0 / (dtype(nx) * dtype(ny));
    dtype** sol_past = sol;
    dtype** sol_new  = buffer;
    dtype err = tol * 10.0;
    int counter = 0;

    while (err*err > tol*tol && counter < iter_max) {
        err = 0.;
#pragma omp parallel for
        //default(none) shared(sol_new,err,err_norm, src, dx2, sol_past) //schedule(static,omp_get_num_threads())  
            for (int i = 0; i < nx; i++){
//            printf("from thread %d, i = %d \n",omp_get_thread_num(),i);
            const int ip = (i + 1) % nx;
            const int im = (i + nx - 1) % nx;
            for (int j = 0; j < ny; j++){
                const int jp = (j + 1) % ny;
                const int jm = (j + ny - 1) % ny;
                sol_new[i][j] = -0.25 * (dx2 * src[i][j] - (sol_past[im][j] + sol_past[ip][j] + sol_past[i][jm] + sol_past[i][jp])) - src[0][0]; //last src[0][0] added to ensure that solution is 0 at src[0][0] boundary condition
#pragma atomic
                err += fabs(sol_new[i][j] -sol_past[i][j]) * err_norm;
            }
        }
        counter++;
        if (counter % 2) {
            sol_past = buffer;
            sol_new  = sol;
        } else {
            sol_past = sol;
            sol_new  = buffer;
        }
    }
    if ((counter % 2) == 1) {
        // copy buffer to sol
        for (int i = 0; i < nx; i++){
            for (int j = 0; j < ny; j++){
                sol[i][j] = buffer[i][j];
            }
        }
    }
printf("Completed the object function with %d iterations! :) \n", counter);

}
#endif
