#include <cstdio>
#include <cassert>
#include <string>
#include "poisson.h"
#include <cmath>
#include <omp.h>
#include <stdio.h>
#define PI 3.14159265358979323846

int getInt(int argc, char* argv[], std::string flag, int default_val){
    for (int i=0; i<argc; i++) {
        if (std::string(argv[i]) == flag) {
            try{
                const int val = atoi(argv[i+1]);
                assert(val > 0);
                return val;
            } catch (...) {}
        }
    }
    return default_val;
}


void writeFieldToFile(const std::string filename, double** field, const int Nx, const int Ny) {
    FILE* file = fopen(filename.c_str(), "w");
    for (int i=0; i<Nx; i++) {
        for (int j=0; j<Ny; j++) {
            fprintf(file, "%23.16E\t", field[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}


double source_function(const double x_loc, const double y_loc,const double Lx, const double Ly){
    return sin(x_loc/Lx * 2 * PI) * cos(y_loc/Ly * 2 * PI);
}


void print_matrix(double** my_matrix, const int Nx, const int Ny){
    for (int j = 0; j < Ny;j++){
        for (int i = 0; i < Nx; i++){
            printf("%23.16e\t", my_matrix[i][j]);
        }
        printf("\n");
    }

}

int main(int argc, char* argv[]) {
//    printf("OpenMP running with %d threads \n", omp_get_max_threads());
// #pragma omp parallel
// {
// printf("Hello world from thread %d \n",omp_get_thread_num());
// }
    const int Nx = getInt(argc, argv, "-N", 16);
    const int Ny = Nx;
    const int iter_max = getInt(argc, argv, "-I", 999);
    const double Lx = 2.0 * PI;
    const double Ly = 2.0 * PI;
    printf("System size: Nx = Ny = %d\n", Nx);

    double** src     = new double* [Nx];
    double** sol_ref = new double* [Nx];
    double** sol_num = new double* [Nx];
    for (int i=0; i<Nx; i++) {
        src    [i] = new double [Ny];
        sol_num[i] = new double [Ny];
        sol_ref[i] = new double [Ny];
    }

    // Set test case
    const double Ax = 2.0 * PI / Lx;
    const double Ay = 2.0 * PI / Ly;
    const double Axy=-1.0/(Ax*Ax+Ay*Ay);
    for (int i=0; i<Nx; i++) {
        const double x = Lx / Nx * i;
        for (int j=0; j<Ny; j++) {
            const double y = Ly / Ny * j;
            src    [i][j] = source_function(x, y, Lx, Ly);
            sol_ref[i][j] = Axy * source_function(x, y, Lx, Ly);
        }
    }

    PoissonSolver<double> poisson(Nx, Ny, Lx, Ly, iter_max);
    printf("Initial condition matrix: \n");
    //print_matrix(sol_num, Nx,Ny);
  
    const double tic = omp_get_wtime();
    poisson(sol_num, src);
    const double toc = omp_get_wtime();
    printf("Elapse = %.5E sec.\n", toc-tic);
    writeFieldToFile("./sol_num.dat", sol_num, Nx, Ny);
    writeFieldToFile("./sol_ref.dat", sol_ref, Nx, Ny);

    //print_matrix(src,Nx,Ny);
    // TODO: Calculate error

    // Free memory
    for (int i=0; i<Nx; i++) {
        delete [] src    [i];
        delete [] sol_ref[i];
        delete [] sol_num[i];
    }
    delete [] src;
    delete [] sol_ref;
    delete [] sol_num;

    return 0;
}
