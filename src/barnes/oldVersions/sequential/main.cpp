#include <math.h>
#include <fstream>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "octree.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"
#include <chrono>


using time_point_t = std::chrono::high_resolution_clock::time_point;

int main(int argc, char** argv) {

    time_point_t prog_start;
    time_point_t prog_end;
    time_point_t io_start;
    time_point_t io_end;
    time_point_t comp_start;
    time_point_t comp_end;
    time_point_t tree_start;
    time_point_t tree_end;
    double  prog_time = 0;
    double io_time = 0;
    double comp_time = 0;
    double tree_time = 0;

    prog_start = std::chrono::high_resolution_clock::now();

    sim::Parameters params;
    readArgs(argc, argv, params);

    const size_t N = params.n;
    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);
 
    io_start = std::chrono::high_resolution_clock::now();
    if (readDataFromFile(params.in_filename, N, m, r, u) == -1) {
        std::cerr << "File " << params.in_filename << " not found!" << std::endl;
        delete[] m;
        delete[] r;
        delete[] u;
        delete[] a;
        return -1;
    }
    params.out_filename = params.in_filename;
    std::ofstream out_file;
    openFileToWrite(out_file, params.out_filename, params.out_dirname);
    writeDataToFile(N, r, u, out_file);

    io_end = std::chrono::high_resolution_clock::now();
    io_time += std::chrono::duration< double >(io_end - io_start).count();

    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(N, r, xc, yc, zc, w2, h2, t2);

    Octree::r = r;
    Octree::a = a;
    Octree::m = m;
    Octree::g = sim::g;
    Octree::theta = params.theta;
    tree_start = std::chrono::high_resolution_clock::now();
    Octree tree = Octree(N, xc, yc, zc, w2, h2, t2);
    tree_end = std::chrono::high_resolution_clock::now();
    tree_time += std::chrono::duration< double >(tree_end - tree_start).count();

// Computation of Initial Energy    
    double initialKEnergy = 0;
    double initialPEnergy = 0;
    double initialEnergy = 0;
    for (int i = 0; i < N; i++){
        initialKEnergy += m[i] * (u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])/2.;
        for (int j = 0; j < i; j++){
            double denominator = sqrt((r[j][0]-r[i][0])*(r[j][0]-r[i][0]) + (r[j][1]-r[i][1])*(r[j][1]-r[i][1]) +
                                      (r[j][2]-r[i][2])*(r[j][2]-r[i][2]));
            initialPEnergy -= sim::g*m[i]*m[j]/denominator;
            }
        }
        initialEnergy = initialKEnergy + initialPEnergy;

    comp_start = std::chrono::high_resolution_clock::now();
    for (size_t j = 0; j < N; j++) {
        tree.computeAcceleration(j);
    }
    comp_end = std::chrono::high_resolution_clock::now();
    comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

    const size_t Ntimesteps = params.s;
    const sim::data_type dt = params.dt;

    for (size_t t = 0; t < Ntimesteps; t++) {
        for (size_t j = 0; j < N; j++) {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
            r[j][0] += u[j][0] * dt;
            r[j][1] += u[j][1] * dt;
            r[j][2] += u[j][2] * dt;

        }

        boxComputation(N, r, xc, yc, zc, w2, h2, t2);
        tree_start = std::chrono::high_resolution_clock::now();
        Octree tree = Octree(N, xc, yc, zc, w2, h2, t2);
        tree_end = std::chrono::high_resolution_clock::now();
        tree_time += std::chrono::duration< double >(tree_end - tree_start).count();
        std::fill(&a[0][0], &a[0][0] + N*3, 0);

        for (size_t j = 0; j < N; j++) {
            comp_start = std::chrono::high_resolution_clock::now();
            tree.computeAcceleration(j);
            comp_end = std::chrono::high_resolution_clock::now();
            comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
        }

        if (t % 200 == 0) {
            writeDataToFile(N, r, u, out_file);
        }
    }

    double energy =0;
    double kineticEnergy = 0;
    double potentialEnergy = 0;
    for (int i = 0; i < N; i++){
        kineticEnergy += m[i] * (u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])/2.;
        for (int j = 0; j < i; j++){
            double denominator = sqrt((r[j][0]-r[i][0])*(r[j][0]-r[i][0]) + (r[j][1]-r[i][1])*(r[j][1]-r[i][1]) +
                                      (r[j][2]-r[i][2])*(r[j][2]-r[i][2]));
            potentialEnergy -= sim::g*m[i]*m[j]/denominator;
            }
        }
        energy = kineticEnergy + potentialEnergy;
    std::cout << "initial energy is = " << initialEnergy << "Error in total energy at the end of simulation = " << (energy - initialEnergy)/initialEnergy*100 << "%" <<  std::endl; 

    prog_end = std::chrono::high_resolution_clock::now();
    prog_time += std::chrono::duration< double >(prog_end - prog_start).count();
    FILE *plotFile;
    plotFile = fopen("plotData.txt", "a");
    fprintf(plotFile, "%lf,  %lf, %lf, %lf \n", prog_time, comp_time, io_time, tree_time);
    fclose(plotFile);

    delete[] m;
    delete[] r;
    delete[] u;
    delete[] a;

    return 0;
}
