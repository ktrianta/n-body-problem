#include <math.h>
#include <fstream>
#include <chrono>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"

using time_point_t = std::chrono::high_resolution_clock::time_point;


void computeAcceleration(const size_t N, sim::data_type (*r)[3], sim::data_type (*a)[3], sim::data_type *m) {
    std::fill(&a[0][0], &a[0][0] + N*3, 0);

    for (size_t i = 0; i < N; i++) {
        sim::data_type a_i0 = 0;  // accumulate accelaration values for particle i and
        sim::data_type a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        sim::data_type a_i2 = 0;

        for (size_t j = i+1; j < N; j++) {
            sim::data_type rji[3];
            rji[0] = r[j][0] - r[i][0];
            rji[1] = r[j][1] - r[i][1];
            rji[2] = r[j][2] - r[i][2];
            sim::data_type r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];
            sim::data_type denom = (r2+sim::e2) * sqrt(r2+sim::e2);
            sim::data_type a_j = - sim::g * m[i] / denom;
            sim::data_type a_i = - sim::g * m[j] / denom;
            a[j][0] += a_j * rji[0];
            a[j][1] += a_j * rji[1];
            a[j][2] += a_j * rji[2];
            a_i0 -= a_i * rji[0];
            a_i1 -= a_i * rji[1];
            a_i2 -= a_i * rji[2];
        }
        a[i][0] += a_i0;  // a(i, 0) and a(i, 1) are accessed once here, avoiding
        a[i][1] += a_i1;  // repeated accesses in the inner loop of j
        a[i][2] += a_i2;
    }
}

int main(int argc, char** argv) {

    time_point_t prog_start;
    time_point_t prog_end;
    time_point_t io_start;
    time_point_t io_end;
    time_point_t comp_start;
    time_point_t comp_end;
    double  prog_time = 0;
    double io_time = 0;
    double comp_time = 0;

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
    writeDataToFile(params.n, r, u, out_file);

    io_end = std::chrono::high_resolution_clock::now();
    io_time += std::chrono::duration< double >(io_end - io_start).count();


    computeAcceleration(params.n, r, a, m);

    const size_t Ntimesteps = params.t / params.dt + 1;
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

      	comp_start = std::chrono::high_resolution_clock::now();
        computeAcceleration(N, r, a, m);
        comp_end = std::chrono::high_resolution_clock::now();
    	comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

        for (size_t j = 0; j < N; j++) {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
        }

        io_start = std::chrono::high_resolution_clock::now();
        if (t % 200 == 0) {
            writeDataToFile(N, r, u, out_file);
        }
        io_end = std::chrono::high_resolution_clock::now();
        io_time += std::chrono::duration< double >(io_end - io_start).count();
    }

    prog_end = std::chrono::high_resolution_clock::now();
    prog_time += std::chrono::duration< double >(prog_end - prog_start).count();
    FILE *plotFile;
    plotFile = fopen("plotData.txt", "a");
    fprintf(plotFile, "%lf,  %lf, %lf \n", prog_time, comp_time, io_time);
    fclose(plotFile);

    delete[] m;
    delete[] r;
    delete[] u;
    delete[] a;

    return 0;
}
