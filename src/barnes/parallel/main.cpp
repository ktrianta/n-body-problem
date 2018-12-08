#include <mpi.h>
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

using namespace std;

int main(int argc, char** argv) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    time_point_t prog_start;
    time_point_t prog_end;
    time_point_t io_start;
    time_point_t io_end;
    time_point_t comp_start;
    time_point_t comp_end;
    time_point_t tree_start;
    time_point_t tree_end;
    time_point_t comm_start;
    time_point_t comm_end;
    double  prog_time = 0;
    double io_time = 0;
    double comp_time = 0;
    double tree_time = 0;
    double comm_time = 0;

    prog_start = std::chrono::high_resolution_clock::now();

    sim::Parameters params;
    readArgs(argc, argv, params);

    const size_t N = params.n;
    const sim::data_type theta = params.theta;
    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);

    io_start = std::chrono::high_resolution_clock::now();
    if (rank == 0) {
        if (!params.in_filename.empty()) {
            if (readDataFromFile(params.in_filename, N, m, r, u) == -1) {
                std::cerr << "File " << params.in_filename << " not found!" << std::endl;
                return -1;
            }
            params.out_filename = params.in_filename;
        } else {
              initializePositionOnSphere(N, r);
        }
    }
    std::ofstream out_file;
    if (rank == 0) {
        openFileToWrite(out_file, params.out_filename, params.out_dirname);
        writeDataToFile(N, r, u, out_file);
    }
    io_end = std::chrono::high_resolution_clock::now();
    io_time += std::chrono::duration< double >(io_end - io_start).count();

    // SEND the position vector r from Process 0 to all processes.
    comm_start = std::chrono::high_resolution_clock::now();
    MPI_Bcast(&r[0][0], N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&u[0][0], N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    comm_end = std::chrono::high_resolution_clock::now();
    comm_time += std::chrono::duration< double >(comm_end - comm_start).count();



    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(N, r, xc, yc, zc, w2, h2, t2);

    tree_start = std::chrono::high_resolution_clock::now();
    Octree tree = Octree(r, m, N, xc, yc, zc, w2, h2, t2);
    tree_end = std::chrono::high_resolution_clock::now();
    tree_time += std::chrono::duration< double >(tree_end - tree_start).count();

    const size_t Ntimesteps = params.s;
    const sim::data_type dt = params.dt;

    // Local Declarations
    size_t *local_N = new size_t[size];
    int *local_Nx3 = new int[size];
    size_t local_N_int = N/size;
    size_t rem = N - local_N_int * size;
    size_t counter = 0;

    for (size_t i = 0; i < size; i++) {
        local_N[i] = local_N_int;
        if (counter < rem) {
            local_N[i] += 1;
            counter ++;
        }
        local_Nx3[i] = local_N[i]*3;
    }

    sim::data_type (*r_local)[3] = new sim::data_type[local_N[rank]][3];

    size_t *offset = new size_t[size]; offset[0] = 0;
    int *offset_x3 = new int[size]; offset_x3[0] = 0;

    for (size_t i = 1; i < size; i++) {
        offset[i] = offset[i-1] + local_N[i-1];
        offset_x3[i] = offset_x3[i-1] + local_Nx3[i-1];
    }

    for (size_t i = 0, j = offset[rank], end = local_N[rank]; i < end; i++, j++) {
        r_local[i][0] = r[j][0];
        r_local[i][1] = r[j][1];
        r_local[i][2] = r[j][2];
    }

    comp_start = std::chrono::high_resolution_clock::now();
    for (size_t j = offset[rank], end = local_N[rank] + offset[rank]; j < end; j++) {
        tree.computeAcceleration(j, r, a, sim::g, theta);
    }
    comp_end = std::chrono::high_resolution_clock::now();
    comp_time += std::chrono::duration< double >(comp_end - comp_start).count();


    for (size_t t = 0; t < Ntimesteps; t++) {
        for (size_t j = 0, idx = offset[rank], end = local_N[rank]; j < end; j++, idx++) {
            u[idx][0] += 0.5 * a[idx][0] * dt;
            u[idx][1] += 0.5 * a[idx][1] * dt;
            u[idx][2] += 0.5 * a[idx][2] * dt;
            r_local[j][0] += u[idx][0] * dt;
            r_local[j][1] += u[idx][1] * dt;
            r_local[j][2] += u[idx][2] * dt;
            a[idx][0] = 0;
            a[idx][1] = 0;
            a[idx][2] = 0;
        }

        comm_start = std::chrono::high_resolution_clock::now();
        MPI_Allgatherv(&(r_local[0][0]), local_N[rank]*3, MPI_DOUBLE, &(r[0][0]),
            local_Nx3, offset_x3, MPI_DOUBLE, MPI_COMM_WORLD);
        comm_end = std::chrono::high_resolution_clock::now();
        comm_time += std::chrono::duration< double >(comm_end - comm_start).count();

        for (size_t idx = offset[rank], end = offset[rank] + local_N[rank]; idx < end; idx++) {
            comp_start = std::chrono::high_resolution_clock::now();
            tree.computeAcceleration(idx, r, a, sim::g, theta);
            comp_end = std::chrono::high_resolution_clock::now();
            comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

            u[idx][0] += 0.5 * a[idx][0] * dt;
            u[idx][1] += 0.5 * a[idx][1] * dt;
            u[idx][2] += 0.5 * a[idx][2] * dt;

        }

        boxComputation(N, r, xc, yc, zc, w2, h2, t2);
        tree_start = std::chrono::high_resolution_clock::now();
        Octree tree = Octree(r, m, N, xc, yc, zc, w2, h2, t2);
        tree_end = std::chrono::high_resolution_clock::now();
        tree_time += std::chrono::duration< double >(tree_end - tree_start).count();

        io_start = std::chrono::high_resolution_clock::now();
        if (t % 200 == 0 && rank == 0) {
            writeDataToFile(N, r, u, out_file);
        }
        io_end = std::chrono::high_resolution_clock::now();
        io_time += std::chrono::duration< double >(io_end - io_start).count();
    }

    prog_end = std::chrono::high_resolution_clock::now();
    prog_time += std::chrono::duration< double >(prog_end - prog_start).count();
    double plotData_comp;
    double plotData_comm;
    double plotData_io;
    double plotData_prog;
    double plotData_tree;
    MPI_Reduce(&comp_time, &plotData_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time, &plotData_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&io_time, &plotData_io, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&prog_time, &plotData_prog, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tree_time, &plotData_tree, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        delete[] r;
        FILE *plotFile;
        plotFile = fopen("plotData.txt", "a");
        fprintf(plotFile, "%lf,  %lf, %lf, %lf, %lf \n", plotData_prog, plotData_comp, plotData_io, plotData_tree, plotData_comm);
        fclose(plotFile);
    }

    delete[] m;
    delete[] r;
    delete[] u;
    delete[] a;
    delete[] local_N;
    delete[] local_Nx3;
    delete[] offset;
    delete[] offset_x3;
    delete[] r_local;

    MPI_Finalize();
    return 0;
}
