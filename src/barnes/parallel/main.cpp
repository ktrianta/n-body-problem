#include <mpi.h>
#include <fstream>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "octree.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"

using namespace std;

int main(int argc, char** argv) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

    // SEND the position vector r from Process 0 to all processes.
    MPI_Bcast(&r[0][0], N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&u[0][0], N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(N, r, xc, yc, zc, w2, h2, t2);

    std::ofstream out_file;
    if (rank == 0) {
        openFileToWrite(out_file, params.out_filename, params.out_dirname);
        writeDataToFile(N, r, u, out_file);
    }

    Octree tree = Octree(r, m, N, xc, yc, zc, w2, h2, t2);

    const size_t Ntimesteps = params.t / params.dt + 1;
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

    for (size_t j = offset[rank], end = local_N[rank] + offset[rank]; j < end; j++) {
        tree.computeAcceleration(j, r, a, sim::g, theta);
    }

    double t_start = 0;
    double t_end = 0;
    double time = 0;

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

        MPI_Allgatherv(&(r_local[0][0]), local_N[rank]*3, MPI_DOUBLE, &(r[0][0]),
            local_Nx3, offset_x3, MPI_DOUBLE, MPI_COMM_WORLD);

        for (size_t idx = offset[rank], end = offset[rank] + local_N[rank]; idx < end; idx++) {
            t_start = MPI_Wtime();
            tree.computeAcceleration(idx, r, a, sim::g, theta);
            t_end = MPI_Wtime();
            time += t_end - t_start;

            u[idx][0] += 0.5 * a[idx][0] * dt;
            u[idx][1] += 0.5 * a[idx][1] * dt;
            u[idx][2] += 0.5 * a[idx][2] * dt;

        }

        boxComputation(N, r, xc, yc, zc, w2, h2, t2);
        Octree tree = Octree(r, m, N, xc, yc, zc, w2, h2, t2);

        if (t % 200 == 0 && rank == 0) {
            writeDataToFile(N, r, u, out_file);
        }
    }

    printf("time= %f \n", time);

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
