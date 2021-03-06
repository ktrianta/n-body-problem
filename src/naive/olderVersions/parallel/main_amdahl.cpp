#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"

using namespace std;



void computeAcceleration(const size_t N, sim::data_type (*r)[3], sim::data_type (*a)[3], sim::data_type *m, const int local_N, const int offset) {
    std::fill(&a[offset][0], &a[offset][0] + local_N*3, 0);

    for (size_t i = offset; i < offset + local_N; i++) {
        sim::data_type a_i0 = 0;  // accumulate accelaration values for particle i and
        sim::data_type a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        sim::data_type a_i2 = 0;

        for (size_t j = 0; j < N; j++) {
            if (i == j) continue;

            sim::data_type rji[3];
            rji[0] = r[j][0] - r[i][0];
            rji[1] = r[j][1] - r[i][1];
            rji[2] = r[j][2] - r[i][2];
            sim::data_type r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];
            sim::data_type denom = (r2+sim::e2) * sqrt(r2+sim::e2);
            sim::data_type a_i = - sim::g * m[j] / denom;
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

    // *** MPI *** // 
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // *** MPI *** // 

    sim::Parameters params;
    readArgs(argc, argv, params);
    for ( int l = 0; l < 100; l++ ) {
        printf("%d\n",l);
        double t_tot_1, t_tot_2, t_tot, t_par_1, t_par_2, t_par, t_comm = 0;
        t_par = 0;
        t_tot_1 = MPI_Wtime();

        const size_t N = params.n;
        sim::data_type *m = new sim::data_type[N];
        sim::data_type (*r)[3] = new sim::data_type[N][3];
        sim::data_type (*u)[3] = new sim::data_type[N][3];
        sim::data_type (*a)[3] = new sim::data_type[N][3];
        std::fill(m, m + N, 1.0 / N);
        std::fill(&u[0][0], &u[0][0] + N * 3, 0);
        std::fill(&a[0][0], &a[0][0] + N * 3, 0);

        t_tot_2 = MPI_Wtime();
        t_tot = t_tot_2 - t_tot_1;

        // PROCESS 0 initialize position vector r.
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
        MPI_Bcast(&r[0][0], N * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&u[0][0], N * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&m[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        std::ofstream out_file;
        if (rank == 0) {
            openFileToWrite(out_file, params.out_filename, params.out_dirname);
            writeDataToFile(params.n, r, u, out_file);
        }

        t_tot_1 = MPI_Wtime();
        const sim::data_type dt = params.dt;
        const size_t timesteps = params.t / dt + 1;

        size_t local_N[size];
        int local_Nx3[size];
        size_t local_N_int = N / size;
        size_t rem = N - local_N_int * size;
        size_t counter = 0;

        for (size_t i = 0; i < size; i++) {
            local_N[i] = local_N_int;
            if (counter < rem) {
                local_N[i] += 1;
                counter++;
            }
            local_Nx3[i] = local_N[i] * 3;
        }

        size_t offset[size];
        offset[0] = 0;
        int offset_x3[size];
        offset_x3[0] = 0;

        for (size_t i = 1; i < size; i++) {
            offset[i] = offset[i - 1] + local_N[i - 1];
            offset_x3[i] = offset_x3[i - 1] + local_Nx3[i - 1];
        }

        t_par_1 = MPI_Wtime();
        sim::data_type (*r_local)[3] = new sim::data_type[local_N[rank]][3];

        for (size_t i = 0, j = offset[rank], end = local_N[rank]; i < end; i++, j++) {
            r_local[i][0] = r[j][0];
            r_local[i][1] = r[j][1];
            r_local[i][2] = r[j][2];
        }

        computeAcceleration(N, r, a, m, local_N[rank], offset[rank]);

        for (size_t t = 0; t < timesteps; t++) {
            for (size_t j = 0, idx = offset[rank]; j < local_N[rank]; j++, idx++) {
                u[idx][0] += 0.5 * a[idx][0] * dt;
                u[idx][1] += 0.5 * a[idx][1] * dt;
                u[idx][2] += 0.5 * a[idx][2] * dt;
                r_local[j][0] += u[idx][0] * dt;
                r_local[j][1] += u[idx][1] * dt;
                r_local[j][2] += u[idx][2] * dt;
            }
            t_tot_2 = MPI_Wtime();
            t_par_2 = t_tot_2;
            t_tot += t_tot_2 - t_tot_1;
            t_par += t_par_2 - t_par_1;
            MPI_Allgatherv(&(r_local[0][0]), local_N[rank] * 3, MPI_DOUBLE, &(r[0][0]),
                           local_Nx3, offset_x3, MPI_DOUBLE, MPI_COMM_WORLD);
            t_tot_1 = MPI_Wtime();
            t_par_1 = t_tot_1;
            t_comm += t_tot_1 - t_par_2;
            computeAcceleration(N, r, a, m, local_N[rank], offset[rank]);

            for (size_t idx = offset[rank], end = offset[rank] + local_N[rank]; idx < end; idx++) {
                u[idx][0] += 0.5 * a[idx][0] * dt;
                u[idx][1] += 0.5 * a[idx][1] * dt;
                u[idx][2] += 0.5 * a[idx][2] * dt;
            }

            t_tot_2 = MPI_Wtime();
            t_par_2 = t_tot_2;
            t_tot += t_tot_2 - t_tot_1;
            t_par += t_par_2 - t_par_1;
            if (rank == 0) {
                if (t % 200 == 0) {
                    writeDataToFile(N, r, u, out_file);
                }
            }
            t_tot_1 = MPI_Wtime();
            t_par_1 = t_tot_1;
        }
        t_par_2 = MPI_Wtime();
        t_tot_2 = t_par_2;
        t_par += t_par_2 - t_par_1;
        t_tot += t_tot_2 - t_tot_1;
        double plot_data_tot;
        double plot_data_par;
        double plot_data_comm;
        MPI_Reduce(&t_par, &plot_data_par, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_tot, &plot_data_tot, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_comm, &plot_data_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            FILE *plotFile;
            plotFile = fopen("plotData.txt", "a");
            fprintf(plotFile, "%d %lf %lf %lf\n", size, plot_data_tot, plot_data_par, plot_data_comm);
        }
    }
        MPI_Finalize();
        return 0;

}
