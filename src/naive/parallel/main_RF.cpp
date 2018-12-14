#include <mpi.h>
#include <cmath>
#include <fstream>
#include <chrono>
#include <iostream>
#include <papi.h>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"

using time_point_t = std::chrono::high_resolution_clock::time_point;


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
    
    // *** PAPI *** //
    int EventSet = PAPI_NULL;
    long_long values[1] = {(long_long) 0};

    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        printf("PAPI library init error!\n");
        exit(1);
    }

    if ( PAPI_create_eventset(&EventSet) != PAPI_OK){
        printf("PAPI library create error!\n");
        exit(1);
    }

    if ( PAPI_add_event(EventSet,PAPI_TOT_CYC ) != PAPI_OK){
        printf("PAPI add PAPI_L1_DCM error!\n");
        exit(1);
    }

    long long int papi_tot = 0;
    // *** PAPI *** //
    

    sim::Parameters params;
    readArgs(argc, argv, params);

    PAPI_start(EventSet);
    const size_t N = params.n;
    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);

    // PROCESS 0 initialize position vector r.
    std::ofstream out_file;
    if (rank == 0) {
        if (readDataFromFile(params.in_filename, N, m, r, u) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            delete[] m;
            delete[] r;
            delete[] u;
            delete[] a;
            return -1;
        }
        params.out_filename = params.in_filename;
        openFileToWrite(out_file, params.out_filename, params.out_dirname);
        writeDataToFile(N, r, out_file);
    }
    PAPI_stop(EventSet, values);
    papi_tot += values[0]; 


    // SEND the position vector r from Process 0 to all processes.
    MPI_Bcast(&r[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&u[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&m[0],N, MPI_DOUBLE,0, MPI_COMM_WORLD);


    PAPI_start(EventSet);
    const sim::data_type dt = params.dt;
    const size_t timesteps = params.s;

    size_t local_N[size];
    int local_Nx3[size];
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

    size_t offset[size]; offset[0] = 0;
    int offset_x3[size]; offset_x3[0] = 0;

    for (size_t i = 1; i < size; i++) {
        offset[i] = offset[i-1] + local_N[i-1];
        offset_x3[i] = offset_x3[i-1] + local_Nx3[i-1];
    }

    for (size_t i = 0, j = offset[rank], end = local_N[rank]; i < end; i++, j++) {
        r_local[i][0] = r[j][0];
        r_local[i][1] = r[j][1];
        r_local[i][2] = r[j][2];
    }

    computeAcceleration(N, r, a, m, local_N[rank], offset[rank]);

    //Start benchmark
    for (size_t t = 0; t < timesteps; t++) {
        for (size_t j = 0, idx = offset[rank]; j < local_N[rank]; j++, idx++) {
            u[idx][0] += 0.5 * a[idx][0] * dt;
            u[idx][1] += 0.5 * a[idx][1] * dt;
            u[idx][2] += 0.5 * a[idx][2] * dt;
            r_local[j][0] += u[idx][0] * dt;
            r_local[j][1] += u[idx][1] * dt;
            r_local[j][2] += u[idx][2] * dt;
        }
        PAPI_stop(EventSet, values);
        papi_tot += values[0]; 

        MPI_Allgatherv(&(r_local[0][0]), local_N[rank]*3, MPI_DOUBLE, &(r[0][0]),
            local_Nx3, offset_x3, MPI_DOUBLE, MPI_COMM_WORLD);

        PAPI_start(EventSet);
        computeAcceleration(N, r, a, m, local_N[rank], offset[rank]);

        for (size_t idx = offset[rank], end = offset[rank] + local_N[rank]; idx < end; idx++) {
            u[idx][0] += 0.5 * a[idx][0] * dt;
            u[idx][1] += 0.5 * a[idx][1] * dt;
            u[idx][2] += 0.5 * a[idx][2] * dt;
	}
        PAPI_stop(EventSet, values);
        papi_tot += values[0]; 

        if (rank == 0) {
            if (t % 200 == 0){
                writeDataToFile(N, r, out_file);
            }   
        PAPI_start(EventSet);
        }
    }
    PAPI_stop(EventSet, values);
    papi_tot += values[0]; 

    printf("%lld\n", papi_tot);

    delete[] m;
    delete[] r;
    delete[] u;
    delete[] a;

    MPI_Finalize();
    return 0;
}
