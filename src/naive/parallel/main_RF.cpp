#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>
#include <papi.h>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"

using namespace std;
using time_point_t = std::chrono::high_resolution_clock::time_point;

void computeAcceleration(const size_t N, sim::data_type (*r)[3], sim::data_type (*a)[3], sim::data_type* m, const int local_N, const int offset) {
    std::fill(&a[offset][0], &a[offset][0] + local_N*3, 0);

    for (int i = offset; i < offset + local_N; i++) {
        sim::data_type a_i0 = 0;  // accumulate accelaration values for particle i and
        sim::data_type a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        sim::data_type a_i2 = 0;

        for (int j = 0; j < N; j++) {
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
    MPI_Win win;
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
        printf("PAPI add event error!\n");
        exit(1);
    }

    long long int papi_tot = 0;
    // *** PAPI *** //

    PAPI_start(EventSet);

    sim::Parameters params;
    readArgs(argc, argv, params);

    const size_t N = params.n;
    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*r1)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&r1[0][0], &r1[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);
    
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
    }

    MPI_Bcast (&r[0][0], 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&u[0][0], 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&m[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    

    MPI_Win_create( &(r[0][0]), sizeof(double)*3*N, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);


    const sim::data_type dt = params.dt;
    const size_t timesteps = params.s;

    int local_N[size];
    int local_Nx2[size];
    int local_N_int = N/size;
    int rem = N - local_N_int * size;
    int counter = 0;

    for (int i = 0; i < size; i++) {
        local_N[i] = local_N_int;
        if (counter < rem) {
                local_N[i] += 1;
                counter ++;
        }
        local_Nx2[i] = local_N[i]*3;
    }

    int offset[size];
    offset[0] = 0;
    int offset_x2[size];
    offset_x2[0] = 0;
    

    for (int i=1; i < size; i++) {
        offset[i] = offset[i-1] + local_N[i-1];
        offset_x2[i] = offset_x2[i-1] + local_Nx2[i-1];
    }

    computeAcceleration(N, r, a, m, local_N[rank], offset[rank]);

    for (int t = 0; t < timesteps; t++) {
        for (int j = offset[rank]; j < offset[rank]+local_N[rank]; j++) {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
            r[j][0] += u[j][0] * dt;
            r[j][1] += u[j][1] * dt;
            r[j][2] += u[j][2] * dt;
        }
        
        MPI_Win_fence(MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE, win);
        for (int i = 0; i < size; i++) {
            MPI_Get(&(r1[offset[i]][0]), local_Nx2[i], MPI_DOUBLE, i, offset_x2[i], local_Nx2[i], MPI_DOUBLE, win);
        }
        MPI_Win_fence(MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOSUCCEED, win);
       
        computeAcceleration(N, r1, a, m, local_N[rank], offset[rank]);


        for (int j = offset[rank]; j < offset[rank]+local_N[rank]; j++)  {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
        }

    }

    PAPI_stop(EventSet, values);
    papi_tot += values[0];

    printf("%lld\n", papi_tot);

    delete[] m;
    delete[] r;
    delete[] u;
    delete[] a;

    MPI_Win_free(&win);
    MPI_Finalize();
    return 0;
}
