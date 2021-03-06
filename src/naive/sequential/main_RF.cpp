#include <math.h>
#include <fstream>
#include <chrono>
#include <iostream>
#include <papi.h>
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
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);



    if (readDataFromFile(params.in_filename, N, m, r, u) == -1) {
        std::cerr << "File " << params.in_filename << " not found!" << std::endl;
        delete[] m;
        delete[] r;
        delete[] u;
        delete[] a;
        return -1;
    }
    params.out_filename = params.in_filename;




    computeAcceleration(N, r, a, m);

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

        computeAcceleration(N, r, a, m);

        for (size_t j = 0; j < N; j++) {
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

    return 0;
}
