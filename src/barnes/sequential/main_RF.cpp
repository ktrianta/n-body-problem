#include <math.h>
#include <fstream>
#include <iostream>
#include <papi.h>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "octree.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"
#include <chrono>


using time_point_t = std::chrono::high_resolution_clock::time_point;

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
        printf("PAPI add PAPI_L1_DCM error!\n");
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
    std::ofstream out_file;
    openFileToWrite(out_file, params.out_filename, params.out_dirname);
    writeDataToFile(N, r, u, out_file);


    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(N, r, xc, yc, zc, w2, h2, t2);

    Octree::r = r;
    Octree::a = a;
    Octree::m = m;
    Octree::g = sim::g;
    Octree::theta = params.theta;
    Octree tree = Octree(N, xc, yc, zc, w2, h2, t2);


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
        Octree tree = Octree(N, xc, yc, zc, w2, h2, t2);
        std::fill(&a[0][0], &a[0][0] + N*3, 0);

        for (size_t j = 0; j < N; j++) {
            tree.computeAcceleration(j);

            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
        }

        if (t % 200 == 0) {
            writeDataToFile(N, r, u, out_file);
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
