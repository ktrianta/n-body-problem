#include <fstream>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "octree.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"


int main(int argc, char** argv) {
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
 
    if (!params.in_filename.empty()) {
        if (readDataFromFile(params.in_filename, N, m, r, u) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
        }
        params.out_filename = params.in_filename;
    } else {
          initializePositionOnSphere(N, r);
    }

    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(N, r, xc, yc, zc, w2, h2, t2);

    std::ofstream out_file;
    openFileToWrite(out_file, params.out_filename, params.out_dirname);
    writeDataToFile(N, r, u, out_file);

    Octree::r = r;
    Octree::a = a;
    Octree::m = m;
    Octree::g = sim::g;
    Octree::theta = theta;
    Octree tree = Octree(N, xc, yc, zc, w2, h2, t2);

    for (size_t j = 0; j < N; j++) {
        tree.computeAcceleration(j);
    }

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

    return 0;
}
