#include <iostream>
#include <fstream>
#include <vector>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"
#include "serialization.hpp"
#include "sort.hpp"
#include <unistd.h>
#include <mpi.h>
#include <math.h>

using namespace std;

int main(int argc, char** argv)
{
    // ----- MPI ----- //
    int size,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    sim::Parameters params;
    readArgs(argc, argv, params);

    const size_t N = params.n;
    const sim::data_type theta = params.theta;
    int local_N = N/size;

    sim::data_type (*r)[7] = new sim::data_type[N][7];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(&a[0][0], &a[0][0] + N*3, 0);

    double initialKEnergy = 0;
    double initialPEnergy = 0;
    double initialEnergy = 0;

    if (rank == 0) {
        r = new sim::data_type[N][7];

        if (readDataFromFile(params.in_filename, N, r) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            return -1;
              p_sort(r, N, size);
        }

        params.out_filename = params.in_filename;

        for (int i = 0; i < N; i++){
            initialKEnergy += r[i][0] * (r[i][1]*r[i][1] + r[i][2]*r[i][2] + r[i][3]*r[i][3])/2.;
            for (int j = 0; j < i; j++){
                double denominator = sqrt((r[j][0]*r[j][1]-r[i][0]*r[i][1])*(r[j][0]*r[j][1]-r[i][0]*r[i][1]) + (r[j][0]*r[j][2]-r[i][0]*r[i][2])*(r[j][0]*r[j][2]-r[i][0]*r[i][2]) +
                                          (r[j][0]*r[j][3]-r[i][0]*r[i][3])*(r[j][0]*r[j][3]-r[i][0]*r[i][3]));
                initialPEnergy -= sim::g*r[i][0]*r[j][0]/denominator;
                }
            }
            initialEnergy = initialKEnergy + initialPEnergy;

    }
    // SEND the position vector r from Process 0 to all processes.
    MPI_Bcast(&r[0][0],N*7, MPI_DOUBLE,0, MPI_COMM_WORLD);

//  the center of the parent node and the half width and height
    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(local_N, rank*local_N, (rank+1)*local_N, r, xc, yc, zc, w2, h2, t2);



    std::ofstream out_file;
    if (rank == 0) {
        openFileToWrite(out_file, params.out_filename, params.out_dirname);
        writeDataToFile(N, r, out_file);
    }
   
    sim::data_type (*a_local)[3] = new sim::data_type[N][3];
    std::fill(&a_local[0][0], &a_local[0][0] + N*3, 0);

    
 
    Serialization *tree = new Serialization(xc, yc, zc, w2, h2, t2);
        for(int j = rank*local_N; j < (rank+1)*local_N; j++)
        {
            tree->insert(j, r[j][1], r[j][2], r[j][3], r[j][0]);
        }


    const size_t Ntimesteps = params.t / params.dt + 1;
    const sim::data_type dt = params.dt;
    

        for (int j = 0; j < N; j++)
        {
        //    tree->computeAcceleration(0, j, r, a_local, sim::g, theta);
            tree->computeAcceleration(j, &r[j][1], a_local[j], sim::g, theta);
        }
//      for (int j = 0; j < N; j++)
//      {
//      //    tree->computeAcceleration(0, j, r, a_local, sim::g, theta);
//      std::cout << a_local[j][0] << std::endl;
//      }
   
 
    MPI_Allreduce(a_local,a,3*N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//      for (int j = 0; j < N; j++)
//      {
//      //    tree->computeAcceleration(0, j, r, a_local, sim::g, theta);
//      std::cout << a[j][0] << std::endl;
//      }


    double start = 0;
    double end = 0;
    double comm = 0;
    double comp = 0;

    for (int t = 0; t < Ntimesteps; t++)
//  for (int t = 0; t < 0; t++)
    {
        for (int j = 0; j < N; j++)
        {
            r[j][4] += 0.5 * a[j][0] * dt;
            r[j][5] += 0.5 * a[j][1] * dt;
            r[j][6] += 0.5 * a[j][2] * dt;
            r[j][1] += r[j][4] * dt;
            r[j][2] += r[j][5] * dt;
            r[j][3] += r[j][6] * dt;

            a_local[j][0] = 0;
            a_local[j][1] = 0;
            a_local[j][2] = 0;
            a[j][0] = 0;
            a[j][1] = 0;
            a[j][2] = 0;
        }
        p_sort(r, N, size);
        boxComputation(local_N, rank*local_N, (rank+1)*local_N, r, xc, yc, zc, w2, h2, t2);

        delete tree;
        tree= new Serialization(xc, yc, zc, w2, h2, t2);
        for(int j = rank*local_N; j < (rank+1)*local_N; j++)
        {
            tree->insert(j, r[j][1], r[j][2], r[j][3], r[j][0]);
        }

        start = MPI_Wtime();
        for (int j = 0; j < N; j++)
        {
        //    tree->computeAcceleration(0, j, r, a_local, sim::g, theta);
            tree->computeAcceleration(j, &r[j][1], a_local[j], sim::g, theta);
        }
        end = MPI_Wtime();
        comp += end - start;
        
        start = MPI_Wtime();
        MPI_Allreduce(a_local,a,3*N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        end = MPI_Wtime();
        comm += end - start;
        for (int j = 0; j < N; j++)
        {

            r[j][4] += 0.5 * a[j][0] * dt;
            r[j][5] += 0.5 * a[j][1] * dt;
            r[j][6] += 0.5 * a[j][2] * dt;
        }


        if (t % 200 == 0 && rank == 0)
        {
            writeDataToFile(N, r, out_file);
        }

    }
    std::cout << comm << "   " << comp << std::endl;
	MPI_Finalize();
    return 0;
}
