#include <math.h> //for sqrt function
#include <stdio.h> //for printf
#include <stdlib.h> //for srand/rand prng
#include <time.h> //for seeding prng
#include <iostream>
#include <fstream>
#include <vector>
#include "../utils/types.hpp"
#include "../utils/initialization.hpp"
#include <unistd.h>
#include <mpi.h>

using namespace std;

int N = 6;      // the number of particles
const sim_data_type g = 1;     // gravitational constant
const sim_data_type epsilon = 0.001;
const sim_data_type epsilon2 = epsilon * epsilon;

void computeAcceleration(sim_data_type (*r)[2], sim_data_type (*a)[2], vector<sim_data_type>& m)
{
    for (int i = 0; i < N; i++)
    {
        a[i][0] = 0;
        a[i][1] = 0;
    }

    for (int i = 0; i < N; i++)
    {
        sim_data_type a_i0 = 0;  // accumulate accelaration values for particle i and
        sim_data_type a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        for (int j = i+1; j < N; j++)
        {
            sim_data_type rji[2];
            rji[0] = r[j][0] - r[i][0];
            rji[1] = r[j][1] - r[i][1];
            sim_data_type r2 = rji[0] * rji[0] + rji[1] * rji[1];
            sim_data_type denom = (r2+epsilon2) * sqrt(r2+epsilon2);
            sim_data_type a_j = -g * m[i] / denom;
            sim_data_type a_i = -g * m[j] / denom;
            a[j][0] += a_j * rji[0];
            a[j][1] += a_j * rji[1];
            a_i0 -= a_i * rji[0];
            a_i1 -= a_i * rji[1];
        }
        a[i][0] += a_i0;  // a(i, 0) and a(i, 1) are accessed once here, avoiding
        a[i][1] += a_i1;  // repeated accesses in the inner loop of j
    }
}

void writeDataToFile(sim_data_type (*r)[2], sim_data_type (*u)[2], ofstream& file)
{
    for (int i = 0; i < N; i++)
    {
        file << r[i][0] << "   "
             << r[i][1] << "   "
             << u[i][0] << "   "
             << u[i][1] << "\n";
    }
}

int main(int argc, char** argv)
{

    // *** MPI *** // 
    int size,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // *** MPI *** // 



    int c;
    sim_data_type T = 10;
    sim_data_type dt = 0.00001;
    string filename;

    while ((c = getopt (argc, argv, "n:t:s:i:")) != -1)
    {
        switch (c)
        {
            case 'n':
                N = atoi(optarg);
                break;
            case 't':
                T = atof(optarg);
                break;
            case 's':
                dt = atof(optarg);
                break;
            case 'i':
                filename = optarg;
                break;
        }
    }

    sim_data_type (*r)[2] = new sim_data_type[N][2];
    sim_data_type (*u)[2] = new sim_data_type[N][2];
    sim_data_type (*a)[2] = new sim_data_type[N][2];
    std::fill(&u[0][0], &u[0][0] + N*2, 0);
    std::fill(&a[0][0], &a[0][0] + N*2, 0);

    vector<sim_data_type> m(N, 1.0/N);

    if (!filename.empty()) {
        ifstream ifile;
        ifile.open(filename);

        for (int i = 0; i < N; i++) {
            ifile >> m[i] >> r[i][0] >> r[i][1] >> u[i][0] >> u[i][0];
        }
    } else {
        initializePositionOnSphere(N, r);
    }

    ofstream file;
    file.open("output.dat");

    writeDataToFile(r, u, file);
    computeAcceleration(r, a, m);
    const int Ntimesteps = T/dt + 1;

    int local_N = N / size; // Asuming size divides N

    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = rank*local_N; j < (rank+1)*local_N ; j++)
        {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            r[j][0] += u[j][0] * dt;
            r[j][1] += u[j][1] * dt;
        }

        MPI_Allgather(&(u[rank*local_N*2][0]),local_N*2,MPI_DOUBLE,&u[0][0],local_N*2,MPI_DOUBLE,MPI_COMM_WORLD);

        computeAcceleration(r, a, m);

        for (int j = rank*local_N; j < (rank+1)*local_N ; j++)
        {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
        }

        MPI_Allgather(&(u[rank*local_N*2][0]),local_N*2,MPI_DOUBLE,&u[0][0],local_N*2,MPI_DOUBLE,MPI_COMM_WORLD);
        if (rank==0)
        {
            if (t % 200 == 0)
             {
                writeDataToFile(r, u, file);
             }   
         }

    }
    MPI_Finalize();
    return 0;
}
