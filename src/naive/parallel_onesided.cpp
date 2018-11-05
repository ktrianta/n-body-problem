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

void computeAcceleration(sim_data_type (*r)[2], sim_data_type (*a)[2], vector<sim_data_type>& m, int rank, int local_N, int offset)
{
    for (int i = 0; i < N; i++)
    {
        a[i][0] = 0;
        a[i][1] = 0;
    }

    for (int i = offset; i < offset + local_N; i++)
    {
        sim_data_type a_i0 = 0;  // accumulate accelaration values for particle i and
        sim_data_type a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        for (int j = 0; j < N; j++)
        {
            if (i==j) continue;
	
            sim_data_type rji[2];
            rji[0] = r[j][0] - r[i][0];
            rji[1] = r[j][1] - r[i][1];
            sim_data_type r2 = rji[0] * rji[0] + rji[1] * rji[1];
            sim_data_type denom = (r2+epsilon2) * sqrt(r2+epsilon2);
            sim_data_type a_i = -g * m[j] / denom;
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
    MPI_Win win;
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
    sim_data_type (*r1)[2] = new sim_data_type[N][2];
    sim_data_type (*u)[2] = new sim_data_type[N][2];
    sim_data_type (*a)[2] = new sim_data_type[N][2];
    std::fill(&u[0][0], &u[0][0] + N*2, 0);
    std::fill(&r1[0][0], &r1[0][0] + N*2, 0);
    std::fill(&a[0][0], &a[0][0] + N*2, 0);
    vector<sim_data_type> m(N, 1.0/N);
    
    if (rank==0)
    {
        if (!filename.empty()) {
        ifstream ifile;
        ifile.open(filename);

            for (int i = 0; i < N; i++) {
                ifile >> m[i] >> r[i][0] >> r[i][1] >> u[i][0] >> u[i][1];
            }
        } else {
            initializePositionOnSphere(N, r);
        }
    }
    MPI_Bcast (&r[0][0], 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
    MPI_Bcast (&u[0][0], 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
    MPI_Bcast (&m[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    

    MPI_Win_create( &(r[0][0]), sizeof(double)*2*N, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    ofstream file;
    file.open("output.dat");

    int local_N[size];
    int local_Nx2[size];
    int local_N_int = N/size;
    int rem = N - local_N_int * size;
    int counter = 0;

    for (int i=0;i<size;i++)
    {
        local_N[i] = local_N_int;
        if (counter < rem) 
        {
                local_N[i] += 1;
                counter ++;
        }
        local_Nx2[i] = local_N[i]*2;
    }
    int offset[size];
    offset[0] = 0;
    int offset_x2[size];
    offset_x2[0] = 0;
    

    for (int i=1; i < size; i++)
    {
        offset[i] = offset[i-1] + local_N[i-1];
        offset_x2[i] = offset_x2[i-1] + local_Nx2[i-1];
    }

    computeAcceleration(r, a, m, rank, local_N[rank], offset[rank]);
    const int Ntimesteps = T/dt + 1;


    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = offset[rank]; j < offset[rank]+local_N[rank]; j++)
        {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            r[j][0] += u[j][0] * dt;
            r[j][1] += u[j][1] * dt;
            //printf("%lf %lf %d \n", r[j][0], r[j][1], rank);
		
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,&(r[0][0]),local_N[rank]*2,MPI_DOUBLE,MPI_COMM_WORLD);
        for (int i = 0; i < size; i++)
        {
            //if( i==rank ) continue;
            //MPI_Win_lock (MPI_LOCK_EXCLUSIVE, i, 0, win);
            MPI_Win_fence(0,win);
            MPI_Get(&(r1[offset[i]][0]), local_Nx2[i], MPI_DOUBLE, i, offset_x2[i], local_Nx2[i], MPI_DOUBLE, win);
            MPI_Win_fence(0,win);
            //MPI_Win_unlock (i, win);
            
        }

        computeAcceleration(r1, a, m, rank, local_N[rank], offset[rank]);

        for (int j = offset[rank]; j < offset[rank]+local_N[rank] ; j++)
        {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
        }
	//Needed for correct u in output.dat
        //MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,&(u[0][0]),local_N*2,MPI_DOUBLE,MPI_COMM_WORLD);
        if (rank==1)
        {
            if (t % 200 == 0)
             {
                writeDataToFile(r1, u, file);
             }   
         }

    }
    MPI_Win_free( &win );
    MPI_Finalize();
    return 0;
}
