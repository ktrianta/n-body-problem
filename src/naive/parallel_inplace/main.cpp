#include <math.h> //for sqrt function
#include <iostream>
#include <fstream>
#include <vector>
#include "types.hpp"
#include "initialization.hpp"
#include <unistd.h>
#include <mpi.h>

using namespace std;

void computeAcceleration(int N, sim::data_type (*r)[3], sim::data_type (*a)[3], vector<sim::data_type>& m, int rank, int local_N)
{
    for (int i = 0; i < N; i++)
    {
        a[i][0] = 0;
        a[i][1] = 0;
        a[i][2] = 0;
    }

    for (int i = rank*local_N; i < (rank+1)*local_N; i++)
    {
        sim::data_type a_i0 = 0;  // accumulate accelaration values for particle i and
        sim::data_type a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        for (int j = 0; j < N; j++)
        {
            if (i==j) continue;
	
            sim::data_type rji[3];
            rji[0] = r[j][0] - r[i][0];
            rji[1] = r[j][1] - r[i][1];
            rji[2] = r[j][2] - r[i][2];
            sim::data_type r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];
            sim::data_type denom = (r2+sim::e2) * sqrt(r2+sim::e2);
            sim::data_type a_i = - sim::g * m[j] / denom;
            a_i0 -= a_i * rji[0];
            a_i1 -= a_i * rji[1];
            a_i1 -= a_i * rji[2];
        
        }
        a[i][0] += a_i0;  // a(i, 0) and a(i, 1) are accessed once here, avoiding
        a[i][1] += a_i1;  // repeated accesses in the inner loop of j
        a[i][2] += a_i1;  // repeated accesses in the inner loop of j
    }
}

void writeDataToFile(int N, sim::data_type (*r)[3], sim::data_type (*u)[3], ofstream& file)
{
    for (int i = 0; i < N; i++)
    {
        file << r[i][0] << "   "
             << r[i][1] << "   "
             << r[i][2] << "   "
             << u[i][0] << "   "
             << u[i][2] << "   "
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
    int N = 6;
    sim::data_type T = 10;
    sim::data_type dt = 0.00001;
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

    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);
    vector<sim::data_type> m(N, 1.0/N);
    
    if (rank==0)
    {
        if (!filename.empty()) {
        ifstream ifile;
        ifile.open(filename);

            for (int i = 0; i < N; i++) {
                ifile >> m[i] >> r[i][0] >> r[i][1] >> r[i][2] >> u[i][0] >> u[i][1] >> u[i][2];
            }
        } else {
            initializePositionOnSphere(N, r);
        }
    }
    MPI_Bcast (&r[0][0], 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
    MPI_Bcast (&u[0][0], 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
    MPI_Bcast (&m[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    

    ofstream file;
    file.open("output.dat");
    int local_N = N / size; // Asuming size divides N

    if( rank == 0)
    {
    	local_N += N-(size)*local_N;
        writeDataToFile(N, r, u, file);
    }
    computeAcceleration(N, r, a, m, rank, local_N);
    const int Ntimesteps = T/dt + 1;


    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = rank*local_N; j < (rank+1)*local_N ; j++)
        {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
            r[j][0] += u[j][0] * dt;
            r[j][1] += u[j][1] * dt;
            r[j][2] += u[j][2] * dt;
            //printf("%lf %lf %d \n", r[j][0], r[j][1], rank);
		
        }

        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,&(r[0][0]),local_N*3,MPI_DOUBLE,MPI_COMM_WORLD);

        computeAcceleration(N, r, a, m, rank, local_N);

        for (int j = rank*local_N; j < (rank+1)*local_N ; j++)
        {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
        }
	//Needed for correct u in output.dat
        //MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,&(u[0][0]),local_N*2,MPI_DOUBLE,MPI_COMM_WORLD);
        if (rank==0)
        {
            if (t % 200 == 0)
             {
                writeDataToFile(N, r, u, file);
             }   
         }

    }
    MPI_Finalize();
    return 0;
}
