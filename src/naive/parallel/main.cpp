#include <math.h> //for sqrt function
#include <iostream>
#include <fstream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include <mpi.h>

using namespace std;

void computeAcceleration(int N, sim::data_type (*r)[3], sim::data_type (*a_local)[3], sim::data_type *m, int local_N, int offset)
{
    for (int i = 0; i < local_N; i++)
    {
        a_local[i][0] = 0;
        a_local[i][1] = 0;
        a_local[i][2] = 0;
    }

    for (int i = 0; i < local_N; i++)
    {
        sim::data_type a_i0 = 0;  // accumulate accelaration values for particle i and
        sim::data_type a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        sim::data_type a_i2 = 0;  // store them at the end of the loop iteration in a(i,x)
        for (int j = 0; j < N; j++)
        {
            sim::data_type rji[3];
            rji[0] = r[j][0] - r[offset+i][0];
            rji[1] = r[j][1] - r[offset+i][1];
            rji[2] = r[j][2] - r[offset+i][2];
            sim::data_type r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];
            sim::data_type denom = (r2+sim::e2) * sqrt(r2+sim::e2);
            sim::data_type a_i = - sim::g * m[j] / denom;
            a_i0 -= a_i * rji[0];
            a_i1 -= a_i * rji[1];
            a_i2 -= a_i * rji[2];
        }
        a_local[i][0] += a_i0;  // a(i, 0) and a(i, 1) are accessed once here, avoiding
        a_local[i][1] += a_i1;  // repeated accesses in the inner loop of j
        a_local[i][2] += a_i2;  // repeated accesses in the inner loop of j
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
    
    sim::Parameters params;
    readArgs(argc, argv, params);

    const int N = params.n;
    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);

    // PROCESS 0 initialize position vector r.
    if (rank == 0)
    {
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
    MPI_Bcast(&r[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&u[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&m[0],N, MPI_DOUBLE,0, MPI_COMM_WORLD);


    std::ofstream out_file;
    if (rank ==0)
    {
        out_file = openFileToWrite(params.out_filename, params.out_dirname);
        writeDataToFile(params.n, r, u, out_file);
    }

    //if (rank == 0)
    //{
    //}

    // SEND the acceleration vector to all processes.
    // ---- NEED TO BI FILLED -----


    const sim::data_type dt = params.dt;
    const int Ntimesteps = params.t/dt + 1;


    // Local Declarations
    int local_N[size];
    int local_Nx3[size];
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
        local_Nx3[i] = local_N[i]*3;
    }
    
    sim::data_type (*u_local)[3] = new sim::data_type[local_N[rank]][3];
    sim::data_type (*r_local)[3] = new sim::data_type[local_N[rank]][3];
    sim::data_type (*a_local)[3] = new sim::data_type[local_N[rank]][3];

    int offset[size];
    offset[0] = 0;
    int offset_x3[size];
    offset_x3[0] = 0;

    for (int i=1; i < size; i++)
    {
        offset[i] = offset[i-1] + local_N[i-1];
        offset_x3[i] = offset_x3[i-1] + local_Nx3[i-1];
    }

    for (int i=0; i < local_N[rank]; i++)
    {
        u_local[i][0] = u[offset[rank]+i][0];
        u_local[i][1] = u[offset[rank]+i][1];
        u_local[i][2] = u[offset[rank]+i][2];
        r_local[i][0] = r[offset[rank]+i][0];
        r_local[i][1] = r[offset[rank]+i][1];
        r_local[i][2] = r[offset[rank]+i][2];
        a_local[i][0] = a[offset[rank]+i][0];
        a_local[i][1] = a[offset[rank]+i][1];
        a_local[i][2] = a[offset[rank]+i][2];
    }

    computeAcceleration(N, r, a_local, m, local_N[rank], offset[rank]); // NEEDS TO BE PARALLELIZED
    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = 0; j < local_N[rank] ; j++)
        {
            u_local[j][0] += 0.5 * a_local[j][0] * dt;
            u_local[j][1] += 0.5 * a_local[j][1] * dt;
            u_local[j][2] += 0.5 * a_local[j][2] * dt;
            r_local[j][0] += u_local[j][0] * dt;
            r_local[j][1] += u_local[j][1] * dt;
            r_local[j][2] += u_local[j][2] * dt;
        }

		//MPI_Allgather(&(u_local[offset*2][0]),local_N*2,MPI_DOUBLE,&u[0][0],local_N*2,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Allgatherv(&(r_local[0][0]),local_N[rank]*3,MPI_DOUBLE,&(r[0][0]),local_Nx3, offset_x3, MPI_DOUBLE,MPI_COMM_WORLD);

	//MPI_Allgather(&(r_local[0][0]),local_N*2,MPI_DOUBLE,&(r[0][0]),local_N[rank]*2,MPI_DOUBLE,MPI_COMM_WORLD);

        computeAcceleration(N, r, a_local, m, local_N[rank], offset[rank]);
		// SEND the acceleration vector to all processes.
		// ---- NEED TO BE FILLED -----
		

        for (int j = 0 ; j < local_N[rank] ; j++)
        {
		
            u_local[j][0] += 0.5 * a_local[j][0] * dt;
            u_local[j][1] += 0.5 * a_local[j][1] * dt;
            u_local[j][2] += 0.5 * a_local[j][2] * dt;
	}
        //for (int j = 0 ; j < N ; j++)
        //{
        //    u[j][0] += 0.5 * a[j][0] * dt;
        //   u[j][1] += 0.5 * a[j][1] * dt;
        //}

        //MPI_Allgather(&(u_local[offset*2][0]),local_N*2,MPI_DOUBLE,&u[0][0],local_N*2,MPI_DOUBLE,MPI_COMM_WORLD);
        
        if (rank==0)
        {
            if (t % 200 == 0)
             {
                writeDataToFile(N, r, u, out_file);
             }   
         }

    }

    MPI_Finalize();
    return 0;
}
