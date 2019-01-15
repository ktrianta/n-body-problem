#include <iostream>
#include <fstream>
#include <vector>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"
#include "serialization.hpp"
#include <unistd.h>
#include <mpi.h>

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
    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);
    vector<vector<int> > rdecomposed;
    rdecomposed.resize(64);
    sim::data_type xcenter[64], ycenter[64], zcenter[64];
    int pNumber[64] = {0};

    if (rank == 0) {
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

//  the center of the parent node and the half width and height
    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(N, r, xc, yc, zc, w2, h2, t2);

    xc = 0.;
    yc = 0;
    zc = 0;
    w2 = 50;
    h2 = 50;
    t2 = 50;

    // split the domain in 64 parts
    for (int k = 0; k < 4; k++)
    { 
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
            xcenter[j+4*i+16*k] = xc + (2.*i-3.)*w2/4.;
            ycenter[j+4*i+16*k] = yc + (2.*j-3.)*h2/4.;
            zcenter[j+4*i+16*k] = zc + (2.*k-3.)*t2/4.;
                for (int n =0; n < N; n++)
                {
                    if (r[n][0] <= (xcenter[j+4*i+16*k] + w2/4) && (r[n][0] > xcenter[j+4*i+16*k] - w2/4) &&
                        r[n][1] <= (ycenter[j+4*i+16*k] + h2/4) && (r[n][1] > ycenter[j+4*i+16*k] - h2/4) &&
                        r[n][2] <= (zcenter[j+4*i+16*k] + t2/4) && (r[n][2] > zcenter[j+4*i+16*k] - t2/4))
                    {
                        rdecomposed[j+4*i+16*k].push_back(n); 
                        pNumber[j+4*i+16*k]++; 
                    }
                }
            }
        }
    }

    

    std::ofstream out_file;
    if (rank == 0) {
        openFileToWrite(out_file, params.out_filename, params.out_dirname);
        writeDataToFile(N, r, u, out_file);
    }
   
    int *localSubdomains = new int[size];
    int  subs = 64/size;
    int counter = 0;
    int res = 64 - subs*size;
    for (int i=0;i<size;i++)
    {
        localSubdomains[i] = subs;
        if (counter < res)
        {
                localSubdomains[i] += 1;
                counter ++;
        }
    }
    int *offset = new int[size];
    offset[0] = 0;
    for (int i = 1; i < size; i++)
    {
        offset[i] = offset[i-1] + localSubdomains[i-1];
    }
    
    //Particles counter
    int *localNumber = new int[size];
    int *offsetParticles = new int[localSubdomains[rank]];
    offsetParticles[0] =0;

    for (int i = 0; i < size; i++)
    {
        localNumber[i] = 0;
    }
    for (int i = 0; i < localSubdomains[rank]; i++)
    {
        localNumber[rank]+= pNumber[i];
    }
    for (int i = 1; i < localSubdomains[rank]; i++)
    {
        offsetParticles[i] = offsetParticles[i-1]+ pNumber[i-1 + offset[rank]];
    }
    
//  sim::data_type (*u_local)[3] = new sim::data_type[N][3];
//  sim::data_type (*r_local)[3] = new sim::data_type[N][3];
    sim::data_type (*a_local)[3] = new sim::data_type[N][3];
    std::fill(&a_local[0][0], &a_local[0][0] + N*3, 0);

    
// Construction of trees for every rank

 
    Serialization *tree[localSubdomains[rank]]; 
    for(int i = offset[rank]; i < localSubdomains[rank] +offset[rank]; i++)
    { 
        tree[i-offset[rank]]= new Serialization(xcenter[i], ycenter[i], zcenter[i], w2/4,h2/4,t2/4);
        for(int j = 0; j < pNumber[i]; j++)
            if(pNumber[i] != 0)
                tree[i-offset[rank]]->insert(rdecomposed[i][j], r[rdecomposed[i][j]][0], r[rdecomposed[i][j]][1], r[rdecomposed[i][j]][2], m[rdecomposed[i][j]]);
    }
  

    const size_t Ntimesteps = params.t / params.dt + 1;
    const sim::data_type dt = params.dt;
    

    for (int i = offset[rank]; i < localSubdomains[rank] +offset[rank]; i++)
    {    
        for (int j = 0; j < N; j++)
        {
            if(pNumber[i] != 0)
            {
                tree[i-offset[rank]]->computeAcceleration(0, j, r, a_local, sim::g, theta);
            }
        }
    }

 
    MPI_Allreduce(a_local,a,3*N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    double start = 0;
    double end = 0;
    double time = 0;
//  /*   
    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = 0; j < N; j++)
        {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
            r[j][0] += u[j][0] * dt;
            r[j][1] += u[j][1] * dt;
            r[j][2] += u[j][2] * dt;

            a_local[j][0] = 0;
            a_local[j][1] = 0;
            a_local[j][2] = 0;
            a[j][0] = 0;
            a[j][1] = 0;
            a[j][2] = 0;
        }


        for (int i = offset[rank]; i < localSubdomains[rank] + offset[rank]; i++)
        {
            start = MPI_Wtime();
            for (int j = 0; j < N; j++)
            {
                if(pNumber[i] != 0){
                    tree[i-offset[rank]]->computeAcceleration(0, j, r, a_local, sim::g, theta);
                    }

            }
            end = MPI_Wtime();
            time += end - start;

         }
        MPI_Allreduce(a_local,a,3*N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int j = 0; j < N; j++)
        {

            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
        }
        for (int i = 0; i < 64; i++)
        {
        pNumber[i] = 0;
        rdecomposed[i].clear();
        
       
        for (int n =0; n < N; n++)
        {
            if (r[n][0] <= (xcenter[i] + w2/4) && (r[n][0] > xcenter[i] - w2/4) &&
                r[n][1] <= (ycenter[i] + h2/4) && (r[n][1] > ycenter[i] - h2/4) &&
                r[n][2] <= (zcenter[i] + t2/4) && (r[n][2] > zcenter[i] - t2/4))
            {
                rdecomposed[i].push_back(n); 
                pNumber[i]++; 
            }
        }
        }
        for(int i = offset[rank]; i < localSubdomains[rank] +offset[rank]; i++)
        {
            delete tree[i-offset[rank]];
            tree[i-offset[rank]]= new Serialization(xcenter[i], ycenter[i], zcenter[i], w2/4,h2/4,t2/4);
            for(int j = 0; j < pNumber[i]; j++)
                if(pNumber[i] != 0){
                    tree[i-offset[rank]]->insert(rdecomposed[i][j], r[rdecomposed[i][j]][0], r[rdecomposed[i][j]][1], r[rdecomposed[i][j]][2], m[rdecomposed[i][j]]);
                    }
        }

        if (t % 200 == 0 && rank == 0)
        {
            writeDataToFile(N, r, u, out_file);
        }

    }
// */

	MPI_Finalize();
    return 0;
}
