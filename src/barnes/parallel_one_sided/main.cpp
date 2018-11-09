#include <iostream>
#include <fstream>
#include <vector>
#include "io.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include "octree.hpp"
#include <unistd.h>
#include <mpi.h>
#

using namespace std;

void writeDataToFile(int N, sim::data_type (*r)[3], sim::data_type (*u)[3], ofstream& file)
{
    for (int i = 0; i < N; i++)
    {
        file << r[i][0] << "   "
             << r[i][1] << "   "
             << r[i][2] << "   "
             << u[i][0] << "   "
             << u[i][1] << "   "
             << u[i][2] << "\n";
    }
}

int main(int argc, char** argv)
{
    // ----- MPI ----- //
    int size,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//  the center of the parent node and the half width and height
    double xc, yc, zc, h2, w2, t2;
// coordinates for tab8096
    xc = 0;
    yc = 0;
    zc = 0;
    w2 = 30;
    h2 = 30;
    t2 = 30;

    int c;
    int N = 5;      // the number of particles
    double theta = 0.0;      // the number of particles
    sim::data_type T = 10;
    sim::data_type dt = 0.00001;
    string filename;

    while ((c = getopt (argc, argv, "n:t:s:i:h:")) != -1)
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
            case 'h':
                theta = atof(optarg);
                break;
        }
    }

    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);
//  sim::data_type *rdecomposed[64];
    vector<vector<int> > rdecomposed;
    rdecomposed.resize(64);
    sim::data_type xcenter[64], ycenter[64], zcenter[64];
    int pNumber[64] = {0};

    if (rank == 0)
    {
        if (!filename.empty()) {
            if (readDataFromFile(filename, N, m, r, u) == -1) {
                std::cerr << "File " << filename << " not found!" << std::endl;
                return -1;
            }
        } else {
              initializePositionOnSphere(N, r);
        }
    }
    // SEND the position vector r from Process 0 to all processes.
    MPI_Bcast(&r[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&u[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&m[0],N, MPI_DOUBLE,0, MPI_COMM_WORLD);


    // split the domain in 64 parts
//  if (rank == 0)
//  {
        for (int k = 0; k < 4; k++)
        { 
            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                xcenter[j+4*i+16*k] = xc + (2.*i-3.)*w2/4.;
                ycenter[j+4*i+16*k] = yc + (2.*j-3.)*h2/4.;
                zcenter[j+4*i+16*k] = zc + (2.*k-3.)*t2/4.;
//              printf(" x = %f y = %f z= %f \n" , xcenter[j+4*i+16*k], ycenter[j+4*i+16*k], zcenter[j+4*i+16*k]);
                    for (int n =0; n < N; n++)
                    {
                        if (r[n][0] <= (xcenter[j+4*i+16*k] + w2/4) && (r[n][0] > xcenter[j+4*i+16*k] - w2/4) &&
                            r[n][1] <= (ycenter[j+4*i+16*k] + h2/4) && (r[n][1] > ycenter[j+4*i+16*k] - h2/4) &&
                            r[n][2] <= (zcenter[j+4*i+16*k] + t2/4) && (r[n][2] > zcenter[j+4*i+16*k] - t2/4))
                        {
                            rdecomposed[j+4*i+16*k].push_back(n); 
//                          printf("particle tad %d and how many particles %d \n", rdecomposed[j+4*i+16*k][pNumber[j+4*i+16*k]],j+4*i+16*k );  
                            pNumber[j+4*i+16*k]++; 
                        }
                    }
//                  printf(" print number in tree %d, %d \n", pNumber[j+4*i+16*k], rank);
                }
            }
        }
//  }

    

    ofstream file;
    if (rank == 0)
    {
        file.open("output.dat");
        writeDataToFile(N, r, u, file);
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
//        printf("%d  %d \n", offsetParticles[i], rank); 
    }
    


    Octree *tree[localSubdomains[rank]]; 
    sim::data_type (*u_local)[3] = new sim::data_type[localNumber[rank]][3];
    sim::data_type (*r_local)[3] = new sim::data_type[localNumber[rank]][3];
    sim::data_type (*a_local)[3] = new sim::data_type[localNumber[rank]][3];

    for(int i = offset[rank]; i < localSubdomains[rank] +offset[rank]; i++)
    {
        tree[i]= new Octree(r, m, pNumber[i], xcenter[i], ycenter[i], zcenter[i], w2/4,h2/4,t2/4);
//        printf(" i = %d %d \n", i, rank); 
    }
   
//    Octree tree = Octree(r, m, N, xc, yc, zc, w2, h2, t2);

    const int Ntimesteps = T/dt + 1;


    for (int i = offset[rank]; i < localSubdomains[rank] +offset[rank]; i++)
    {
        for (int j = 0; j < pNumber[i]; j++)
        {
            printf("%d , %d , %d  \n", rank, j+offsetParticles[i-offset[rank]], rdecomposed[i][j]);
//          printf("%d , %d , %d %d \n", rank, i, j, offsetParticles[i]);
            u_local[j+offsetParticles[i-offset[rank]]][0] = u[rdecomposed[i][j]][0];
            u_local[j+offsetParticles[i-offset[rank]]][1] = u[rdecomposed[i][j]][1];
            u_local[j+offsetParticles[i-offset[rank]]][2] = u[rdecomposed[i][j]][2];
            r_local[j+offsetParticles[i-offset[rank]]][0] = r[rdecomposed[i][j]][0];
            r_local[j+offsetParticles[i-offset[rank]]][1] = r[rdecomposed[i][j]][1];
            r_local[j+offsetParticles[i-offset[rank]]][2] = r[rdecomposed[i][j]][2];
            a_local[j+offsetParticles[i-offset[rank]]][0] = a[rdecomposed[i][j]][0];
            a_local[j+offsetParticles[i-offset[rank]]][1] = a[rdecomposed[i][j]][1];
            a_local[j+offsetParticles[i-offset[rank]]][2] = a[rdecomposed[i][j]][2];
        }    
    }

//  for (int j = offset[rank]; j < local_N[rank] + offset[rank]; j++)
//  {
//      int jLocal=j-offset[rank];
//      tree.computeAcceleration(j,jLocal, r, a_local, sim::g, theta); // NEEDS TO BE PARALLELIZED
//  }

/*
    double start = 0;
    double end = 0;
    double time = 0;

    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = offset[rank]; j < local_N[rank] + offset[rank]; j++)
        {
            int jLocal=j-offset[rank];
            u_local[jLocal][0] += 0.5 * a_local[jLocal][0] * dt;
            u_local[jLocal][1] += 0.5 * a_local[jLocal][1] * dt;
            u_local[jLocal][2] += 0.5 * a_local[jLocal][2] * dt;
            r_local[jLocal][0] += u_local[jLocal][0] * dt;
            r_local[jLocal][1] += u_local[jLocal][1] * dt;
            r_local[jLocal][2] += u_local[jLocal][2] * dt;
            a_local[jLocal][0] = 0;
            a_local[jLocal][1] = 0;
            a_local[jLocal][2] = 0;
        }

        MPI_Allgatherv(&(r_local[0][0]),local_N[rank]*3,MPI_DOUBLE,&(r[0][0]),local_Nx3, offset_x3, MPI_DOUBLE,MPI_COMM_WORLD);

        for (int j = offset[rank]; j < local_N[rank] + offset[rank]; j++)
        {
            int jLocal=j-offset[rank];
            start = MPI_Wtime();
            tree.computeAcceleration(j, jLocal, r, a_local, sim::g, theta);
            end = MPI_Wtime();
            time += end - start;

            u_local[jLocal][0] += 0.5 * a_local[jLocal][0] * dt;
            u_local[jLocal][1] += 0.5 * a_local[jLocal][1] * dt;
            u_local[jLocal][2] += 0.5 * a_local[jLocal][2] * dt;
        }

        Octree tree = Octree(r, m, N, xc, yc, zc, w2, h2, t2);
        if (t % 200 == 0 && rank == 0)
        {
            writeDataToFile(N, r, u, file);
        }
    }

    printf("time= %f \n", time);
    */
    MPI_Finalize();
    return 0;
}
