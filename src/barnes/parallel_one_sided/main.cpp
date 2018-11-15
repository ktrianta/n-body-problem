#include <iostream>
#include <fstream>
#include <vector>
#include "io.hpp"
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

	Treenode defTree;
	MPI_Datatype treeNodeStruct;
	int blocklength[] = {10,1,1,2};
	MPI_Datatype old_types[] = {MPI_DOUBLE,MPI_INT,MPI_C_BOOL,MPI_UNSIGNED_LONG};
	MPI_Aint baseaddr,a1,a2,a3,a4;
	MPI_Get_address(&defTree,&baseaddr);
	MPI_Get_address(&defTree.x, &a1);
	MPI_Get_address(&defTree.index,&a2);
	MPI_Get_address(&defTree.leaf,&a3);
	MPI_Get_address(&defTree.cum_size,&a4);
	MPI_Aint indices[] = {a1-baseaddr,a2-baseaddr,a3-baseaddr,a4-baseaddr};
	MPI_Type_create_struct(4,blocklength,indices,old_types,&treeNodeStruct);
	MPI_Type_commit(&treeNodeStruct);


	MPI_Win win;
	MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &win);

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

//  the center of the parent node and the half width and height
    double xc, yc, zc, h2, w2, t2;
 
    boxComputation(N, r, xc, yc, zc, w2, h2, t2);

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
    }
    
    sim::data_type (*u_local)[3] = new sim::data_type[localNumber[rank]][3];
    sim::data_type (*r_local)[3] = new sim::data_type[localNumber[rank]][3];
    sim::data_type (*a_local)[3] = new sim::data_type[localNumber[rank]][3];


// Construction of trees for every rank

    double sum = 0; // is the number of nodes in all trees of a processor
 
    Serialization *tree[localSubdomains[rank]]; 
    for(int i = offset[rank]; i < localSubdomains[rank] +offset[rank]; i++)
    {
        tree[i]= new Serialization(xcenter[i], ycenter[i], zcenter[i], w2/4,h2/4,t2/4);
        for(int j = 0; j < pNumber[i]; j++)
            tree[i]->insert(rdecomposed[i][j], r[j][0], r[j][1], r[j][2], m[j]);
            sum += tree[i]->size;
    }
  

    MPI_Aint my_displ[localSubdomains[rank]];
    MPI_Aint ** target_disp = new MPI_Aint*[localSubdomains[rank]];

    for (int i = 0; i < localSubdomains[rank] ;i++){
        target_disp[i] = new MPI_Aint[tree[i]->size];
        }

    for(int i = 0; i < localSubdomains[rank]; i++)
    {
//      MPI_Alloc_mem(tree[i]->size*sizeof(struct Treenode), MPI_INFO_NULL, &tree[i]->treeArray);
        MPI_Win_attach(win,tree[i]->treeArray,tree[i]->size*sizeof(struct Treenode));

        MPI_Get_address(tree[i]->treeArray, &my_displ[i]); 
        MPI_Allgather(&my_displ[i], localSubdomains[rank], MPI_AINT, target_disp, localSubdomains[rank], MPI_AINT, MPI_COMM_WORLD);
    }   

    const int Ntimesteps = T/dt + 1;


    for (int i = offset[rank]; i < localSubdomains[rank] +offset[rank]; i++)
    {
        for (int j = 0; j < pNumber[i]; j++)
        {
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
    for(int i = 0; i < localSubdomains[rank]; i++)
    {
        MPI_Win_detach(win,tree[i]->treeArray);
        free(tree[i]->treeArray);
    }   
    free(target_disp);
	MPI_Win_free(&win);
	MPI_Type_free(&treeNodeStruct);
	MPI_Finalize();
    return 0;
}
