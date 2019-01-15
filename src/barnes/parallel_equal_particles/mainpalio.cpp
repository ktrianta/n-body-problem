#include <iostream>
#include <fstream>
#include <vector>
#include "io.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"
#include "serialization.hpp"
#include "sort.hpp"
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
    size_t N = 5;      // the number of particles
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

    int local_N = N/size;
    int *offset = new int[size];
    int *Nloc = new int[size];
    for (int i = 0; i < N; i++)
    {
        offset[i]=i*local_N;
        Nloc[i]=local_N;
    }
    sim::data_type (*r)[3];
    sim::data_type (*u)[3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    sim::data_type *m;
    sim::data_type (*r_local)[3] = new sim::data_type[local_N][3];
    sim::data_type (*u_local)[3] = new sim::data_type[local_N][3];
    sim::data_type (*a_local)[3] = new sim::data_type[local_N][3];
    std::fill(&a_local[0][0], &a_local[0][0] + local_N*3, 0);
    sim::data_type (*m_local) = new sim::data_type[local_N];

    if (rank == 0)
    {
        r = new sim::data_type[N][3];
        m = new sim::data_type[N];
        u = new sim::data_type[N][3];
        if (!filename.empty()) {
            if (readDataFromFile(filename, N, m, r, u) == -1) {
                std::cerr << "File " << filename << " not found!" << std::endl;
                return -1;
            }
        } 
        else {
            initializePositionOnSphere(N, r);
            std::fill(&u[0][0], &u[0][0] + N*3, 0);
            std::fill(m_local, m+local_N, 1.0/local_N);
        }
        p_sort(r, N, size);
        MPI_Scatter(&r[0][0], local_N*3, MPI_DOUBLE, &r_local[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&u[0][0], local_N*3, MPI_DOUBLE, &u_local[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&m[0], local_N, MPI_DOUBLE, &m_local[0], local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Scatter(NULL, local_N*3, MPI_DOUBLE, &r_local[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(NULL, local_N*3, MPI_DOUBLE, &u_local[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(NULL, local_N, MPI_DOUBLE, &m_local[0], local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }


//  the center of the parent node and the half width and height
    double xc, yc, zc, h2, w2, t2;
 
    boxComputation(local_N, r_local, xc, yc, zc, w2, h2, t2);
    std::cout << xc << "   " << w2 << "   "
              << yc << "   " << h2 << "   "
              << zc << "   " << t2 << std::endl;
     xc = 0;
     yc = 0;
     zc = 0;
     w2 = 1;
     h2 = 1;
     t2 = 1;

    ofstream file;
    if (rank == 0)
    {
        file.open("output.dat");
        writeDataToFile(N, r, u, file);
    }
    Serialization *tree[size];
    tree[rank]= new Serialization(xc, yc, zc, w2, h2, t2); 
        for(int j = 0; j < local_N; j++)
            tree[rank]->insert(j, r_local[j][0], r_local[j][1], r_local[j][2], m_local[j]);
    int tSize = tree[rank]->size;
    int treeSize[size];
    MPI_Allgather(&tSize, 1, MPI_INT, &treeSize, 1, MPI_INT, MPI_COMM_WORLD);


    for (int i = 0; i < size; i++)
    {
        if (i !=rank){
            tree[i] = new Serialization();
            MPI_Alloc_mem(treeSize[i]*sizeof(struct Treenode), MPI_INFO_NULL, &(tree[i]->treeArray));
            }
    }
    MPI_Win_attach(win, tree[rank]->treeArray, treeSize[rank]*sizeof(struct Treenode));




    MPI_Aint my_disp;
    MPI_Aint *target_disp = new MPI_Aint[size];
    MPI_Get_address(tree[rank]->treeArray, &my_disp);
    printf("%td \n", my_disp);

    MPI_Allgather(&my_disp, 1, MPI_AINT, &target_disp[0], 1, MPI_AINT, MPI_COMM_WORLD);
    printf("%td , %td \n", target_disp[0], target_disp[1]);


    for (int i = 0; i < size; i++)
    {
        if (i !=rank){
            MPI_Win_fence(0,win);
            MPI_Get(tree[i]->treeArray, treeSize[i], treeNodeStruct, i, target_disp[i], treeSize[i], treeNodeStruct, win);
            MPI_Win_fence(0,win);
           }

        for (int j = 0; j < local_N; j++)
        {
            std::cout << a_local[j][0] << "   " << rank << "  " << i << std::endl;
            tree[i]->computeAcceleration(0, j, r_local, a_local, sim::g, theta);
            std::cout << a_local[j][0] <<  "   " << rank << "  " << i << std::endl;
        }

    }

    MPI_Gather(&r_local[0][0], local_N*3, MPI_DOUBLE, &r[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    const int Ntimesteps = T/dt + 1;

  

 

///*

    double start = 0;
    double end = 0;
    double time = 0;

    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = 0; j < local_N; j++)
        {
            u_local[j][0] += 0.5 * a_local[j][0] * dt;
            u_local[j][1] += 0.5 * a_local[j][1] * dt;
            u_local[j][2] += 0.5 * a_local[j][2] * dt;
            r_local[j][0] += u_local[j][0] * dt;
            r_local[j][1] += u_local[j][1] * dt;
            r_local[j][2] += u_local[j][2] * dt;

            a_local[j][0] = 0;
            a_local[j][1] = 0;
            a_local[j][2] = 0;
        }


        start = MPI_Wtime();
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < local_N; j++)
                {
                    tree[i]->computeAcceleration(0, j, r_local, a_local, sim::g, theta);
                }
        }
            end = MPI_Wtime();
            time += end - start;

            for (int j = 0; j < local_N; j++)
            {

                u_local[j][0] += 0.5 * a_local[j][0] * dt;
                u_local[j][1] += 0.5 * a_local[j][1] * dt;
                u_local[j][2] += 0.5 * a_local[j][2] * dt;
            }

        MPI_Gatherv(&(r_local[0][0]), local_N*3, MPI_DOUBLE, &(r[0][0]), Nloc, offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);


        if (rank == 0)
        {
            p_sort(r, N, size);
            MPI_Scatter(&r[0][0], local_N*3, MPI_DOUBLE, &r_local[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Scatter(&u[0][0], local_N*3, MPI_DOUBLE, &u_local[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Scatter(&m[0], local_N, MPI_DOUBLE, &m_local[0], local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        else{
            MPI_Scatter(NULL, local_N*3, MPI_DOUBLE, &r_local[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Scatter(NULL, local_N*3, MPI_DOUBLE, &u_local[0][0], local_N*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Scatter(NULL, local_N, MPI_DOUBLE, &m_local[0], local_N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
     
        boxComputation(local_N, r_local, xc, yc, zc, w2, h2, t2);
        for (int j = 0; j < size; j++)
        {
            MPI_Win_detach(win,tree[j]->treeArray);
        }
        for (int i = 0; i < size; i++)
        {
            delete tree[i];
        }
        tree[rank]= new Serialization(xc, yc, zc, w2, h2, t2); 
            for(int j = 0; j < local_N; j++)
                tree[rank]->insert(j, r_local[j][0], r_local[j][1], r_local[j][2], m_local[j]);
        tSize = tree[rank]->size;
        MPI_Allgather(&tSize, 1, MPI_INT, &treeSize, 1, MPI_INT, MPI_COMM_WORLD);


        for (int i = 0; i < size; i++)
        {
            if (i !=rank){
                MPI_Alloc_mem(treeSize[i]*sizeof(struct Treenode), MPI_INFO_NULL, &tree[i]->treeArray);
                MPI_Win_attach(win,tree[i]->treeArray,treeSize[i]*sizeof(struct Treenode));
                }
        }



        if (t % 200 == 0 && rank == 0)
        {
            writeDataToFile(N, r, u, file);
        }

    }
    for (int j = 0; j < size; j++)
    {
    MPI_Win_detach(win,tree[j]->treeArray);
    }
    for (int i = 0; i < size; i++)
    {
        delete tree[i];
    }

//*/


    MPI_Win_free(&win);
	MPI_Finalize();
    return 0;
}
