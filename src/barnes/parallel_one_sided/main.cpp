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

using namespace std;

int main(int argc, char** argv) {
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

    sim::Parameters params;
    readArgs(argc, argv, params);

    const size_t N = params.n;
    const sim::data_type theta = params.theta;
    int local_N = N/size;

    MPI_Datatype vtype;
    MPI_Type_vector(1, 4, 7, MPI_DOUBLE, &vtype);
    MPI_Type_commit(&vtype);

    sim::data_type (*r)[7];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    sim::data_type (*r_local)[7] = new sim::data_type[local_N][7];
    sim::data_type (*a_local)[3] = new sim::data_type[local_N][3];
    std::fill(&a_local[0][0], &a_local[0][0] + local_N*3, 0);

    if (rank == 0) {
        r = new sim::data_type[N][7];

        if (readDataFromFile(params.in_filename, N, r) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            return -1;
        }

        params.out_filename = params.in_filename;
        MPI_Scatter(&r[0][0], local_N*7, MPI_DOUBLE, &r_local[0][0], local_N*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else{
        MPI_Scatter(NULL, local_N*7, MPI_DOUBLE, &r_local[0][0], local_N*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(local_N, r_local, xc, yc, zc, w2, h2, t2);
     
    std::ofstream out_file;
    if (rank == 0) {
        openFileToWrite(out_file, params.out_filename, params.out_dirname);
        writeDataToFile(N, r, out_file);
    }

    Serialization *tree[size];
    tree[rank] = new Serialization(xc, yc, zc, w2, h2, t2);

    for (int j = 0; j < local_N; j++) {
        tree[rank]->insert(rank*local_N+j, r_local[j][1], r_local[j][2], r_local[j][3], r_local[j][0]);
    }

    int tSize = tree[rank]->position+1;
    int treeSize[size];
    MPI_Allgather(&tSize, 1, MPI_INT, &treeSize, 1, MPI_INT, MPI_COMM_WORLD);

    for (int i = 0; i < size; i++) {
        if (i != rank) {
            tree[i] = new Serialization();
            MPI_Alloc_mem(treeSize[i]*sizeof(struct Treenode), MPI_INFO_NULL, &(tree[i]->treeArray));
        }
    }

    MPI_Win_attach(win, tree[rank]->treeArray, treeSize[rank]*sizeof(struct Treenode));
    MPI_Aint my_disp;
    MPI_Aint *target_disp = new MPI_Aint[size];
    MPI_Get_address(tree[rank]->treeArray, &my_disp);
    MPI_Allgather(&my_disp, 1, MPI_AINT, &target_disp[0], 1, MPI_AINT, MPI_COMM_WORLD);

    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Win_fence(0,win);
            MPI_Get(tree[i]->treeArray, treeSize[i], treeNodeStruct, i, target_disp[i], treeSize[i], treeNodeStruct, win);
            MPI_Win_fence(0,win);
        }

        for (int j = 0; j < local_N; j++) {
            tree[i]->computeAcceleration(rank*local_N+j, &r_local[j][1], a_local[j], sim::g, theta);
        }
    }

    const size_t Ntimesteps = params.t / params.dt + 1;
    const sim::data_type dt = params.dt;

    for (int t = 0; t < Ntimesteps; t++) {
        for (int j = 0; j < local_N; j++) {
            r_local[j][4] += 0.5 * a_local[j][0] * dt;
            r_local[j][5] += 0.5 * a_local[j][1] * dt;
            r_local[j][6] += 0.5 * a_local[j][2] * dt;
            r_local[j][1] += r_local[j][4] * dt;
            r_local[j][2] += r_local[j][5] * dt;
            r_local[j][3] += r_local[j][6] * dt;

            a_local[j][0] = 0;
            a_local[j][1] = 0;
            a_local[j][2] = 0;
        }

        MPI_Gather(r_local, local_N*7, MPI_DOUBLE, r, local_N*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            p_sort(r, N, size);
            MPI_Scatter(&r[0][0], local_N*7, MPI_DOUBLE, &r_local[0][0], local_N*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else{
            MPI_Scatter(NULL, local_N*7, MPI_DOUBLE, &r_local[0][0], local_N*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
     
        MPI_Win_detach(win,tree[rank]->treeArray);
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                MPI_Free_mem(tree[i]->treeArray);
                tree[i]->treeArray = NULL;
            }
            delete tree[i];
        }

        boxComputation(local_N, r_local, xc, yc, zc, w2, h2, t2);
        tree[rank] = new Serialization(xc, yc, zc, w2, h2, t2);
        for(int j = 0; j < local_N; j++) {
            tree[rank]->insert(rank*local_N+j, r_local[j][1], r_local[j][2], r_local[j][3], r_local[j][0]);
        }

        tSize = tree[rank]->size;
        MPI_Allgather(&tSize, 1, MPI_INT, &treeSize, 1, MPI_INT, MPI_COMM_WORLD);

        for (int i = 0; i < size; i++) {
            if (i !=rank) {
                tree[i] = new Serialization();
                MPI_Alloc_mem(treeSize[i] * sizeof(struct Treenode), MPI_INFO_NULL, &(tree[i]->treeArray));
            }
        }

        MPI_Win_attach(win, tree[rank]->treeArray, treeSize[rank]*sizeof(struct Treenode));
        MPI_Get_address(tree[rank]->treeArray, &my_disp);
        MPI_Allgather(&my_disp, 1, MPI_AINT, target_disp, 1, MPI_AINT, MPI_COMM_WORLD);

        for (int i = rank; i < rank+size; i++) {
            MPI_Win_fence(0, win);
            if (i+1 != rank) {
                MPI_Get(tree[(i+1)%size]->treeArray, treeSize[(i+1)%size], treeNodeStruct, (i+1)%size, target_disp[(i+1)%size], treeSize[(i+1)%size], treeNodeStruct, win);
            }

            for (int j = 0; j < local_N; j++) {
                tree[i%size]->computeAcceleration(rank*local_N+j, &r_local[j][1], a_local[j], sim::g, theta);
            }
            MPI_Win_fence(0, win);
        }

        for (int j = 0; j < local_N; j++) {
            r_local[j][4] += 0.5 * a_local[j][0] * dt;
            r_local[j][5] += 0.5 * a_local[j][1] * dt;
            r_local[j][6] += 0.5 * a_local[j][2] * dt;
        }

        if (rank == 0 && t % 200 == 0) {
            writeDataToFile(N, r, out_file);
        }

    }

    MPI_Win_detach(win,tree[rank]->treeArray);
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Free_mem(tree[i]->treeArray);
            tree[i]->treeArray = NULL;
        }
        delete tree[i];
    }

    MPI_Win_free(&win);
	MPI_Finalize();
    return 0;
}
