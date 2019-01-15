#include <iostream>
#include <vector>
#include "serialization.hpp"
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

    const size_t N =1000;
    int local_N = N/size;

    MPI_Datatype vtype;
    MPI_Type_vector(1, 4, 7, MPI_DOUBLE, &vtype);
    MPI_Type_commit(&vtype);

    sim::data_type (*r)[7];
    sim::data_type (*r_local)[7] = new sim::data_type[local_N][7];

    if (rank == 0) {
        r = new sim::data_type[N][7];
        MPI_Scatter(&r[0][0], local_N*7, MPI_DOUBLE, &r_local[0][0], local_N*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else{
        MPI_Scatter(NULL, local_N*7, MPI_DOUBLE, &r_local[0][0], local_N*7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    Serialization *tree[size];
    int tSize = 65000;
    for (int i = 0; i < size; i++) {
            tree[i] = new Serialization();
            MPI_Alloc_mem(tSize*sizeof(struct Treenode), MPI_INFO_NULL, &(tree[i]->treeArray));
    }

    MPI_Win_attach(win, tree[rank]->treeArray, tSize*sizeof(struct Treenode));
    MPI_Aint my_disp;
    MPI_Aint *target_disp = new MPI_Aint[size];
    MPI_Get_address(tree[rank]->treeArray, &my_disp);
    MPI_Allgather(&my_disp, 1, MPI_AINT, &target_disp[0], 1, MPI_AINT, MPI_COMM_WORLD);

    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Win_fence(0,win);
            MPI_Get(tree[i]->treeArray, tSize, treeNodeStruct, i, target_disp[i], tSize, treeNodeStruct, win);
            MPI_Win_fence(0,win);
        }

    }


    double t1_start;
    double t1_end;
    double time1;

    for (int t = 0; t < 100; t++) {
        MPI_Win_detach(win,tree[rank]->treeArray);
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                MPI_Free_mem(tree[i]->treeArray);
                tree[i]->treeArray = NULL;
            }
            delete tree[i];
        }

        for (int i = 0; i < size; i++) {
                tree[i] = new Serialization();
                MPI_Alloc_mem(tSize * sizeof(struct Treenode), MPI_INFO_NULL, &(tree[i]->treeArray));
        }

        MPI_Win_attach(win, tree[rank]->treeArray, tSize*sizeof(struct Treenode));
        MPI_Get_address(tree[rank]->treeArray, &my_disp);
        MPI_Allgather(&my_disp, 1, MPI_AINT, target_disp, 1, MPI_AINT, MPI_COMM_WORLD);

        t1_start = MPI_Wtime();
        MPI_Win_fence(0, win);
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                MPI_Get(tree[i]->treeArray, tSize, treeNodeStruct, i, target_disp[i], tSize, treeNodeStruct, win);
            }
        }
        MPI_Win_fence(0, win);
        t1_end = MPI_Wtime();
        time1 += t1_end - t1_start;
    }

    MPI_Win_detach(win,tree[rank]->treeArray);
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Free_mem(tree[i]->treeArray);
            tree[i]->treeArray = NULL;
        }
        delete tree[i];
    }

    std::cout << time1  << std::endl;
    MPI_Win_free(&win);
	MPI_Finalize();
    return 0;
}
