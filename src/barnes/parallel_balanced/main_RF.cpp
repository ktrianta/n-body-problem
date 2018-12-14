#include <iostream>
#include <fstream>
#include <vector>
#include <papi.h>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"
#include "serialization.hpp"
#include "sort.hpp"
#include <unistd.h>
#include <mpi.h>
#include <math.h>
#include <chrono>

using time_point_t = std::chrono::high_resolution_clock::time_point;

using namespace std;

int main(int argc, char** argv) {

    // ----- MPI ----- //
    int size,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    Treenode defTree;
    MPI_Datatype treeNodeStruct;
    int blocklength[] = {11,1,1,2};
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

    // *** PAPI *** //
    int EventSet = PAPI_NULL;
    long_long values[1] = {(long_long) 0};

    int retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        printf("PAPI library init error!\n");
        exit(1);
    }

    if ( PAPI_create_eventset(&EventSet) != PAPI_OK){
        printf("PAPI library create error!\n");
        exit(1);
    }

    if ( PAPI_add_event(EventSet,PAPI_TOT_CYC ) != PAPI_OK){
        printf("PAPI add event error!\n");
        exit(1);
    }

    long long int papi_tot = 0;
    // *** PAPI *** //

    PAPI_start(EventSet);

    sim::Parameters params;
    readArgs(argc, argv, params);


    size_t N0;
    if (params.n % size != 0){
        N0 = params.n + size - params.n % size;
    } 
    else { 
        N0 = params.n;
    }
    const size_t N = N0;

    const sim::data_type theta = params.theta;
    int* local_N = new int[size];
    int* local_Nb = new int[size];
    int* offset = new int[size];

    sim::data_type (*r)[8];
    sim::data_type (*r_local)[8] = new sim::data_type[N][8];
    sim::data_type (*a_local)[3] = new sim::data_type[N][3];
    std::fill(&a_local[0][0], &a_local[0][0] + N*3, 0);

    if (rank == 0) {
        r = new sim::data_type[N][8];

        if (readDataFromFile(params.in_filename, params.n, r) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            return -1;
        }
        if(params.n %size != 0){
            for (int i = N - size + params.n % size; i < N; i++){
                r[i][0]=0;
                r[i][1]=r[i-1][1]+0.003;
                r[i][2]=r[i-1][2]+0.003;
                r[i][3]=r[i-1][3]+0.003;
                r[i][4]=0;
                r[i][5]=0;
                r[i][6]=0;
            }
        }

        params.out_filename = params.in_filename;


        p_sort(r, N, size, local_N);
        local_Nb[0] = local_N[0] * 8; offset[0] = 0;
        for (int i = 1; i < size; i++) {
            local_Nb[i] = local_N[i] * 8;
            offset[i] = offset[i-1] + local_Nb[i-1];
        }
    }
    
    MPI_Bcast(local_N, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(local_Nb, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(offset, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(r, local_Nb, offset, MPI_DOUBLE, &r_local[0][0], local_Nb[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);


    sim::data_type xc, yc, zc, h2, w2, t2;

    boxComputation(local_N[rank], r_local, xc, yc, zc, w2, h2, t2);
    Serialization *tree[size];
    tree[rank] = new Serialization(xc, yc, zc, w2, h2, t2);

    for (int j = 0; j < local_N[rank]; j++) {
        tree[rank]->insert(rank*local_N[rank]+j, r_local[j][1], r_local[j][2], r_local[j][3], r_local[j][0]);
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

    MPI_Win_fence(0,win);
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Get(tree[i]->treeArray, treeSize[i], treeNodeStruct, i, target_disp[i], treeSize[i], treeNodeStruct, win);
        }
    }
    MPI_Win_fence(0,win);

    for (int j = 0; j < local_N[rank]; j++) {
        r_local[j][7] = 0;
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < local_N[rank]; j++) {
            tree[i]->computeAcceleration(rank*local_N[rank]+j, &r_local[j][1], a_local[j], sim::g, theta, r_local[j][7]);
        }
    }


    const size_t Ntimesteps = params.s;
    const sim::data_type dt = params.dt;


    for (int t = 0; t < Ntimesteps; t++) {
        for (int j = 0; j < local_N[rank]; j++) {
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

        double work = 0;
        for (int j = 0; j < local_N[rank]; j++) {
            work += r_local[j][7];
        }

        MPI_Gatherv(r_local, local_Nb[rank], MPI_DOUBLE, r, local_Nb, offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            p_sort(r, N, size, local_N);
            local_Nb[0] = local_N[0] * 8; offset[0] = 0;
            for (int i = 1; i < size; i++) {
                local_Nb[i] = local_N[i] * 8;
                offset[i] = offset[i-1] + local_Nb[i-1];
            }
        }
        MPI_Bcast(local_N, size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(local_Nb, size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(offset, size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatterv(r, local_Nb, offset, MPI_DOUBLE, &r_local[0][0], local_Nb[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Win_detach(win,tree[rank]->treeArray);
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                MPI_Free_mem(tree[i]->treeArray);
                tree[i]->treeArray = NULL;
            }
            delete tree[i];
        }

        boxComputation(local_N[rank], r_local, xc, yc, zc, w2, h2, t2);
        tree[rank] = new Serialization(xc, yc, zc, w2, h2, t2);
        for(int j = 0; j < local_N[rank]; j++) {
            tree[rank]->insert(rank*local_N[rank]+j, r_local[j][1], r_local[j][2], r_local[j][3], r_local[j][0]);
        }

        tSize = tree[rank]->position+1;
        MPI_Allgather(&tSize, 1, MPI_INT, &treeSize, 1, MPI_INT, MPI_COMM_WORLD);
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                tree[i] = new Serialization();
                MPI_Alloc_mem(treeSize[i] * sizeof(struct Treenode), MPI_INFO_NULL, &(tree[i]->treeArray));
            }
        }

        MPI_Win_attach(win, tree[rank]->treeArray, treeSize[rank]*sizeof(struct Treenode));
        MPI_Get_address(tree[rank]->treeArray, &my_disp);
        MPI_Allgather(&my_disp, 1, MPI_AINT, target_disp, 1, MPI_AINT, MPI_COMM_WORLD);


        MPI_Win_fence(0, win);
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                MPI_Get(tree[i]->treeArray, treeSize[i], treeNodeStruct, i, target_disp[i], treeSize[i], treeNodeStruct, win);
            }
        }
        MPI_Win_fence(0, win);

        for (int j = 0; j < local_N[rank]; j++) {
            r_local[j][7] = 0;
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < local_N[rank]; j++) {
                tree[i]->computeAcceleration(rank*local_N[rank]+j, &r_local[j][1], a_local[j], sim::g, theta, r_local[j][7]);
            }
        }

        for (int j = 0; j < local_N[rank]; j++) {
            r_local[j][4] += 0.5 * a_local[j][0] * dt;
            r_local[j][5] += 0.5 * a_local[j][1] * dt;
            r_local[j][6] += 0.5 * a_local[j][2] * dt;
        }

    }

    PAPI_stop(EventSet, values);
    papi_tot += values[0];

    printf("%lld\n", papi_tot);

    MPI_Win_detach(win,tree[rank]->treeArray);
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Free_mem(tree[i]->treeArray);
            tree[i]->treeArray = NULL;
        }
        delete tree[i];
    }


    delete[] r_local;
    delete[] a_local;
    delete[] target_disp;
    MPI_Finalize();

    return 0;
}
