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

    time_point_t prog_start;
    time_point_t prog_end;
    time_point_t io_start;
    time_point_t io_end;
    time_point_t comp_start;
    time_point_t comp_end;
    time_point_t tree_start;
    time_point_t tree_end;
    time_point_t comm_start;
    time_point_t comm_end;
    time_point_t gath_start;
    time_point_t gath_end;
    time_point_t aloc_start;
    time_point_t aloc_end;
    double  prog_time = 0;
    double io_time = 0;
    double comp_time = 0;
    double tree_time = 0;
    double comm_time = 0;
    double gath_time = 0;
    double aloc_time = 0;

    prog_start = std::chrono::high_resolution_clock::now();

  	MPI_Win win;
  	MPI_Win_create_dynamic(MPI_INFO_NULL, MPI_COMM_WORLD, &win);

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

    io_start = std::chrono::high_resolution_clock::now();
    if (rank == 0) {
        r = new sim::data_type[N][8];

        if (readDataFromFile(params.in_filename, N - size + params.n % size, r) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            return -1;
        }
        for (int i = N - size + params.n % size; i < N; i++){
            r[i][0]=0;
            r[i][1]=r[i-1][1]+0.003;
            r[i][2]=r[i-1][2]+0.003;
            r[i][3]=r[i-1][3]+0.003;
            r[i][4]=0;
            r[i][5]=0;
            r[i][6]=0;
        }

        params.out_filename = params.in_filename;

        io_end = std::chrono::high_resolution_clock::now();
        io_time += std::chrono::duration< double >(io_end - io_start).count();

        p_sort(r, N, size, local_N);
        local_Nb[0] = local_N[0] * 8; offset[0] = 0;
        for (int i = 1; i < size; i++) {
            local_Nb[i] = local_N[i] * 8;
            offset[i] = offset[i-1] + local_Nb[i];
        }
    }
    MPI_Bcast(local_N, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(local_Nb, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(offset, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(r, local_Nb, offset, MPI_DOUBLE, &r_local[0][0], local_Nb[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    io_start = std::chrono::high_resolution_clock::now();
    std::ofstream out_file;
    if (rank == 0) {
        openFileToWrite(out_file, params.out_filename, params.out_dirname);
        writeDataToFile(N, r, out_file);
    }
    io_end = std::chrono::high_resolution_clock::now();
    io_time += std::chrono::duration< double >(io_end - io_start).count();

    sim::data_type xc, yc, zc, h2, w2, t2;
    boxComputation(local_N[rank], r_local, xc, yc, zc, w2, h2, t2);

    tree_start = std::chrono::high_resolution_clock::now();
    Serialization *tree[size];
    tree[rank] = new Serialization(xc, yc, zc, w2, h2, t2);

    for (int j = 0; j < local_N[rank]; j++) {
        tree[rank]->insert(rank*local_N[rank]+j, r_local[j][1], r_local[j][2], r_local[j][3], r_local[j][0]);
    }
    tree_end = std::chrono::high_resolution_clock::now();
    tree_time += std::chrono::duration< double >(tree_end - tree_start).count();

    aloc_start = std::chrono::high_resolution_clock::now();
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
    aloc_end = std::chrono::high_resolution_clock::now();
    aloc_time += std::chrono::duration< double >(aloc_end - aloc_start).count();

    comm_start = std::chrono::high_resolution_clock::now();
    MPI_Win_fence(0,win);
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Get(tree[i]->treeArray, treeSize[i], treeNodeStruct, i, target_disp[i], treeSize[i], treeNodeStruct, win);
        }
    }
    MPI_Win_fence(0,win);
    comm_end = std::chrono::high_resolution_clock::now();
    comm_time += std::chrono::duration< double >(comm_end - comm_start).count();

    comp_start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < local_N[rank]; j++) {
        r_local[j][7] = 0;
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < local_N[rank]; j++) {
            tree[i]->computeAcceleration(rank*local_N[rank]+j, &r_local[j][1], a_local[j], sim::g, theta, r_local[j][7]);
        }
    }
    comp_end = std::chrono::high_resolution_clock::now();
    comp_time += std::chrono::duration< double >(comp_end - comp_start).count();


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

        std::cout << rank << " work " << work << std::endl;
        MPI_Gatherv(r_local, local_Nb[rank], MPI_DOUBLE, r, local_Nb, offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            p_sort(r, N, size, local_N);
            local_Nb[0] = local_N[0] * 8; offset[0] = 0;
            for (int i = 1; i < size; i++) {
                local_Nb[i] = local_N[i] * 8;
                offset[i] = offset[i-1] + local_Nb[i];
            }
        }
        MPI_Bcast(local_N, size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(local_Nb, size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(offset, size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatterv(r, local_Nb, offset, MPI_DOUBLE, &r_local[0][0], local_Nb[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        aloc_start = std::chrono::high_resolution_clock::now();
        MPI_Win_detach(win,tree[rank]->treeArray);
        for (int i = 0; i < size; i++) {
            if (i != rank) {
                MPI_Free_mem(tree[i]->treeArray);
                tree[i]->treeArray = NULL;
            }
            delete tree[i];
        }
        aloc_end = std::chrono::high_resolution_clock::now();
        aloc_time += std::chrono::duration< double >(aloc_end - aloc_start).count();

        boxComputation(local_N[rank], r_local, xc, yc, zc, w2, h2, t2);
        tree[rank] = new Serialization(xc, yc, zc, w2, h2, t2);
        for(int j = 0; j < local_N[rank]; j++) {
            tree[rank]->insert(rank*local_N[rank]+j, r_local[j][1], r_local[j][2], r_local[j][3], r_local[j][0]);
        }

        aloc_start = std::chrono::high_resolution_clock::now();
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
        aloc_end = std::chrono::high_resolution_clock::now();
        aloc_time += std::chrono::duration< double >(aloc_end - aloc_start).count();


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
        comp_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < local_N[rank]; j++) {
                tree[i]->computeAcceleration(rank*local_N[rank]+j, &r_local[j][1], a_local[j], sim::g, theta, r_local[j][7]);
            }
        }
        comp_end = std::chrono::high_resolution_clock::now();
        comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

        for (int j = 0; j < local_N[rank]; j++) {
            r_local[j][4] += 0.5 * a_local[j][0] * dt;
            r_local[j][5] += 0.5 * a_local[j][1] * dt;
            r_local[j][6] += 0.5 * a_local[j][2] * dt;
        }
    }

  //if (rank == 0){
  //    double energy =0;
  //    double kineticEnergy = 0;
  //    double potentialEnergy = 0;
  //    for (int i = 0; i < N; i++){
  //        kineticEnergy += r[i][0] * (r[i][4]*r[i][4] + r[i][5]*r[i][5] + r[i][6]*r[i][6])/2.;
  //        for (int j = 0; j < i; j++){
  //            double denominator = sqrt((r[j][1]-r[i][1])*(r[j][1]-r[i][1]) + (r[j][2]-r[i][2])*(r[j][2]-r[i][2]) +
  //                                      (r[j][3]-r[i][3])*(r[j][3]-r[i][3]));
  //            potentialEnergy -= sim::g*r[i][0]*r[j][0]/denominator;
  //            }
  //        }
  //        energy = kineticEnergy + potentialEnergy;
  //    std::cout << "initial energy is = " << initialEnergy << "Error in total energy at the end of simulation = " << (energy - initialEnergy)/initialEnergy*100 << "%" <<  endl; 
  //}


    MPI_Win_detach(win,tree[rank]->treeArray);
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            MPI_Free_mem(tree[i]->treeArray);
            tree[i]->treeArray = NULL;
        }
        delete tree[i];
    }


    prog_end = std::chrono::high_resolution_clock::now();
    prog_time += std::chrono::duration< double >(prog_end - prog_start).count();
//  printf("%lf,  %lf, %lf, %lf, %lf \n", prog_time, comp_time, io_time, tree_time, comm_time);
    double plotData_comp;
    double plotData_comm;
    double plotData_io;
    double plotData_prog;
    double plotData_tree;
    double plotData_gath;
    double plotData_aloc;
    MPI_Reduce(&comp_time, &plotData_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time, &plotData_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&io_time, &plotData_io, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&prog_time, &plotData_prog, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tree_time, &plotData_tree, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&gath_time, &plotData_gath, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&aloc_time, &plotData_aloc, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Win_free(&win);
	MPI_Finalize();
    if (rank == 0) {
        delete[] r;
        FILE *plotFile;
        plotFile = fopen("plotData.txt", "a");
        fprintf(plotFile, "%lf,  %lf, %lf, %lf, %lf, %lf, %lf \n", plotData_prog, plotData_comp, plotData_io, plotData_tree, plotData_comm, plotData_gath, plotData_aloc);
        fclose(plotFile);
    }
    delete[] r_local;
    delete[] a_local;
    delete[] target_disp;

    return 0;
}
