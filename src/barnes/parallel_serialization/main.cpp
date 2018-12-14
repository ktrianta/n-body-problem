#include <iostream>
#include <fstream>
#include <vector>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include "boxComputation.hpp"
#include "serialization.hpp"
#include "energy.hpp"
#include "sort.hpp"
#include <unistd.h>
#include <mpi.h>
#include <math.h>
#include <chrono>

using time_point_t = std::chrono::high_resolution_clock::time_point;

using namespace std;

int main(int argc, char** argv)
{
    // ----- MPI ----- //
    int size,rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

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
    double  prog_time = 0;
    double io_time = 0;
    double comp_time = 0;
    double tree_time = 0;
    double comm_time = 0;

    prog_start = std::chrono::high_resolution_clock::now();

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
    int local_N = N/size;

    sim::data_type (*r)[7] = new sim::data_type[N][7];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(&a[0][0], &a[0][0] + N*3, 0);

    sim::data_type initialEnergy = 0;

    if (rank == 0) {
        io_start = std::chrono::high_resolution_clock::now();
        r = new sim::data_type[N][7];

        if (readDataFromFile(params.in_filename, params.n , r) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            return -1;
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
        io_end = std::chrono::high_resolution_clock::now();
        io_time += std::chrono::duration< double >(io_end - io_start).count();
        p_sort(r, N, size);
        }

        params.out_filename = params.in_filename;
        
//      initialEnergy = energy(N, r);

    }
    // SEND the position vector r from Process 0 to all processes.
//    comm_start = std::chrono::high_resolution_clock::now();
    MPI_Bcast(&r[0][0],N*7, MPI_DOUBLE,0, MPI_COMM_WORLD);
//    comm_end = std::chrono::high_resolution_clock::now();
//    comm_time += std::chrono::duration< double >(comm_end - comm_start).count();

//  the center of the parent node and the half width and height
    sim::data_type xc, yc, zc, h2, w2, t2;



//  io_start = std::chrono::high_resolution_clock::now();
//  std::ofstream out_file;
//  if (rank == 0) {
//      openFileToWrite(out_file, params.out_filename, params.out_dirname);
//      writeDataToFile(N, r, out_file);
//  }
//  io_end = std::chrono::high_resolution_clock::now();
//  io_time += std::chrono::duration< double >(io_end - io_start).count();
   
    sim::data_type (*a_local)[3] = new sim::data_type[N][3];
    std::fill(&a_local[0][0], &a_local[0][0] + N*3, 0);

    
 
    tree_start = std::chrono::high_resolution_clock::now();
    boxComputation(local_N, rank*local_N, (rank+1)*local_N, r, xc, yc, zc, w2, h2, t2);
    Serialization *tree = new Serialization(xc, yc, zc, w2, h2, t2);
    for(int j = rank*local_N; j < (rank+1)*local_N; j++)
    {
        tree->insert(j, r[j][1], r[j][2], r[j][3], r[j][0]);
    }
    tree_end = std::chrono::high_resolution_clock::now();
    tree_time += std::chrono::duration< double >(tree_end - tree_start).count();

    const size_t Ntimesteps = params.s;
    const sim::data_type dt = params.dt;
    

    comp_start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < N; j++)
    {
        tree->computeAcceleration(j, &r[j][1], a_local[j], sim::g, theta);
    }
    comp_end = std::chrono::high_resolution_clock::now();
    comp_time += std::chrono::duration< double >(comp_end - comp_start).count();
   
 
    comm_start = std::chrono::high_resolution_clock::now();
    //    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(a_local,a,3*N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //    MPI_Barrier(MPI_COMM_WORLD);
    comm_end = std::chrono::high_resolution_clock::now();
    comm_time += std::chrono::duration< double >(comm_end - comm_start).count();



    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = 0; j < N; j++)
        {
            r[j][4] += 0.5 * a[j][0] * dt;
            r[j][5] += 0.5 * a[j][1] * dt;
            r[j][6] += 0.5 * a[j][2] * dt;
            r[j][1] += r[j][4] * dt;
            r[j][2] += r[j][5] * dt;
            r[j][3] += r[j][6] * dt;

            a_local[j][0] = 0;
            a_local[j][1] = 0;
            a_local[j][2] = 0;
            a[j][0] = 0;
            a[j][1] = 0;
            a[j][2] = 0;
        }
        p_sort(r, N, size);

        tree_start = std::chrono::high_resolution_clock::now();
        boxComputation(local_N, rank*local_N, (rank+1)*local_N, r, xc, yc, zc, w2, h2, t2);
        delete tree;
        tree= new Serialization(xc, yc, zc, w2, h2, t2);
        for(int j = rank*local_N; j < (rank+1)*local_N; j++)
        {
            tree->insert(j, r[j][1], r[j][2], r[j][3], r[j][0]);
        }
        tree_end = std::chrono::high_resolution_clock::now();
        tree_time += std::chrono::duration< double >(tree_end - tree_start).count();

        comp_start = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < N; j++)
        {
            tree->computeAcceleration(j, &r[j][1], a_local[j], sim::g, theta);
        }
        comp_end = std::chrono::high_resolution_clock::now();
        comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

      //  MPI_Barrier(MPI_COMM_WORLD);
        comm_start = std::chrono::high_resolution_clock::now();
        MPI_Allreduce(a_local,a,3*N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        comm_end = std::chrono::high_resolution_clock::now();
        comm_time += std::chrono::duration< double >(comm_end - comm_start).count();
     //   MPI_Barrier(MPI_COMM_WORLD);
        for (int j = 0; j < N; j++)
        {

            r[j][4] += 0.5 * a[j][0] * dt;
            r[j][5] += 0.5 * a[j][1] * dt;
            r[j][6] += 0.5 * a[j][2] * dt;
        }


//      io_start = std::chrono::high_resolution_clock::now();
//      if (t % 200 == 0 && rank == 0)
//      {
//          writeDataToFile(N, r, out_file);
//      }
//      io_end = std::chrono::high_resolution_clock::now();
//      io_time += std::chrono::duration< double >(io_end - io_start).count();

    }
//  if (rank == 0){
//      sim::data_type finalEnergy;
//      finalEnergy = energy(N, r);
//      printEnergy(finalEnergy, initialEnergy);
//  }

    prog_end = std::chrono::high_resolution_clock::now();
    prog_time += std::chrono::duration< double >(prog_end - prog_start).count();
    double plotData_comp;
    double plotData_comm;
    double plotData_io;
    double plotData_prog;
    double plotData_tree;
    MPI_Reduce(&comp_time, &plotData_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time, &plotData_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&io_time, &plotData_io, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&prog_time, &plotData_prog, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tree_time, &plotData_tree, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        delete[] r;
        FILE *plotFile;
        plotFile = fopen("plotData.txt", "a");
        fprintf(plotFile, "%lf,  %lf, %lf, %lf, %lf \n", plotData_prog, plotData_comp, plotData_io, plotData_tree, plotData_comm);
        fclose(plotFile);
    }

	MPI_Finalize();
    return 0;
}
