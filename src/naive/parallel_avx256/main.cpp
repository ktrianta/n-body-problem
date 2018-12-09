#include <mpi.h>
#include <immintrin.h>
#include <cmath>
#include <fstream>
#include <chrono>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"

using time_point_t = std::chrono::high_resolution_clock::time_point;
void* operator new(size_t size) {
    void *p;

    if ((p = _mm_malloc(size, 32)) == 0) {
        static const std::bad_alloc nomem;
        throw (nomem);
    }

    return (p);
}

void operator delete(void* ptr) noexcept {
    _mm_free(ptr);
}


void computeAcceleration(const size_t N,sim::data_type (*rx),sim::data_type (*ry),sim::data_type (*rz),sim::data_type (*ax),
	       				sim::data_type (*ay),sim::data_type (*az),sim::data_type (*m), const int local_N, const int offset) {
    
//	std::fill(&a[offset][0], &a[offset][0] + local_N*3, 0);

	std::fill(&ax[offset], &ax[offset] + local_N, 0);
	std::fill(&ay[offset], &ay[offset] + local_N, 0);
	std::fill(&az[offset], &az[offset] + local_N, 0);
	
	sim::data_type* a_ix = new sim::data_type[4];
	sim::data_type* a_iy = new sim::data_type[4];
	sim::data_type* a_iz = new sim::data_type[4];

       	for (size_t i = offset; i < offset + local_N; i++) {

		__m256d aix = _mm256_setzero_pd();
		__m256d aiy = _mm256_setzero_pd();
        	__m256d aiz = _mm256_setzero_pd();

		__m256 rix = _mm256_set1_pd( rx[i] );
		__m256 riy = _mm256_set1_pd( ry[i] );
		__m256 riz = _mm256_set1_pd( rz[i] );

	for (size_t j = 0; j < N; j+=4 ) {
            if (i == j) continue;
		
	    __m256d rjx = _mm256_load_pd(rx+j);
	    __m256d rjy = _mm256_load_pd(ry+j);
	    __m256d rjz = _mm256_load_pd(rz+j);
            
	    __m256d subx = _mm256_sub_pd(rjx,rix);
	    __m256d suby = _mm256_sub_pd(rjy,riy);
	    __m256d subz = _mm256_sub_pd(rjz,riz);
		
	    __m256d multx = _mm256_mul_pd(subx,subx);
	    __m256d multy = _mm256_mul_pd(suby,suby);
	    __m256d multz = _mm256_mul_pd(subz,subz);

	    __m256d add = _mm256_add_pd(_mm256_add_pd(multx,multy),multz);
	    __m256d denom = _mm256_mul_pd(_mm256_sqrt_pd(add),add);
	    __m256d zeros = _mm256_set1_pd(0);

	    __m256i cmp_res = _mm256_cmpeq_epi64(_mm256_castpd_si256(zeros), _mm256_castpd_si256(denom) );
	    __m256d masses = _mm256_load_pd(m+j);

	   __m256d a_i = _mm256_div_pd(_mm256_mul_pd(_mm256_set1_pd(-sim::g),masses),denom);
	   a_i = _mm256_andnot_pd(_mm256_castsi256_pd(cmp_res), a_i);

	   aix = _mm256_sub_pd(aix,_mm256_mul_pd(a_i,subx));
	   aiy = _mm256_sub_pd(aiy,_mm256_mul_pd(a_i,suby));
	   aiz = _mm256_sub_pd(aiz,_mm256_mul_pd(a_i,subz));

	  _mm256_store_pd(a_ix,aix);
	  _mm256_store_pd(a_iy,aiy);
	  _mm256_store_pd(a_iz,aiz);
	}
        ax[i] = a_ix[0] + a_ix[1] + a_ix[2] + a_ix[3]; 
        ay[i] = a_iy[0] + a_iy[1] + a_iy[2] + a_iy[3];
        az[i] = a_iz[0] + a_iz[1] + a_iz[2] + a_iz[3];
    }
	delete[] a_ix;
	delete[] a_iy;
	delete[] a_iz;
}

int main(int argc, char** argv) {

    // *** MPI *** // 
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    // *** MPI *** // 
    
    time_point_t prog_start;
    time_point_t prog_end;
    time_point_t io_start;
    time_point_t io_end;
    time_point_t comp_start;
    time_point_t comp_end;
    time_point_t comm_start;
    time_point_t comm_end;
    double  prog_time = 0;
    double io_time = 0;
    double comp_time = 0;
    double comm_time = 0;

    prog_start = std::chrono::high_resolution_clock::now();

    sim::Parameters params;
    readArgs(argc, argv, params);

    const size_t N = params.n;
    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*rx) = new sim::data_type[N];
    sim::data_type (*ry) = new sim::data_type[N];
    sim::data_type (*rz) = new sim::data_type[N];
    
    sim::data_type (*u)[3] = new sim::data_type[N][3];
   
    sim::data_type (*ax) = new sim::data_type[N];
    sim::data_type (*ay) = new sim::data_type[N];
    sim::data_type (*az) = new sim::data_type[N];
    std::fill(m, m+N, 1.0/N); 
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    
    
    std::fill(&ax[0], &ax[0] + N, 0);
    std::fill(&ay[0], &ax[0] + N, 0);
    std::fill(&az[0], &ay[0] + N, 0);

    // PROCESS 0 initialize position vector r.
    io_start = std::chrono::high_resolution_clock::now();
    std::ofstream out_file;
    if (rank == 0) {
        if (readDataFromFile(params.in_filename, N, m, rx, ry, rz, u) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            delete[] m;
            delete[] rx;
            delete[] ry;
            delete[] rz;
            delete[] u;
            delete[] ax;
            delete[] ay;
            delete[] az;
            return -1;
        }
        params.out_filename = params.in_filename;
        openFileToWrite(out_file, params.out_filename, params.out_dirname);
        writeDataToFile(N, rx, ry, rz, out_file);
    }
    io_end = std::chrono::high_resolution_clock::now();
    io_time += std::chrono::duration< double >(io_end - io_start).count();

// Computation of Initial Energy    
    double initialKEnergy = 0;
    double initialPEnergy = 0;
    double initialEnergy = 0;
    if (rank == 0){
        for (int i = 0; i < N; i++){
            initialKEnergy += m[i] * (u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])/2.;
            for (int j = 0; j < i; j++){
                double denominator = sqrt((rx[j]-rx[i])*(rx[j]-rx[i]) + (ry[j]-ry[i])*(ry[j]-ry[i]) +
                                          (rz[j]-rz[i])*(rz[j]-rz[i]));
                initialPEnergy -= sim::g*m[i]*m[j]/denominator;
                }
            }
            initialEnergy = initialKEnergy + initialPEnergy;
    }

    comm_start = std::chrono::high_resolution_clock::now();
    // SEND the position vector r from Process 0 to all processes.
    MPI_Bcast(&rx[0],N, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&ry[0],N, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&rz[0],N, MPI_DOUBLE,0, MPI_COMM_WORLD);
    
    MPI_Bcast(&u[0][0],N*3, MPI_DOUBLE,0, MPI_COMM_WORLD);
    MPI_Bcast(&m[0],N, MPI_DOUBLE,0, MPI_COMM_WORLD);
    comm_end = std::chrono::high_resolution_clock::now();
    comm_time += std::chrono::duration< double >(comm_end - comm_start).count();

//  std::ofstream out_file;
//  if (rank ==0) {
//      openFileToWrite(out_file, params.out_filename, params.out_dirname);
//      writeDataToFile(params.n, r, u, out_file);
//  }

    const sim::data_type dt = params.dt;
    const size_t timesteps = params.s;

    size_t local_N[size];
    int local_Nx3[size];
    size_t local_N_int = N/size;
    size_t rem = N - local_N_int * size;
    size_t counter = 0;

    for (size_t i = 0; i < size; i++) {
        local_N[i] = local_N_int;
        if (counter < rem) {
            local_N[i] += 1;
            counter ++;
        }
        local_Nx3[i] = local_N[i]*3;
    }
    
    sim::data_type (*rx_local) = new sim::data_type[local_N[rank]];
    sim::data_type (*ry_local) = new sim::data_type[local_N[rank]];
    sim::data_type (*rz_local) = new sim::data_type[local_N[rank]];

    size_t offset[size]; offset[0] = 0;
    int offset_x3[size]; offset_x3[0] = 0;

    for (size_t i = 1; i < size; i++) {
        offset[i] = offset[i-1] + local_N[i-1];
        offset_x3[i] = offset_x3[i-1] + local_Nx3[i-1];
    }

    for (size_t i = 0, j = offset[rank], end = local_N[rank]; i < end; i++, j++) {
        rx_local[i] = rx[j];
        ry_local[i] = ry[j];
        rz_local[i] = rz[j];
    }

    comp_start = std::chrono::high_resolution_clock::now();
    computeAcceleration(N, rx,ry,rz, ax,ay,az, m, local_N[rank], offset[rank]);
    comp_end = std::chrono::high_resolution_clock::now();
    comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

    //Start benchmark
//  double t1_comp, t2_comp, t_tot_comp, t1_comm, t2_comm, t_tot_comm;
    for (size_t t = 0; t < timesteps; t++) {
        for (size_t j = 0, idx = offset[rank]; j < local_N[rank]; j++, idx++) {
            u[idx][0] += 0.5 * ax[idx] * dt;
            u[idx][1] += 0.5 * ay[idx] * dt;
            u[idx][2] += 0.5 * az[idx] * dt;
            rx_local[j] += u[idx][0] * dt;
            ry_local[j] += u[idx][1] * dt;
            rz_local[j] += u[idx][2] * dt;
        }
//        t1_comm = MPI_Wtime();
        comm_start = std::chrono::high_resolution_clock::now();


	MPI_Allgatherv( &(rx_local[0]), local_N[rank], MPI_DOUBLE, &(rx[0]), local_Nx3, offset_x3, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv( &(ry_local[0]), local_N[rank], MPI_DOUBLE, &(ry[0]), local_Nx3, offset_x3, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv( &(rz_local[0]), local_N[rank], MPI_DOUBLE, &(rz[0]), local_Nx3, offset_x3, MPI_DOUBLE, MPI_COMM_WORLD);
        
	comm_end = std::chrono::high_resolution_clock::now();
        comm_time += std::chrono::duration< double >(comm_end - comm_start).count();
//        t2_comm  = MPI_Wtime();
//        t_tot_comm += (t2_comm-t1_comm);

//        t1_comp = MPI_Wtime();
        comp_start = std::chrono::high_resolution_clock::now();
        computeAcceleration(N, rx,ry,rz, ax,ay,az, m, local_N[rank], offset[rank]);
        comp_end = std::chrono::high_resolution_clock::now();
        comp_time += std::chrono::duration< double >(comp_end - comp_start).count();
//        t2_comp  = MPI_Wtime();
//        t_tot_comp += (t2_comp-t1_comp);

        for (size_t idx = offset[rank], end = offset[rank] + local_N[rank]; idx < end; idx++) {
            u[idx][0] += 0.5 * ax[idx] * dt;
            u[idx][1] += 0.5 * ay[idx] * dt;
            u[idx][2] += 0.5 * az[idx] * dt;
	    }

        io_start = std::chrono::high_resolution_clock::now();
        if (rank == 0) {
            if (t % 200 == 0){
                writeDataToFile(N, rx,ry,rz, out_file);
            }   
        }
        io_end = std::chrono::high_resolution_clock::now();
        io_time += std::chrono::duration< double >(io_end - io_start).count();
    }

//Computation of final energy and estimate error
    if (rank == 0){
        double energy =0;
        double kineticEnergy = 0;
        double potentialEnergy = 0;
        for (int i = 0; i < N; i++){
            kineticEnergy += m[i] * (u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])/2.;
            for (int j = 0; j < i; j++){
                double denominator = sqrt((rx[j]-rx[i])*(rx[j]-rx[i]) + (ry[j]-ry[i])*(ry[j]-ry[i]) +
                                          (rz[j]-rz[i])*(rz[j]-rz[i]));
                potentialEnergy -= sim::g*m[i]*m[j]/denominator;
                }
            }
            energy = kineticEnergy + potentialEnergy;
        std::cout << "initial energy is = " << initialEnergy << "Error in total energy at the end of simulation = " << (energy - initialEnergy)/initialEnergy*100 << "%" <<  std::endl; 
    }

    prog_end = std::chrono::high_resolution_clock::now();
    prog_time += std::chrono::duration< double >(prog_end - prog_start).count();
    double plotData_comp;
    double plotData_comm;
    double plotData_io;
    double plotData_prog;
    MPI_Reduce(&comp_time, &plotData_comp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time, &plotData_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&io_time, &plotData_io, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&prog_time, &plotData_prog, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if( rank== 0 ) {
        FILE *plotFile;
        plotFile = fopen("plotData.txt", "a");
        fprintf(plotFile, "%lf, %lf, %lf, %lf\n", plotData_prog, plotData_comp, plotData_io, plotData_comm);
        fclose(plotFile);
    }

    delete[] m;
    delete[] rx;
    delete[] ry;
    delete[] rz;
    delete[] u;
    delete[] ax;
    delete[] ay;
    delete[] az;

    MPI_Finalize();
    return 0;
}
