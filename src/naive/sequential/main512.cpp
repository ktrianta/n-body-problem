#include <math.h>
#include <fstream>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include <immintrin.h>
#include <chrono>

using time_point_t = std::chrono::high_resolution_clock::time_point;
 
void * operator new(size_t size) throw(std::bad_alloc){
   void *p;
   if ((p = _mm_malloc(size,64)) == 0){
        static const std::bad_alloc nomem;
        throw (nomem);
   }
     return (p);
}
 void operator delete(void* ptr) noexcept{
	_mm_free(ptr);
}
 void computeAcceleration(const size_t N, sim::data_type (*rx), sim::data_type (*ry), sim::data_type (*rz), 
		 			sim::data_type (*ax), sim::data_type (*ay), sim::data_type (*az), sim::data_type *m) {
    
     std::fill(&ax[0], &ax[0] + N, 0);
     std::fill(&ay[0], &ay[0] + N, 0);
     std::fill(&az[0], &az[0] + N, 0);


     sim::data_type* a_ix = new sim::data_type[8]; 
     sim::data_type* a_iy = new sim::data_type[8];
     sim::data_type* a_iz = new sim::data_type[8];    
     
     for (size_t i = 0; i < N; i++) {
	__m512d aix = _mm512_setzero_pd();
	__m512d aiy = _mm512_setzero_pd();
	__m512d aiz = _mm512_setzero_pd();
 	
	__m512d rix   = _mm512_set1_pd( rx[i] );
	__m512d riy   = _mm512_set1_pd( ry[i] );
	__m512d riz   = _mm512_set1_pd( rz[i] );
 		
	for (int j = 0; j<N ; j+= 8){
		__512d rjx   = _mm512_load_pd( rx+j );
		__512d rjy   = _mm512_load_pd( ry+j );
		__512d rjz   = _mm512_load_pd( rz+j );

 		__512d subx  = _mm512_sub_pd( rjx, rix );
		__512d suby  = _mm512_sub_pd( rjy, riy );
		__512d subz  = _mm512_sub_pd( rjz, riz );
 		
		__512d multx = _mm512_mul_pd( subx, subx );
		__512d multy = _mm512_mul_pd( suby, suby );
		__512d multz = _mm512_mul_pd( subz, subz );  
 	
		__512d add      = _mm512_add_pd(_mm512_add_pd(multx, multy), multz); 
 		__512d denom    = _mm512_mul_pd(_mm512_sqrt_pd(add), add);        
		__512d zeros = _mm512_set1_pd(0);
		
		__m512i cmp_res = _mm512_cmpeq_epi64(_mm512_castpd_si512(zeros),_mm512_castpd_si512(denom));	
		__m512d masses = _mm512_load_pd(m+j);
 		
		__m512d a_i = _mm512_div_pd(_mm512_mul_pd(_mm512_set1_pd(-sim::g),masses),denom);	
//
//		__m512d a_i = _mm256_mul_pd( _mm256_mul_pd(_mm256_set1_pd(-sim::g),masses),
//			     	_mm512_rcp14_pd(denom) );	


		a_i = _mm512_andnot_si512(_mm512_castsi512_pd(cmp_res), a_i);
 	
		aix = _mm512_sub_pd(aix,_mm512_mul_pd(a_i,subx));
		aiy = _mm512_sub_pd(aiy,_mm512_mul_pd(a_i,suby));
		aiz = _mm512_sub_pd(aiz,_mm512_mul_pd(a_i,subz));
		
		
		_mm512_store_pd(a_ix,aix);
		_mm512_store_pd(a_iy,aiy);
		_mm512_store_pd(a_iz,aiz);
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
    std::fill(&ay[0], &ay[0] + N, 0);
    std::fill(&az[0], &az[0] + N, 0);

     if (!params.in_filename.empty()) {
        if (readDataFromFile(params.in_filename, N, m, rx, ry, rz, u) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            return -1;
        }
        params.out_filename = params.in_filename;
     }

    std::ofstream out_file;
    openFileToWrite(out_file, params.out_filename, params.out_dirname);
    writeDataToFile(params.n, rx, ry, rz, u, out_file);
    
    time_point_t tstart_;
    time_point_t tend_;
    double T = 0;

    tstart_ = std::chrono::high_resolution_clock::now();
    computeAcceleration(params.n, rx, ry, rz, ax, ay, az,  m);
    tend_ = std::chrono::high_resolution_clock::now();
    T += std::chrono::duration< double >(tend_ - tstart_).count();


    const size_t Ntimesteps = params.t / params.dt + 1;
    const sim::data_type dt = params.dt;
     for (size_t t = 0; t < Ntimesteps; t++) {
        for (size_t j = 0; j < N; j++) {
            u[j][0] += 0.5 * ax[j] * dt;
            u[j][1] += 0.5 * ay[j] * dt;
            u[j][2] += 0.5 * az[j] * dt;
            rx[j] += u[j][0] * dt;
            ry[j] += u[j][1] * dt;
            rz[j] += u[j][2] * dt;
        }

      	tstart_ = std::chrono::high_resolution_clock::now();
        computeAcceleration(N, rx, ry, rz, ax, ay, az,  m);
        tend_ = std::chrono::high_resolution_clock::now();
    	T += std::chrono::duration< double >(tend_ - tstart_).count();

         for (size_t j = 0; j < N; j++) {
            u[j][0] += 0.5 * ax[j] * dt;
            u[j][1] += 0.5 * ay[j] * dt;
            u[j][2] += 0.5 * az[j] * dt;
        }
         if (t % 200 == 0) {
            writeDataToFile(N, rx, ry, rz, u, out_file);
        }
    }


    std::cout << T << std::endl;


        return 0;
}
