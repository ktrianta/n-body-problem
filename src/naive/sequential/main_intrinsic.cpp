#include <math.h>
#include <fstream>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include <immintrin.h>

void * operator new(size_t size) throw(std::bad_alloc){
   void *p;
   if ((p = _mm_malloc(size,16)) == 0){
        static const std::bad_alloc nomem;
        throw (nomem);
   }
     return (p);
}

void operator delete(void* ptr) noexcept{
	_mm_free(ptr);
}


void computeAcceleration(const size_t N, sim::data_type (*r)[3], sim::data_type (*a)[3], sim::data_type *m) {
    std::fill(&a[0][0], &a[0][0] + N*3, 0);

    for (size_t i = 0; i < N; i++) {
        sim::data_type a_i0 = 0;  // accumulate accelaration values for particle i and
        sim::data_type a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        sim::data_type a_i2 = 0;

	__m256 ai0 = _mm256_setzero_pd();
	__m256 ai1 = _mm256_setzero_pd();
	__m256 ai2 = _mm256_setzero_pd();

	 sim::data_type a_data[12];

   sim::data_type r_data_j[12];
   sim::data_type r_data_i[12];

        for (size_t j = 0; j < N; j+=4) { /* Get contribution of all j particules on particule i. */

		for (int k=0; k<4; k++){
			/* Particule i loaded */
			r_data_i[k]   = r[i][0];
			r_data_i[k+4] = r[i][1];
			r_data_i[k+8] = r[i][2];

			/* Load other particules j */
			r_data_j[k]   = r[j+k][0];
			r_data_j[k+4] = r[j+k][1];
			r_data_j[k+8] = r[j+k][2];
		}

		__m256 rj0   = _mm256_load_pd( r_data_j	    );
		__m256 rj1   = _mm256_load_pd( &r_data_j[4] );
		__m256 rj2   = _mm256_load_pd( &r_data_j[8] );

		__m256 ri0   = _mm256_load_pd( r_data_i     );
		__m256 ri1   = _mm256_load_pd( &r_data_i[4] );
		__m256 ri2   = _mm256_load_pd( &r_data_i[8] );

		__m256 sub0  = _mm256_sub_pd( rj0, ri0 );
		__m256 sub1  = _mm256_sub_pd( rj1, ri1 );
		__m256 sub2  = _mm256_sub_pd( rj2, ri2 );

		__m256 mult0 = _mm256_mul_pd( sub0, sub0 );	// delta_x squared
		__m256 mult1 = _mm256_mul_pd( sub1, sub1 );	// delta_y squared
		__m256 mult2 = _mm256_mul_pd( sub2, sub2 );     // delta_z squared

		__m256 add      = _mm256_add_pd( _mm256_add_pd(mult0, mult1), mult2  ); 	// delta_x squared + delta_y squared + delta_z squared

		__m256 denom    = _mm256_mul_pd(_mm256_sqrt_pd(add), add);         // sqrt(expr) * expr , expr = above.
		__m256 zeros = _mm256_set1_pd(0);
		__m256 cmp_res = _mm256_cmpeq_epi64(zeros,denom);
	
		__m256 masses = _mm256_load_pd(m+j);
		__m256 grav = _mm256_set1_pd(-sim::g);

		__m256 a_i = _mm256_div_pd(_mm256_mul_pd(grav,masses),denom);
		a_i = _mm256_andnot_si256(cmp_res, a_i);

		ai0 = _mm256_sub_pd(ai0,_mm256_mul_pd(a_i,sub0));
		ai1 = _mm256_sub_pd(ai1,_mm256_mul_pd(a_i,sub1));
		ai2 = _mm256_sub_pd(ai2,_mm256_mul_pd(a_i,sub2));

		_mm256_store_pd(a_data,ai0);
		_mm256_store_pd(&a_data[4],ai1);
		_mm256_store_pd(&a_data[8],ai2);

        }
		for (int k=0; k<4; k++){
			a[i][0] += a_data[k];
			a[i][1] += a_data[k+4];
			a[i][2] += a_data[k+8];
		}


    }
}

int main(int argc, char** argv) {
    sim::Parameters params;
    readArgs(argc, argv, params);

    const size_t N = params.n;
    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);


    if (!params.in_filename.empty()) {
        if (readDataFromFile(params.in_filename, N, m, r, u) == -1) {
            std::cerr << "File " << params.in_filename << " not found!" << std::endl;
            return -1;
        }
        params.out_filename = params.in_filename;
    } else {
        initializePositionOnSphere(N, r);
    }

    std::ofstream out_file;
    openFileToWrite(out_file, params.out_filename, params.out_dirname);
    writeDataToFile(params.n, r, u, out_file);

    computeAcceleration(params.n, r, a, m);

    const size_t Ntimesteps = params.t / params.dt + 1;
    const sim::data_type dt = params.dt;

    for (size_t t = 0; t < Ntimesteps; t++) {
        for (size_t j = 0; j < N; j++) {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
            r[j][0] += u[j][0] * dt;
            r[j][1] += u[j][1] * dt;
            r[j][2] += u[j][2] * dt;
        }

        computeAcceleration(N, r, a, m);

        for (size_t j = 0; j < N; j++) {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
        }

        if (t % 200 == 0) {
            writeDataToFile(N, r, u, out_file);
        }
    }
    return 0;
}
