#include <math.h>
#include <fstream>
#include <iostream>
#include "io.hpp"
#include "args.hpp"
#include "types.hpp"
#include "energy.hpp"
#include "initialization.hpp"
#include <immintrin.h>
#include <chrono>

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

void computeAcceleration(const size_t N, sim::data_type (*rx), sim::data_type (*ry), sim::data_type (*rz),
        sim::data_type (*ax), sim::data_type (*ay), sim::data_type (*az), sim::data_type *m) {
    std::fill(&ax[0], &ax[0] + N, 0);
    std::fill(&ay[0], &ay[0] + N, 0);
    std::fill(&az[0], &az[0] + N, 0);

    sim::data_type* a_ix = new sim::data_type[4];
    sim::data_type* a_iy = new sim::data_type[4];
    sim::data_type* a_iz = new sim::data_type[4];

    for (size_t i = 0; i < N; i++) {
        __m256d aix = _mm256_setzero_pd();
        __m256d aiy = _mm256_setzero_pd();
        __m256d aiz = _mm256_setzero_pd();

        __m256d rix   = _mm256_set1_pd( rx[i] );
        __m256d riy   = _mm256_set1_pd( ry[i] );
        __m256d riz   = _mm256_set1_pd( rz[i] );

        for (int j = 0; j < N ; j += 4) {
            __m256d rjx   = _mm256_load_pd( rx+j );
            __m256d rjy   = _mm256_load_pd( ry+j );
            __m256d rjz   = _mm256_load_pd( rz+j );

            __m256d subx  = _mm256_sub_pd( rjx, rix );
            __m256d suby  = _mm256_sub_pd( rjy, riy );
            __m256d subz  = _mm256_sub_pd( rjz, riz );

            __m256d multx = _mm256_mul_pd( subx, subx );
            __m256d multy = _mm256_mul_pd( suby, suby );
            __m256d multz = _mm256_mul_pd( subz, subz );

            __m256d add   = _mm256_add_pd(_mm256_add_pd(multx, multy), multz);
            __m256d denom = _mm256_mul_pd(_mm256_sqrt_pd(add), add);
            __m256d zeros = _mm256_set1_pd(0);

            __m256i cmp_res = _mm256_cmpeq_epi64(_mm256_castpd_si256(zeros), _mm256_castpd_si256(denom));
            __m256d masses = _mm256_load_pd(m+j);

            __m256d a_i = _mm256_div_pd(_mm256_mul_pd(_mm256_set1_pd(-sim::g),masses),denom);

            //		__m256d a_i = _mm256_mul_pd( _mm256_mul_pd(_mm256_set1_pd(-sim::g),masses),
            //			     	_mm256_rcp_ps(denom) );

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
    time_point_t prog_start;
    time_point_t prog_end;
    time_point_t io_start;
    time_point_t io_end;
    time_point_t comp_start;
    time_point_t comp_end;
    double  prog_time = 0;
    double io_time = 0;
    double comp_time = 0;

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
    std::fill(&ay[0], &ay[0] + N, 0);
    std::fill(&az[0], &az[0] + N, 0);

    io_start = std::chrono::high_resolution_clock::now();
    if (readDataFromFile(params.in_filename, N, m, rx, ry, rz, u) == -1) {
        std::cerr << "File " << params.in_filename << " not found!" << std::endl;
        return -1;
    }
    params.out_filename = params.in_filename;
    std::ofstream out_file;
    openFileToWrite(out_file, params.out_filename, params.out_dirname);

    sim::data_type initialEnergy;
    if (params.en_comp == true)
       initialEnergy = energy(N, rx, ry, rz, u , m);

    comp_start = std::chrono::high_resolution_clock::now();
    computeAcceleration(params.n, rx, ry, rz, ax, ay, az,  m);
    comp_end = std::chrono::high_resolution_clock::now();
    comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

    const size_t Ntimesteps = params.s;
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

        comp_start = std::chrono::high_resolution_clock::now();
        computeAcceleration(N, rx, ry, rz, ax, ay, az,  m);
        comp_end = std::chrono::high_resolution_clock::now();
        comp_time += std::chrono::duration< double >(comp_end - comp_start).count();

        for (size_t j = 0; j < N; j++) {
            u[j][0] += 0.5 * ax[j] * dt;
            u[j][1] += 0.5 * ay[j] * dt;
            u[j][2] += 0.5 * az[j] * dt;
        }

        if (params.wr_data == true)
        {
            io_start = std::chrono::high_resolution_clock::now();
            if (t % 200 == 0) {
                writeDataToFile(N, rx, ry, rz, out_file);
            }
            io_end = std::chrono::high_resolution_clock::now();
            io_time += std::chrono::duration< double >(io_end - io_start).count();
        }
    }

    if (params.en_comp == true)
    {
        sim::data_type finalEnergy = energy(N, rx, ry, rz, u, m);
        printEnergy(finalEnergy, initialEnergy);
    }

    prog_end = std::chrono::high_resolution_clock::now();
    prog_time += std::chrono::duration< double >(prog_end - prog_start).count();
    FILE *plotFile;
    plotFile = fopen("plotData.txt", "a");
    fprintf(plotFile, "%lf,  %lf, %lf \n", prog_time, comp_time, io_time);
    fclose(plotFile);

    delete[] m;
    delete[] rx;
    delete[] ry;
    delete[] rz;
    delete[] u;
    delete[] ax;
    delete[] ay;
    delete[] az;

    return 0;
}
