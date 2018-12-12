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

void* operator new(size_t size) {
    void *p;

    if ((p = _mm_malloc(size,64)) == 0) {
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

    sim::data_type* a_ix = new sim::data_type[8];
    sim::data_type* a_iy = new sim::data_type[8];
    sim::data_type* a_iz = new sim::data_type[8];

    __m512d g = _mm512_set1_pd(-sim::g);
    __m512d zeros = _mm512_set1_pd(0);

    for (size_t i = 0; i < N; i++) {
        __m512d aix = _mm512_setzero_pd();
        __m512d aiy = _mm512_setzero_pd();
        __m512d aiz = _mm512_setzero_pd();

        __m512d rix   = _mm512_set1_pd( rx[i] );
        __m512d riy   = _mm512_set1_pd( ry[i] );
        __m512d riz   = _mm512_set1_pd( rz[i] );

        for (int j = 0; j < N; j += 8) {
            __m512d rjx   = _mm512_load_pd( rx+j );
            __m512d rjy   = _mm512_load_pd( ry+j );
            __m512d rjz   = _mm512_load_pd( rz+j );
            __m512d masses = _mm512_load_pd(m+j);
            __m512d mg = _mm512_mul_pd(g,masses);

            __m512d subx  = _mm512_sub_pd( rjx, rix );
            __m512d suby  = _mm512_sub_pd( rjy, riy );
            __m512d subz  = _mm512_sub_pd( rjz, riz );

            __m512d multx = _mm512_mul_pd( subx, subx );
            __m512d multy = _mm512_mul_pd( suby, suby );
            __m512d multz = _mm512_mul_pd( subz, subz );

            __m512d add   = _mm512_add_pd(_mm512_add_pd(multx, multy), multz);
            __m512d denom = _mm512_mul_pd(_mm512_sqrt_pd(add), add);
            __mmask8 cmp_res = _mm512_cmp_pd_mask(zeros, denom, _CMP_EQ_OQ);
            __mmask8 notmask = _knot_mask8(cmp_res);
            //		__m512d a_i = _mm512_mul_pd(_mm512_mul_pd(_mm512_set1_pd(-sim::g),masses),_mm512_rcp14_pd(denom));
            //__m512d a_i = _mm512_div_pd(mg,denom);
            __m512d a_i = _mm512_mul_pd(mg, _mm512_rcp14_pd(denom));

            a_i = _mm512_maskz_andnot_pd(notmask, zeros, a_i);
            //		a_i = _mm512_andnot_si512(_mm512_castsi512_pd(cmp_res), a_i);

            aix = _mm512_sub_pd(aix,_mm512_mul_pd(a_i,subx));
            aiy = _mm512_sub_pd(aiy,_mm512_mul_pd(a_i,suby));
            aiz = _mm512_sub_pd(aiz,_mm512_mul_pd(a_i,subz));

            _mm512_store_pd(a_ix,aix);
            _mm512_store_pd(a_iy,aiy);
            _mm512_store_pd(a_iz,aiz);
        }

        ax[i] = a_ix[0] + a_ix[1] + a_ix[2] + a_ix[3] + a_ix[4] + a_ix[5] + a_ix[6] + a_ix[7];
        ay[i] = a_iy[0] + a_iy[1] + a_iy[2] + a_iy[3] + a_iy[4] + a_iy[5] + a_iy[6] + a_iy[7];
        az[i] = a_iz[0] + a_iz[1] + a_iz[2] + a_iz[3] + a_iz[4] + a_iz[5] + a_iz[6] + a_iz[7];
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
    writeDataToFile(params.n, rx, ry, rz, out_file);
    io_end = std::chrono::high_resolution_clock::now();
    io_time += std::chrono::duration< double >(io_end - io_start).count();

//Computation of Initial Energy
    double initialKEnergy = 0;
    double initialPEnergy = 0;
    double initialEnergy = 0;
    for (int i = 0; i < N; i++){
        initialKEnergy += m[i] * (u[i][0]*u[i][0] + u[i][1]*u[i][1] + u[i][2]*u[i][2])/2.;
        for (int j = 0; j < i; j++){
            double denominator = sqrt((rx[j]-rx[i])*(rx[j]-rx[i]) + (ry[j]-ry[i])*(ry[j]-ry[i]) +
                                      (rz[j]-rz[i])*(rz[j]-rz[i]));
            initialPEnergy -= sim::g*m[i]*m[j]/denominator;
            }
        }
        initialEnergy = initialKEnergy + initialPEnergy;

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

        io_start = std::chrono::high_resolution_clock::now();
        if (t % 200 == 0) {
            writeDataToFile(N, rx, ry, rz, out_file);
        }
        io_end = std::chrono::high_resolution_clock::now();
        io_time += std::chrono::duration< double >(io_end - io_start).count();
    }

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
