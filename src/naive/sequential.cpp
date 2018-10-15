#include <math.h> //for sqrt function
#include <stdio.h> //for printf
#include <stdlib.h> //for srand/rand prng
#include <time.h> //for seeding prng
#include <iostream>
#include <fstream>
#include <vector>
#include "array2d.hpp"
#include "initialization.hpp"

using namespace std;

const int N = 5;      // the number of particles
const float g = 1 ;     // gravitational constant
const float epsilon = 0.001;
const float epsilon2 = epsilon * epsilon;

void computeAcceleration(Array2D<float>& r, Array2D<float>& a, vector<float>& m)
{
    for (int i = 0; i < N; i++)
    {
        float a_i0 = 0;  // accumulate accelaration values for particle i and
        float a_i1 = 0;  // store them at the end of the loop iteration in a(i,x)
        for (int j = i+1; j < N; j++)
        {
            float rji[2];
            rji[0] = r(j, 0) - r(i, 0);
            rji[1] = r(j, 1) - r(i, 1);
            float r2 = rji[0] * rji[0] + rji[1] * rji[1];
            float denom = (r2+epsilon2) * sqrt(r2+epsilon2);
            float fi = -g * m[j] / denom;
            float fj = -g * m[i] / denom;
            a(j, 0) += fi * rji[0];
            a(j, 1) += fi * rji[1];
            a_i0 -= fj * rji[0];
            a_i1 -= fj * rji[1];
        }
        a(i, 0) += a_i0;  // a(i, 0) and a(i, 1) are accessed once here, avoiding
        a(i, 1) += a_i1;  // repeated accesses in the inner loop of j
    }
}

void writeDataToFile(Array2D<float>& r, Array2D<float>& u, ofstream& file)
{
    for (int i = 0; i < N; i++)
    {
        file << r(i, 0) << "   "
             << r(i, 1) << "   "
             << u(i, 0) << "   "
             << u(i, 1) << "\n";
    }
}

int main(int argc,char** argv)
{
    Array2D<float> r(N, 2);
    Array2D<float> u(N, 2, 0);
    Array2D<float> a(N, 2, 0);
    vector<float> m(N, 1.0/N);
    ofstream file;
    file.open("output.dat");

    initializePositionOnSphere(N, r);
    writeDataToFile(r, u, file);

    const float T = 10;
    const float dt = 0.00001;
    const int Ntimesteps = T/dt + 1;

    computeAcceleration(r, a, m);

    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = 0; j < N; j++)
        {
            u(j, 0) += 0.5 * a(j, 0) * dt;
            u(j, 1) += 0.5 * a(j, 1) * dt;
            r(j, 0) += u(j, 0) * dt;
            r(j, 1) += u(j, 1) * dt;
        }
        for (int j = 0; j < N; j++)
        {
            a(j, 0) = 0;
            a(j, 1) = 0;
        }

        computeAcceleration(r, a, m);

        for (int j = 0; j < N; j++)
        {
            u(j, 0) += 0.5 * a(j, 0) * dt;
            u(j, 1) += 0.5 * a(j, 1) * dt;
        }

        if (t % 200 == 0)
        {
            writeDataToFile(r, u, file);
        }
    }
    return 0;
}
