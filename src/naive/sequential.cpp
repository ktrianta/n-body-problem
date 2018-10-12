#include <math.h> //for sqrt function
#include <stdio.h> //for printf
#include <stdlib.h> //for srand/rand prng
#include <time.h> //for seeding prng
#include <iostream>
#include <fstream>
#include "initialization.h"

using namespace std;

const int N = 5;      // the number of particles
const float g = 1 ;     // gravitational constant
const float epsilon = 0.001;
const float epsilon2 = epsilon * epsilon;

void ComputeF(int j, float (*r)[2], float *m, float (*a)[2])
{
    for (int i = j+1; i < N; i++)
    {
        float rij[2];
        rij[0] = r[i][0] - r[j][0];
        rij[1] = r[i][1] - r[j][1];
        float r2 = rij[0] * rij[0] + rij[1] * rij[1];
        float denom = (r2+epsilon2) * sqrt(r2+epsilon2);
        float fi = -g * m[j] / denom;
        float fj = -g * m[i] / denom;
        a[i][0] += fi * rij[0];
        a[i][1] += fi * rij[1];
        a[j][0] -= fj * rij[0];
        a[j][1] -= fj * rij[1];
    }
}

int main(int argc,char** argv)
{
    float r[N][2], u[N][2], a[N][2], m[N];
    ofstream myfile;
    myfile.open ("output.dat");

    initializeOnSphere(N, r, u, m);
    for (int i = 0; i < N; i++)
    { 
        u[i][0] = 0;
        u[i][1] = 0;
        a[i][0] = 0;
        a[i][1] = 0;
        m[i] = 1;
        myfile << r[i][0] << "   "
               << r[i][1] << "   "
               << u[i][0] << "   "
               << u[i][1] << "\n";
    }

    const float T = 10;
    const float dt = 0.00001;
    const int Ntimesteps = T/dt + 1;

    for (int k = 0; k < N; k++)
    {
          ComputeF(k, r, m, a);
    }

    for (int t = 0; t < Ntimesteps; t++)
    {
        for (int j = 0; j < N; j++)
        {
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            r[j][0] += u[j][0] * dt;
            r[j][1] += u[j][1] * dt;
        }
        for (int j = 0; j < N; j++)
        {
            a[j][0] = 0;
            a[j][1] = 0;
        }
        for (int j = 0; j < N; j++)
        {
            ComputeF(j, r, m, a);
            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            if (t % 200 == 0)
            {
                myfile << r[j][0] << "   "
                       << r[j][1] << "   "
                       << u[j][0] << "   "
                       << u[j][1] << "\n";
            }
        }
    }
    return 0;
}
