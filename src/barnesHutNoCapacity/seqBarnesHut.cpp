#include <math.h> //for sqrt function
#include <stdio.h> //for printf
#include <stdlib.h> //for srand/rand prng
#include <time.h> //for seeding prng
#include <iostream>
#include <fstream>
#include <vector>
#include "../utils/types.hpp"
#include "../utils/initialization.hpp"
#include "quadtree.hpp"
#include <unistd.h>

using namespace std;

int N = 5;      // the number of particles
double theta = 0.0;      // the number of particles
const sim_data_type g = 1;     // gravitational constant
const sim_data_type epsilon = 0.001;
const sim_data_type epsilon2 = epsilon * epsilon;


void writeDataToFile(sim_data_type (*r)[3], sim_data_type (*u)[3], ofstream& file)
{
    for (int i = 0; i < N; i++)
    {
        file << r[i][0] << "   "
             << r[i][1] << "   "
             << r[i][2] << "   "
             << u[i][0] << "   "
             << u[i][1] << "   "
             << u[i][2] << "\n";
    }
}
int main(int argc, char** argv)
{


//  the center of the parent node and the half width and height
    double xc, yc, zc, h2, w2, t2;
    xc=0.;
    yc=0.;
    zc=0.;
    w2=4.;
    h2=4.;
    t2=4.;

    int c;
    sim_data_type T = 10;
    sim_data_type dt = 0.00001;
    string filename;

    while ((c = getopt (argc, argv, "n:t:s:i:h:")) != -1)
    {
        switch (c)
        {
            case 'n':
                N = atoi(optarg);
                break;
            case 't':
                T = atof(optarg);
                break;
            case 's':
                dt = atof(optarg);
                break;
            case 'i':
                filename = optarg;
                break;
            case 'h':
                theta = atof(optarg);
                break;
        }
    }

    sim_data_type (*r)[3] = new sim_data_type[N][3];
    sim_data_type (*u)[3] = new sim_data_type[N][3];
    sim_data_type (*a)[3] = new sim_data_type[N][3];
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);

    vector<sim_data_type> m(N, 1.0/N);

    if (!filename.empty()) {
        ifstream ifile;
        ifile.open(filename);

        for (int i = 0; i < N; i++) {
            ifile >> m[i] >> r[i][0] >> r[i][1] >> r[i][2] >> u[i][0] >> u[i][1] >> u[i][2];
        }
    } else {
          initializePositionOnSphere(N, r);
    }

    ofstream file;
    file.open("output.dat");

    writeDataToFile(r, u, file);
    
    QuadTree tree = QuadTree(r, m, N, xc, yc, zc, w2, h2, t2);
    for (int j = 0; j < N; j++)
    {
        tree.computeAcceleration(j, r, a, g, theta);
    }
    const int Ntimesteps = T/dt + 1;

    for (int t = 0; t < Ntimesteps; t++)                                                         
    {                         
        for (int j = 0; j < N; j++)                                                              
        {                                                                                    
            u[j][0] += 0.5 * a[j][0] * dt;                                                  
            u[j][1] += 0.5 * a[j][1] * dt;                                                        
            u[j][2] += 0.5 * a[j][2] * dt;                                                        
            r[j][0] += u[j][0] * dt;                                                              
            r[j][1] += u[j][1] * dt;                                                      
            r[j][2] += u[j][2] * dt;                                                      
                                                                                                 
            tree.computeAcceleration(j, r, a, g, theta);                                                             

                                                                                         
            u[j][0] += 0.5 * a[j][0] * dt;                           
            u[j][1] += 0.5 * a[j][1] * dt;                       
            u[j][2] += 0.5 * a[j][2] * dt;                       
            a[j][0] = 0;
            a[j][1] = 0;
            a[j][2] = 0;
        }                                                                                         
        QuadTree tree = QuadTree(r, m, N, xc, yc, zc, w2, h2, t2);
                                                                                                  
        if (t % 600 == 0)                                                                         
        {                                                       
            writeDataToFile(r, u, file);                       
        }                                                                                         
    }       


//  tree.print();

    return 0;
}
