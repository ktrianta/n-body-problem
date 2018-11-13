#include <iostream>
#include <fstream>
#include <vector>
#include "io.hpp"
#include "types.hpp"
#include "initialization.hpp"
#include "serialization.hpp"
#include <unistd.h>

using namespace std;

int main(int argc, char** argv)
{

    int c;
    int N = 5;
    sim::data_type theta = 0.0;
    sim::data_type T = 10;
    sim::data_type dt = 0.00001;
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

    sim::data_type *m = new sim::data_type[N];
    sim::data_type (*r)[3] = new sim::data_type[N][3];
    sim::data_type (*u)[3] = new sim::data_type[N][3];
    sim::data_type (*a)[3] = new sim::data_type[N][3];
    std::fill(m, m+N, 1.0/N);
    std::fill(&u[0][0], &u[0][0] + N*3, 0);
    std::fill(&a[0][0], &a[0][0] + N*3, 0);
 
    if (!filename.empty()) {
        if (readDataFromFile(filename, N, m, r, u) == -1) {
            std::cerr << "File " << filename << " not found!" << std::endl;
        }
    } else {
          initializePositionOnSphere(N, r);
    }

//  the center of the parent node and the half width and height
    double xc, yc, zc, h2, w2, t2;
// coordinates for tab8096
    double maxX, minX, maxY, minY, maxZ, minZ;
    
    maxX = r[0][0];
    minX = r[0][0];
    maxY = r[0][1];
    minY = r[0][1];
    maxZ = r[0][2];
    minZ = r[0][2];
    for (int i = 1; i < N; i++) 
    {
        if (r[i][0] > maxX)
            maxX = r[i][0];
        if (r[i][0] < minX)
            minX = r[i][0];
        if (r[i][1] > maxY)
            maxY = r[i][1];
        if (r[i][1] < minY)
            minY = r[i][1];
        if (r[i][2] > maxZ)
            maxZ = r[i][2];
        if (r[i][2] < minZ)
            minZ = r[i][2];
    }
    w2 = (maxX-minX+0.05)/2.;
    xc = (maxX+minX)/2.;
    h2 = (maxY-minY+0.05)/2.;
    yc = (maxY+minY)/2.;
    t2 = (maxZ-minZ+0.05)/2.;
    zc = (maxZ+minZ)/2.;

    ofstream file;
    file.open("output.dat");

    writeDataToFile(N, r, u, file);


    Serialization* tree = new Serialization(0, 0, 0, 50, 50, 50);

    for (int i = 0; i < N; i++)
    {
        tree->insert(i, r[i][0], r[i][1], r[i][2], m[i]);
    }

    for (int j = 0; j < N; j++)
    {
        tree->computeAcceleration(0, j, r, a, sim::g, theta);
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

            a[j][0] = 0;
            a[j][1] = 0;
            a[j][2] = 0;
         
            tree->computeAcceleration(0, j, r, a, sim::g, theta);

            u[j][0] += 0.5 * a[j][0] * dt;
            u[j][1] += 0.5 * a[j][1] * dt;
            u[j][2] += 0.5 * a[j][2] * dt;
        }

        delete tree;
        tree = new Serialization(0, 0, 0, 50, 50, 50);

        //std::cout << "HEY" << std::endl;
        for (int i = 0; i < N; i++) {
            tree->insert(i, r[i][0], r[i][1], r[i][2], m[i]);
        }

        if (t % 200 == 0)
        {
            writeDataToFile(N, r, u, file);
        }
    }

    return 0;
}
