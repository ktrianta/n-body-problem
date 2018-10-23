#include <math.h> //for sqrt function
#include <stdio.h> //for printf
#include <stdlib.h> //for srand/rand prng
#include <time.h> //for seeding prng
#include <iostream>
#include <fstream>
#include <vector>
#include "quadtree.hpp"
#include <unistd.h>

using namespace std;

int N =5;      // the number of particles


int main(int argc, char** argv)
{
    double xc,yc,h2,w2;
    xc=0.;
    yc=0.;
    w2=1.;
    h2=1.;
    double (*r)[2] = new double[N][2];
    r[0][0]=0.8;
    r[0][1]=0.8;
    r[1][0]=0.8;
    r[1][1]=-0.8;
    r[2][0]=-0.8;
    r[2][1]=-0.8;
    r[3][0]=-0.8;
    r[3][1]=0.8;
    r[4][0]=-0.4;
    r[4][1]=0.5;
    r[5][0]=-0.6;
    r[5][1]=0.1;
    QuadTree tree = QuadTree(r,N,xc,yc,w2,h2);
    printf("tasos");

    return 0;
}
