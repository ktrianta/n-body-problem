#include <cmath>
#include <ctime>
#include <cstdlib>
#include <boxComputation.hpp>




void boxComputation(int N, sim::data_type (*r)[3], sim::data_type &x, sim::data_type &y, sim::data_type &z, sim::data_type &w, sim::data_type &h, sim::data_type &t){


    sim::data_type xmax, xmin, ymax, ymin, zmax, zmin;
    xmax = r[0][0];
    xmin = r[0][0];
    ymax = r[0][1];
    ymin = r[0][1];
    zmax = r[0][2];
    zmin = r[0][2];
    for (int i = 1; i < N; i++) 
    {
        if (r[i][0] > xmax)
            xmax = r[i][0];
        if (r[i][0] < xmin)
            xmin = r[i][0];
        if (r[i][1] > ymax)
            ymax = r[i][1];
        if (r[i][1] < ymin)
            ymin = r[i][1];
        if (r[i][2] > zmax)
            zmax = r[i][2];
        if (r[i][2] < zmin)
            zmin = r[i][2];
    }
    w = (xmax-xmin+0.05)/2.;
    x = (xmax+xmin)/2.;
    h = (ymax-ymin+0.05)/2.;
    y = (ymax+ymin)/2.;
    t = (zmax-zmin+0.05)/2.;
    z = (zmax+zmin)/2.;
}
