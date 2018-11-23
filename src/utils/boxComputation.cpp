#include <boxComputation.hpp>

void boxComputation(const size_t N, sim::data_type (*r)[7], sim::data_type &x, sim::data_type &y, sim::data_type &z, sim::data_type &w, sim::data_type &h, sim::data_type &t) {
    sim::data_type xmax = r[0][1];
    sim::data_type xmin = r[0][1];
    sim::data_type ymax = r[0][2];
    sim::data_type ymin = r[0][2];
    sim::data_type zmax = r[0][3];
    sim::data_type zmin = r[0][3];

    for (int i = 1; i < N; i++)  {
        xmax = std::max(xmax, r[i][1]);
        xmin = std::min(xmin, r[i][1]);
        ymax = std::max(ymax, r[i][2]);
        ymin = std::min(ymin, r[i][2]);
        zmax = std::max(zmax, r[i][3]);
        zmin = std::min(zmin, r[i][3]);
    }

    w = (xmax-xmin+0.05) / 2;
    x = (xmax+xmin) / 2; 
    h = (ymax-ymin+0.05) / 2;
    y = (ymax+ymin) / 2;
    t = (zmax-zmin+0.05) / 2;
    z = (zmax+zmin) / 2;
}


void boxComputation(const size_t N, sim::data_type (*r)[3], sim::data_type &x, sim::data_type &y, sim::data_type &z, sim::data_type &w, sim::data_type &h, sim::data_type &t) {
    sim::data_type xmax = r[0][0];
    sim::data_type xmin = r[0][0];
    sim::data_type ymax = r[0][1];
    sim::data_type ymin = r[0][1];
    sim::data_type zmax = r[0][2];
    sim::data_type zmin = r[0][2];

    for (int i = 1; i < N; i++)  {
        xmax = std::max(xmax, r[i][0]);
        xmin = std::min(xmin, r[i][0]);
        ymax = std::max(ymax, r[i][1]);
        ymin = std::min(ymin, r[i][1]);
        zmax = std::max(zmax, r[i][2]);
        zmin = std::min(zmin, r[i][2]);
    }

    w = (xmax-xmin+0.05) / 2;
    x = (xmax+xmin) / 2; 
    h = (ymax-ymin+0.05) / 2;
    y = (ymax+ymin) / 2;
    t = (zmax-zmin+0.05) / 2;
    z = (zmax+zmin) / 2;
}
