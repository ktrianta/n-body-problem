#include <cmath>
#include <ctime>
#include <cstdlib>
#include <stdio.h>
#include "octree.hpp"

Rectangle::Rectangle(sim::data_type px, sim::data_type py, sim::data_type pz,
                     sim::data_type pw, sim::data_type ph, sim::data_type pt)
    : x(px), y(py), z(pz), w(pw), h(ph), t(pt)
{
}

bool Rectangle::containPoint(const sim::data_type point[]) {
    if (x - w > point[0]) return false;
    if (x + w < point[0]) return false;
    if (y - h > point[1]) return false;
    if (y + h < point[1]) return false;
    if (z - t > point[2]) return false;
    if (z + t < point[2]) return false;
    return true;
}


// Initialize Octree static fields
sim::data_type (*Octree::r)[3] = NULL;
sim::data_type (*Octree::a)[3] = NULL;
sim::data_type *Octree::m = NULL;
sim::data_type Octree::g = 0;
sim::data_type Octree::theta = 0;


Octree::Octree(size_t N, sim::data_type px, sim::data_type py, sim::data_type pz,
               sim::data_type pw, sim::data_type ph, sim::data_type pt)
    : Octree(px, py, pz, pw, ph, pt)
{
    fill(N);
}

Octree::Octree(sim::data_type px, sim::data_type py, sim::data_type pz, sim::data_type pw,
               sim::data_type ph, sim::data_type pt)
    : boundary(px, py, pz, pw, ph, pt), massCenter { 0, 0, 0 }, mass(0), index(-1), leaf(true), cumSize(0),
      fnorthWest(NULL), fnorthEast(NULL), fsouthWest(NULL), fsouthEast(NULL),
      bnorthWest(NULL), bnorthEast(NULL), bsouthWest(NULL), bsouthEast(NULL)
{
    massCenter[0] = 0;
    massCenter[1] = 0;
    massCenter[2] = 0;
}

Octree::~Octree() {
    delete fnorthWest;
    delete fnorthEast;
    delete fsouthWest;
    delete fsouthEast;
    delete bnorthWest;
    delete bnorthEast;
    delete bsouthWest;
    delete bsouthEast;
}

bool Octree::insert(int new_index) {
    if (!boundary.containPoint(r[new_index])) return false;
    cumSize++;

    sim::data_type mult1 = (cumSize - 1) / (sim::data_type) cumSize;
    sim::data_type mult2 = 1.0 / cumSize;
    for (size_t d = 0; d < dimension; d++) massCenter[d] *= mult1;
    for (size_t d = 0; d < dimension; d++) massCenter[d] += mult2 * r[new_index][d];
    mass += m[new_index];

    if (leaf && cumSize == 1) {
        index = new_index;
        return true;
    }

    if (leaf) subdivide();

    if (fnorthWest->insert(new_index)) return true;
    if (fnorthEast->insert(new_index)) return true;
    if (fsouthWest->insert(new_index)) return true;
    if (fsouthEast->insert(new_index)) return true;
    if (bnorthWest->insert(new_index)) return true;
    if (bnorthEast->insert(new_index)) return true;
    if (bsouthWest->insert(new_index)) return true;
    if (bsouthEast->insert(new_index)) return true;

    return false;
}

void Octree::subdivide() {
    fnorthWest = new Octree(boundary.x - .5 * boundary.w, boundary.y - .5 * boundary.h,
        boundary.z + .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    fnorthEast = new Octree(boundary.x + .5 * boundary.w, boundary.y - .5 * boundary.h,
        boundary.z + .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    fsouthWest = new Octree(boundary.x - .5 * boundary.w, boundary.y + .5 * boundary.h,
        boundary.z + .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    fsouthEast = new Octree(boundary.x + .5 * boundary.w, boundary.y + .5 * boundary.h,
        boundary.z + .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    bnorthWest = new Octree(boundary.x - .5 * boundary.w, boundary.y - .5 * boundary.h,
        boundary.z - .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    bnorthEast = new Octree(boundary.x + .5 * boundary.w, boundary.y - .5 * boundary.h,
        boundary.z - .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    bsouthWest = new Octree(boundary.x - .5 * boundary.w, boundary.y + .5 * boundary.h,
        boundary.z - .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    bsouthEast = new Octree(boundary.x + .5 * boundary.w, boundary.y + .5 * boundary.h,
        boundary.z - .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    
    bool success = false;
    if (!success) success = fnorthWest->insert(index);
    if (!success) success = fnorthEast->insert(index);
    if (!success) success = fsouthWest->insert(index);
    if (!success) success = fsouthEast->insert(index);
    if (!success) success = bnorthWest->insert(index);
    if (!success) success = bnorthEast->insert(index);
    if (!success) success = bsouthWest->insert(index);
    if (!success) success = bsouthEast->insert(index);
    
    // Empty parent node
    index = -1;
    leaf = false;
}

void Octree::fill(size_t N) {
    for (size_t i = 0; i < N; i++) {
        insert(i);
    }
}

// Print out tree
void Octree::print()  {
    if(cumSize == 0) {
        printf("Empty node\n");
        return;
    }

    if (leaf) {
        printf("Leaf node; r = [");
        sim::data_type point[3];
        point[0] = r[index][0];
        point[1] = r[index][1];
        point[2] = r[index][2];

        for(int d = 0; d < dimension; d++) printf("%f, ", point[d]);
        printf("]\n");
    } else {
        printf("Intersection node with center-of-m = [");
        for(int d = 0; d < dimension; d++) printf("%f, ", massCenter[d]);
        printf("]; children are:\n");
        fnorthEast->print();
        fnorthWest->print();
        fsouthEast->print();
        fsouthWest->print();
        bnorthEast->print();
        bnorthWest->print();
        bsouthEast->print();
        bsouthWest->print();
    }
}


void Octree::computeAcceleration(int idx) {
    if (cumSize == 0 || index == idx) {
        return;
    }

    sim::data_type rji[3];
    rji[0] = massCenter[0] - r[idx][0];
    rji[1] = massCenter[1] - r[idx][1];
    rji[2] = massCenter[2] - r[idx][2];
    sim::data_type r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];
    sim::data_type d = sqrt(r2+sim::e2);

    if (leaf == true) {
        sim::data_type denom = (r2+sim::e2) * d;
        sim::data_type a_i = -g * mass / denom;
        a[idx][0] -= a_i * rji[0];
        a[idx][1] -= a_i * rji[1];
        a[idx][2] -= a_i * rji[2];
        return;
    }

    if (2.*boundary.w/d <= theta) {
        sim::data_type denom = (r2+sim::e2) * d;
        sim::data_type a_i = -g * mass / denom;
        a[idx][0] -= a_i * rji[0];
        a[idx][1] -= a_i * rji[1];
        a[idx][2] -= a_i * rji[2];
        return;
    } else {
        fnorthWest->computeAcceleration(idx);
        fnorthEast->computeAcceleration(idx);
        fsouthWest->computeAcceleration(idx);
        fsouthEast->computeAcceleration(idx);
        bnorthWest->computeAcceleration(idx);
        bnorthEast->computeAcceleration(idx);
        bsouthWest->computeAcceleration(idx);
        bsouthEast->computeAcceleration(idx);
    }
}
