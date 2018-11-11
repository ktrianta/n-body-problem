#include <cmath>
#include <ctime>
#include <cstdlib>
#include "octree.hpp"
#include <stdio.h> //for printf


bool Rectangle::containPoint(double point[])
{
    if (x - w > point[0]) return false;
    if (x + w < point[0]) return false;
    if (y - h > point[1]) return false;
    if (y + h < point[1]) return false;
    if (z - t > point[2]) return false;
    if (z + t < point[2]) return false;
    return true;
}



// Constructor


// Constructor for quadtree with particular size and parent -- build the tree, too!
Octree::Octree(double (*r)[3], double* m, int N, double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht)
{
    initialization(NULL, r, m, inp_x, inp_y, inp_z, inp_hw, inp_hh, inp_ht);
    fill(N);
}

// Constructor for quadtree with particular size and parent -- build the tree, too!
Octree::Octree(Octree* inp_parent, double (*r)[3], double* m, int N, double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht)
{
    initialization(inp_parent, r, m, inp_x, inp_y, inp_z, inp_hw, inp_hh, inp_ht);
    fill(N);
}

/* we dont need these two constructors
// Constructor for quadtree with particular size (do not fill the tree)
Octree::Octree(double (*r)[2], double inp_x, double inp_y, double inp_hw, double inp_hh)
{
    initialization(NULL, data, inp_x, inp_y, inp_hw, inp_hh);
}
*/
// Constructor for quadtree with particular size and parent (do not fill the tree)
Octree::Octree(Octree* inp_parent, double (*r)[3], double* m, double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht)
{
    initialization(inp_parent, r, m, inp_x, inp_y, inp_z, inp_hw, inp_hh, inp_ht);
}

Octree::~Octree()
{
    delete fnorthWest;
    delete fnorthEast;
    delete fsouthWest;
    delete fsouthEast;
    delete bnorthWest;
    delete bnorthEast;
    delete bsouthWest;
    delete bsouthEast;
}

void Octree::initialization(Octree* inp_parent, double (*r)[3], double* m, double inp_x, double inp_y, double inp_z, double inp_w, double inp_h, double inp_t){

    parent = inp_parent;
    data = r;
    mass = m;
    leaf = true;
    size = 0;
    cum_size = 0;
    boundary.x = inp_x;
    boundary.y = inp_y;
    boundary.z = inp_z;
    boundary.w = inp_w;
    boundary.h = inp_h;
    boundary.t = inp_t;
    fnorthWest = NULL;
    fnorthEast = NULL;
    fsouthWest = NULL;
    fsouthEast = NULL;
    bnorthWest = NULL;
    bnorthEast = NULL;
    bsouthWest = NULL;
    bsouthEast = NULL;
    centerMass[0] = 0.;
    centerMass[1] = 0.;
    centerMass[2] = 0.;
    cum_Mass = 0;
}

bool Octree::insert(int new_index){

    double point[3];
    point[0] = data[new_index][0];
    point[1] = data[new_index][1];
    point[2] = data[new_index][2];
    if (!boundary.containPoint(point)) return false;

    cum_size++;
    double mult1 = (double) (cum_size - 1) / (double) cum_size;
    double mult2 = 1.0 / (double) cum_size;
    for (int d = 0; d < dimension; d++) centerMass[d] *= mult1;
    for (int d = 0; d < dimension; d++) centerMass[d] += mult2 * point[d];
    cum_Mass += mass[new_index];
    
    if (leaf && size==0){
        index=new_index;
        size++;
        return true;
    }

    if (leaf) subdivide();

    if(fnorthWest->insert(new_index)) return true;
    if(fnorthEast->insert(new_index)) return true;
    if(fsouthWest->insert(new_index)) return true;
    if(fsouthEast->insert(new_index)) return true;
    if(bnorthWest->insert(new_index)) return true;
    if(bnorthEast->insert(new_index)) return true;
    if(bsouthWest->insert(new_index)) return true;
    if(bsouthEast->insert(new_index)) return true;

    return false;
}

void Octree::subdivide(){

    fnorthWest = new Octree(this, data, mass, boundary.x - .5 * boundary.w, boundary.y - .5 * boundary.h, boundary.z + .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    fnorthEast = new Octree(this, data, mass, boundary.x + .5 * boundary.w, boundary.y - .5 * boundary.h, boundary.z + .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    fsouthWest = new Octree(this, data, mass, boundary.x - .5 * boundary.w, boundary.y + .5 * boundary.h, boundary.z + .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    fsouthEast = new Octree(this, data, mass, boundary.x + .5 * boundary.w, boundary.y + .5 * boundary.h, boundary.z + .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    bnorthWest = new Octree(this, data, mass, boundary.x - .5 * boundary.w, boundary.y - .5 * boundary.h, boundary.z - .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    bnorthEast = new Octree(this, data, mass, boundary.x + .5 * boundary.w, boundary.y - .5 * boundary.h, boundary.z - .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    bsouthWest = new Octree(this, data, mass, boundary.x - .5 * boundary.w, boundary.y + .5 * boundary.h, boundary.z - .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    bsouthEast = new Octree(this, data, mass, boundary.x + .5 * boundary.w, boundary.y + .5 * boundary.h, boundary.z - .5 * boundary.t, .5 * boundary.w, .5 * boundary.h, .5 * boundary.t);
    
    // Move existing points to correct children
//    for(int i = 0; i < size; i++) {
        bool success = false;
        if (!success) success = fnorthWest->insert(index);
        if (!success) success = fnorthEast->insert(index);
        if (!success) success = fsouthWest->insert(index);
        if (!success) success = fsouthEast->insert(index);
        if (!success) success = bnorthWest->insert(index);
        if (!success) success = bnorthEast->insert(index);
        if (!success) success = bsouthWest->insert(index);
        if (!success) success = bsouthEast->insert(index);
        index = -1;
//    }
    
    // Empty parent node
    size = 0;
    leaf = false;

}

void Octree::fill(int N)
{
    for (int i = 0; i < N; i++) insert(i);
}

// Print out tree
void Octree::print() 
{
    if(cum_size == 0) {
        printf("Empty node\n");
        return;
    }

    if(leaf) {
        printf("Leaf node; data = [");
//        for(int i = 0; i < size; i++) {
            double point[3];
            point[0] = data[index][0];
            point[1] = data[index][1];
            point[2] = data[index][2];

            for(int d = 0; d < dimension; d++) printf("%f, ", point[d]);
//            if(i < size - 1) printf("\n");
            printf("]\n");
//        }        
    }
    else {
        printf("Intersection node with center-of-mass = [");
        for(int d = 0; d < dimension; d++) printf("%f, ", centerMass[d]);
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


void Octree::computeAcceleration(int idx, double (*r)[3], double (*a)[3], double g, double theta) {
    if (cum_size == 0 || (leaf == true && size == 1 && index == idx)) {
        return;
    }

    
    double d = sqrt( (centerMass[0] - r[idx][0]) * (centerMass[0] - r[idx][0])
                   + (centerMass[1] - r[idx][1]) * (centerMass[1] - r[idx][1])
                   + (centerMass[2] - r[idx][2]) * (centerMass[2] - r[idx][2]));

    if (2.*boundary.w/d <= theta) {
        
        double rji[3];
        rji[0] = centerMass[0] - r[idx][0];
        rji[1] = centerMass[1] - r[idx][1];
        rji[2] = centerMass[2] - r[idx][2];
        double r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];
        double denom = r2 * sqrt(r2);
        double a_i = -g * cum_Mass / denom;
        a[idx][0] -= a_i * rji[0];
        a[idx][1] -= a_i * rji[1];
        a[idx][2] -= a_i * rji[2];
        return;
    }
    if (leaf == true) {
//        for(int i = 0; i < size; i++) {
            int idx_i = index;
            if (idx == idx_i) return;

            double rji[3];
            rji[0] = r[idx_i][0] - r[idx][0];
            rji[1] = r[idx_i][1] - r[idx][1];
            rji[2] = r[idx_i][2] - r[idx][2];
            double r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];
            double denom = r2 * sqrt(r2);
            double a_i = -g * mass[idx_i] / denom;
            a[idx][0] -= a_i * rji[0];
            a[idx][1] -= a_i * rji[1];
            a[idx][2] -= a_i * rji[2];
//        }
    } else {
        fnorthWest->computeAcceleration(idx, r, a, g, theta);
        fnorthEast->computeAcceleration(idx, r, a, g, theta);
        fsouthWest->computeAcceleration(idx, r, a, g, theta);
        fsouthEast->computeAcceleration(idx, r, a, g, theta);
        bnorthWest->computeAcceleration(idx, r, a, g, theta);
        bnorthEast->computeAcceleration(idx, r, a, g, theta);
        bsouthWest->computeAcceleration(idx, r, a, g, theta);
        bsouthEast->computeAcceleration(idx, r, a, g, theta);
    }
}

