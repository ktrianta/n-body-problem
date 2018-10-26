#include <cmath>
#include <ctime>
#include <cstdlib>
#include "quadtree.hpp"
#include <stdio.h> //for printf


bool Rectangle::containPoint(double point[])
{
    if (x - w > point[0]) return false;
    if (x + w < point[0]) return false;
    if (y - h > point[1]) return false;
    if (y + h < point[1]) return false;
    return true;
}



// Constructor


// Constructor for quadtree with particular size and parent -- build the tree, too!
QuadTree::QuadTree(double (*r)[2], vector<double>& m, int N, double inp_x, double inp_y, double inp_hw, double inp_hh)
{
    initialization(NULL, r, m, inp_x, inp_y, inp_hw, inp_hh);
    fill(N);
}

// Constructor for quadtree with particular size and parent -- build the tree, too!
QuadTree::QuadTree(QuadTree* inp_parent, double (*r)[2], vector<double>& m, int N, double inp_x, double inp_y, double inp_hw, double inp_hh)
{
    initialization(inp_parent, r, m, inp_x, inp_y, inp_hw, inp_hh);
    fill(N);
}

/* we dont need these two constructors
// Constructor for quadtree with particular size (do not fill the tree)
QuadTree::QuadTree(double (*r)[2], double inp_x, double inp_y, double inp_hw, double inp_hh)
{
    initialization(NULL, data, inp_x, inp_y, inp_hw, inp_hh);
}
*/
// Constructor for quadtree with particular size and parent (do not fill the tree)
QuadTree::QuadTree(QuadTree* inp_parent, double (*r)[2], vector<double>& m, double inp_x, double inp_y, double inp_hw, double inp_hh)
{
    initialization(inp_parent, r, m, inp_x, inp_y, inp_hw, inp_hh);
}

QuadTree::~QuadTree()
{
    delete northWest;
    delete northEast;
    delete southWest;
    delete southEast;
}

void QuadTree::initialization(QuadTree* inp_parent, double (*r)[2], vector<double>& m, double inp_x, double inp_y, double inp_w, double inp_h){

    parent = inp_parent;
    data = r;
    mass = m;
    leaf = true;
    size = 0;
    cum_size = 0;
    boundary.x = inp_x;
    boundary.y = inp_y;
    boundary.w = inp_w;
    boundary.h = inp_h;
    northWest = NULL;
    northEast = NULL;
    southWest = NULL;
    southEast = NULL;
    centerMass[0] = 0.;
    centerMass[1] = 0.;
    cum_Mass = 0;
}

bool QuadTree::insert(int new_index){

    double point[2];
    point[0] = data[new_index][0];
    point[1] = data[new_index][1];
    if (!boundary.containPoint(point)) return false;

    cum_size++;
    double mult1 = (double) (cum_size - 1) / (double) cum_size;
    double mult2 = 1.0 / (double) cum_size;
    for (int d = 0; d < dimension; d++) centerMass[d] *= mult1;
    for (int d = 0; d < dimension; d++) centerMass[d] += mult2 * point[d];
    for (int k = 0; k < cum_size;  k++) cum_Mass += mass[k];
    
    if (leaf && size<capacity){
        index[size]=new_index;
        size++;
        return true;
    }

    if (leaf) subdivide();

    if(northWest->insert(new_index)) return true;
    if(northEast->insert(new_index)) return true;
    if(southWest->insert(new_index)) return true;
    if(southEast->insert(new_index)) return true;

    return false;
}

void QuadTree::subdivide(){

    northWest = new QuadTree(this, data, mass, boundary.x - .5 * boundary.w, boundary.y - .5 * boundary.h, .5 * boundary.w, .5 * boundary.h);
    northEast = new QuadTree(this, data, mass, boundary.x + .5 * boundary.w, boundary.y - .5 * boundary.h, .5 * boundary.w, .5 * boundary.h);
    southWest = new QuadTree(this, data, mass, boundary.x - .5 * boundary.w, boundary.y + .5 * boundary.h, .5 * boundary.w, .5 * boundary.h);
    southEast = new QuadTree(this, data, mass, boundary.x + .5 * boundary.w, boundary.y + .5 * boundary.h, .5 * boundary.w, .5 * boundary.h);
    
    // Move existing points to correct children
    for(int i = 0; i < size; i++) {
        bool success = false;
        if (!success) success = northWest->insert(index[i]);
        if (!success) success = northEast->insert(index[i]);
        if (!success) success = southWest->insert(index[i]);
        if (!success) success = southEast->insert(index[i]);
        index[i] = -1;
    }
    
    // Empty parent node
    size = 0;
    leaf = false;

}

void QuadTree::fill(int N)
{
    for (int i = 0; i < N; i++) insert(i);
}

// Print out tree
void QuadTree::print() 
{
    if(cum_size == 0) {
        printf("Empty node\n");
        return;
    }

    if(leaf) {
        printf("Leaf node; data = [");
        for(int i = 0; i < size; i++) {
            double point[2];
            point[0] = data[index[i]][0];
            point[1] = data[index[i]][1];

            for(int d = 0; d < dimension; d++) printf("%f, ", point[d]);
            if(i < size - 1) printf("\n");
            else printf("]\n");
        }        
    }
    else {
        printf("Intersection node with center-of-mass = [");
        for(int d = 0; d < dimension; d++) printf("%f, ", centerMass[d]);
        printf("]; children are:\n");
        northEast->print();
        northWest->print();
        southEast->print();
        southWest->print();
    }
}


void QuadTree::computeAcceleration(int idx, double (*r)[2], double (*a)[2], double g, double theta) {
    if (cum_size == 0 || (leaf == true && size == 1 && index[0] == idx)) {
        return;
    }

    
    double d = sqrt((centerMass[0]-r[idx][0])*(centerMass[0]-r[idx][0])+(centerMass[1]-r[idx][1])*(centerMass[1]-r[idx][1]));

    if (2.*boundary.w/d <= theta) {
        
        double rji[2];
        rji[0] = centerMass[0] - r[idx][0];
        rji[1] = centerMass[1] - r[idx][1];
        double r2 = rji[0] * rji[0] + rji[1] * rji[1];
        double denom = r2 * sqrt(r2);
        double a_i = -g * cum_Mass / denom;
        a[idx][0] -= a_i * rji[0];
        a[idx][1] -= a_i * rji[1];
        return;
    }
    if (leaf == true) {
        for(int i = 0; i < size; i++) {
            int idx_i = index[i];
            if (idx == idx_i) continue;

            double rji[2];
            rji[0] = r[idx_i][0] - r[idx][0];
            rji[1] = r[idx_i][1] - r[idx][1];
            double r2 = rji[0] * rji[0] + rji[1] * rji[1];
            double denom = r2 * sqrt(r2);
            double a_i = -g * mass[idx_i] / denom;
            a[idx][0] -= a_i * rji[0];
            a[idx][1] -= a_i * rji[1];
        }
    } else {
        northWest->computeAcceleration(idx, r, a, g, theta);
        northEast->computeAcceleration(idx, r, a, g, theta);
        southWest->computeAcceleration(idx, r, a, g, theta);
        southEast->computeAcceleration(idx, r, a, g, theta);
    }
}

