#include <cmath>
#include <ctime>
#include <cstdlib>
#include "serialization.hpp"
#include <stdio.h> //for printf
#include <queue>
#include <iostream>

Treenode::Treenode()
    : x(0), y(0), z(0), w(0), h(0), t(0), mass(0), index(-1), leaf(true), cum_size(0), child(0)
    {}

    
void Treenode::set(double px, double py, double pz, double pw, double ph, double pt) {
    x = px;
    y = py;
    z = pz;
    w = pw;
    h = ph;
    t = pt;
    mass = 0;
    index = -1;
    leaf = true;
    cum_size = 0;
    child = 0;
    massCenter[0] = -1;
    massCenter[0] = -1;
    massCenter[0] = -1;
}


Serialization::Serialization(double px, double py, double pz, double pw, double ph, double pt)
    : size(64), treeArray(new Treenode[size])
{
    treeArray[0].x = px;
    treeArray[0].y = py;
    treeArray[0].z = pz;
    treeArray[0].w = pw;
    treeArray[0].h = ph;
    treeArray[0].t = pt;
    position = 0;
}

Serialization::Serialization()
    : size(0), treeArray(NULL)
    {}


Serialization::~Serialization() {
    delete[] treeArray;
}

void Serialization::insert(int j, double rx, double ry, double rz, double m) {
    queue<size_t> list;
    list.push(0);

    while (!list.empty()) {
        size_t f = list.front();
        list.pop();


        const double x = treeArray[f].x;
        const double y = treeArray[f].y;
        const double z = treeArray[f].z;
        const double w = treeArray[f].w;
        const double h = treeArray[f].h;
        const double t = treeArray[f].t;

        if (x - w > rx) continue;
        if (x + w < rx) continue;
        if (y - h > ry) continue;
        if (y + h < ry) continue;
        if (z - t > rz) continue;
        if (z + t < rz) continue;

        size_t cum_size = ++treeArray[f].cum_size;

        if (treeArray[f].leaf && cum_size == 1) {
            treeArray[f].index = j;
            treeArray[f].massCenter[0] = rx;
            treeArray[f].massCenter[1] = ry;
            treeArray[f].massCenter[2] = rz;
            treeArray[f].mass = m;
            return;
        }

        if (treeArray[f].leaf) {
            subdivide(f);
            double mult1 = (double) (cum_size - 1) / cum_size;
            double mult2 = 1.0 / (double) cum_size;
            for (int d = 0; d < 3; d++) treeArray[f].massCenter[d] *= mult1;
            treeArray[f].massCenter[0] += mult2 * rx;
            treeArray[f].massCenter[1] += mult2 * ry;
            treeArray[f].massCenter[2] += mult2 * rz;
            treeArray[f].mass += m;
        }

        if (treeArray[f].child != 0) {
            list.push(treeArray[f].child);
            list.push(treeArray[f].child+1);
            list.push(treeArray[f].child+2);
            list.push(treeArray[f].child+3);
            list.push(treeArray[f].child+4);
            list.push(treeArray[f].child+5);
            list.push(treeArray[f].child+6);
            list.push(treeArray[f].child+7);
        }
    }
}

bool Serialization::insertInLeaf(int current, int index, double rx, double ry, double rz, double m) {
    const double x = treeArray[current].x;
    const double y = treeArray[current].y;
    const double z = treeArray[current].z;
    const double w = treeArray[current].w;
    const double h = treeArray[current].h;
    const double t = treeArray[current].t;

    if (x - w > rx) return false;
    if (x + w < rx) return false;
    if (y - h > ry) return false;
    if (y + h < ry) return false;
    if (z - t > rz) return false;
    if (z + t < rz) return false;

    treeArray[current].cum_size = 1;
    treeArray[current].massCenter[0] = rx;
    treeArray[current].massCenter[1] = ry;
    treeArray[current].massCenter[2] = rz;
    treeArray[current].mass = m;
    treeArray[current].index = index;
    return true;
}

void Serialization::subdivide(int current) {
    if (position + 8 >= size) {
        size_t old_size = size;
        size *= 2;
        Treenode *temp = new Treenode[size];
        std::copy(treeArray, treeArray + old_size, temp);
        delete[] treeArray;
        treeArray = temp;
    }

    const double x = treeArray[current].x;
    const double y = treeArray[current].y;
    const double z = treeArray[current].z;
    const double w = treeArray[current].w;
    const double h = treeArray[current].h;
    const double t = treeArray[current].t;
    treeArray[position+1].set(x - 0.5*w, y - 0.5*h, z + 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+2].set(x + 0.5*w, y - 0.5*h, z + 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+3].set(x - 0.5*w, y + 0.5*h, z + 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+4].set(x + 0.5*w, y + 0.5*h, z + 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+5].set(x - 0.5*w, y - 0.5*h, z - 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+6].set(x + 0.5*w, y - 0.5*h, z - 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+7].set(x - 0.5*w, y + 0.5*h, z - 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+8].set(x + 0.5*w, y + 0.5*h, z - 0.5*t, 0.5*w, 0.5*h, 0.5*t);  


    int index = treeArray[current].index;
    double rx = treeArray[current].massCenter[0];
    double ry = treeArray[current].massCenter[1];
    double rz = treeArray[current].massCenter[2];

    treeArray[current].leaf = false;
    treeArray[current].index = -1;
    treeArray[current].child = position + 1;
    position += 8;

    if (insertInLeaf(position-7, index, rx, ry, rz, treeArray[current].mass)) return;
    if (insertInLeaf(position-6, index, rx, ry, rz, treeArray[current].mass)) return;
    if (insertInLeaf(position-5, index, rx, ry, rz, treeArray[current].mass)) return;
    if (insertInLeaf(position-4, index, rx, ry, rz, treeArray[current].mass)) return;
    if (insertInLeaf(position-3, index, rx, ry, rz, treeArray[current].mass)) return;
    if (insertInLeaf(position-2, index, rx, ry, rz, treeArray[current].mass)) return;
    if (insertInLeaf(position-1, index, rx, ry, rz, treeArray[current].mass)) return;
    if (insertInLeaf(position, index, rx, ry, rz, treeArray[current].mass)) return;
}

/*
void Serialization::computeAcceleration(int current, int idx, double (*r)[7], double (*a)[3], double g, double theta) {
    if (treeArray[current].cum_size == 0 || (treeArray[current].leaf == true && treeArray[current].cum_size == 1 && treeArray[current].index == idx)) {                           
        return;                                                                                   
    }                                                                                             
    
    double centerMass[3];
    centerMass[0] = treeArray[current].massCenter[0];
    centerMass[1] = treeArray[current].massCenter[1];
    centerMass[2] = treeArray[current].massCenter[2];

    double d = sqrt( (centerMass[0] - r[idx][1]) * (centerMass[0] - r[idx][1])                    
                   + (centerMass[1] - r[idx][2]) * (centerMass[1] - r[idx][2])
                   + (centerMass[2] - r[idx][3]) * (centerMass[2] - r[idx][3]));                  
    
    if (2.*treeArray[current].w/d <= theta) {                                                               
        double rji[3];
        rji[0] = centerMass[0] - r[idx][1];                                                       
        rji[1] = centerMass[1] - r[idx][2];                                                       
        rji[2] = centerMass[2] - r[idx][3];
        double r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];                          
        double denom = r2+sim::e2 * sqrt(r2+sim::e2);
        double a_i = -g * treeArray[current].mass / denom;                                                       
        a[idx][0] -= a_i * rji[0];                                                                
        a[idx][1] -= a_i * rji[1];                                                                
        a[idx][2] -= a_i * rji[2];                                                                
        return;                                                                                   
    }
    if (treeArray[current].leaf == true) {
            int idx_i = treeArray[current].index;
            if (idx == idx_i) return;                                                             
            
            double rji[3];
            rji[0] = r[idx_i][1] - r[idx][1];                                                     
            rji[1] = r[idx_i][2] - r[idx][2];                                                     
            rji[2] = r[idx_i][3] - r[idx][3];
            double r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];                      
            double denom = r2+sim::e2 * sqrt(r2+sim::e2); 
            double a_i = -g * treeArray[current].mass / denom;
            a[idx][0] -= a_i * rji[0];
            a[idx][1] -= a_i * rji[1];
            a[idx][2] -= a_i * rji[2];
    } else {
        current = treeArray[current].child;
        computeAcceleration(current, idx, r, a, g, theta);
        computeAcceleration(current+1, idx, r, a, g, theta);
        computeAcceleration(current+2, idx, r, a, g, theta);
        computeAcceleration(current+3, idx, r, a, g, theta);
        computeAcceleration(current+4, idx, r, a, g, theta);
        computeAcceleration(current+5, idx, r, a, g, theta);
        computeAcceleration(current+6, idx, r, a, g, theta);
        computeAcceleration(current+7, idx, r, a, g, theta);
    }
}
*/

void Serialization::computeAcceleration(int idx, sim::data_type r[3], sim::data_type a[3], sim::data_type g, sim::data_type theta) {
    queue<size_t> list; list.push(0);
    sim::data_type rji[3];
    while (!list.empty()) {
        size_t c = list.front(); list.pop();
        const Treenode& current = treeArray[c];

        if (current.cum_size == 0 || current.index == idx) {
            continue;
        }

        rji[0] = current.massCenter[0] - r[0];
        rji[1] = current.massCenter[1] - r[1];
        rji[2] = current.massCenter[2] - r[2];
        sim::data_type r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];
        sim::data_type d = sqrt(r2+sim::e2);

        if (current.leaf) {
            sim::data_type denom = (r2+sim::e2) * d;
            sim::data_type a_i = -g * current.mass / denom;
            a[0] -= a_i * rji[0];
            a[1] -= a_i * rji[1];
            a[2] -= a_i * rji[2];
            continue;
        }

        if (2 * current.w / d <= theta) {
            sim::data_type denom = (r2+sim::e2) * d;
            sim::data_type a_i = -g * current.mass / denom;
            a[0] -= a_i * rji[0];
            a[1] -= a_i * rji[1];
            a[2] -= a_i * rji[2];
        } else {
            list.push(current.child);
            list.push(current.child+1);
            list.push(current.child+2);
            list.push(current.child+3);
            list.push(current.child+4);
            list.push(current.child+5);
            list.push(current.child+6);
            list.push(current.child+7);
        }
    }
}
