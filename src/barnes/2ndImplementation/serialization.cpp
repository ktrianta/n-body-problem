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

Serialization::~Serialization() {
    delete[] treeArray;
}

void Serialization::insert(int j, double rx, double ry, double rz, double m) {
    queue<size_t> list;
    list.push(0);

    while (!list.empty()) {
        size_t f = list.front();
        std::cout << f << std::endl;
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
        double mult1 = (double) (cum_size - 1) / (double) cum_size;
        double mult2 = 1.0 / (double) cum_size;
        for (int d = 0; d < 3; d++) treeArray[f].massCenter[d] *= mult1;
        treeArray[f].massCenter[0] += mult2 * rx;
        treeArray[f].massCenter[1] += mult2 * ry;
        treeArray[f].massCenter[2] += mult2 * rz;
        treeArray[f].mass += m;

        if (treeArray[f].leaf && treeArray[f].cum_size == 1) {
            treeArray[f].index = j;
            std::cout << "leaf " << j << std::endl;
            return;
        }

        if (treeArray[f].leaf) {
            std::cout << "sub " << f << std::endl;
           subdivide(f);
           treeArray[f].child = position+1;
        }

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

bool Serialization::insertInLeaf(int current, int index, double rx, double ry, double rz, double m) {
    std::cout << "!!!" << std::endl;
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

    size_t cum_size = 1;
    treeArray[current].massCenter[0] = rx;
    treeArray[current].massCenter[1] = ry;
    treeArray[current].massCenter[2] = rz;
    treeArray[current].mass = m;
    treeArray[current].index = index;
    return true;
}

void Serialization::subdivide(int current) {
    if (position + 8 > size) {
        size *= 2;
        Treenode *temp = new Treenode[size];
        std::copy(treeArray, treeArray + size/2, temp);
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

    insertInLeaf(position+1, treeArray[current].index, x, y, z, treeArray[current].mass);
    insertInLeaf(position+2, treeArray[current].index, x, y, z, treeArray[current].mass);
    insertInLeaf(position+3, treeArray[current].index, x, y, z, treeArray[current].mass);
    insertInLeaf(position+4, treeArray[current].index, x, y, z, treeArray[current].mass);
    insertInLeaf(position+5, treeArray[current].index, x, y, z, treeArray[current].mass);
    insertInLeaf(position+6, treeArray[current].index, x, y, z, treeArray[current].mass);
    insertInLeaf(position+7, treeArray[current].index, x, y, z, treeArray[current].mass);
    insertInLeaf(position+8, treeArray[current].index, x, y, z, treeArray[current].mass);

    treeArray[current].leaf = false;
    treeArray[current].index = -1;
    position += 8;
}

void Serialization::computeAcceleration(int current, int idx, double (*r)[3], double (*a)[3], double g, double theta) {
    if (treeArray[current].cum_size == 0 || (treeArray[current].leaf == true && treeArray[current].cum_size == 1 && treeArray[current].index == idx)) {                           
        return;                                                                                   
    }                                                                                             

    //std::cout << current << std::endl;                                                                                      
    double centerMass[3];
    centerMass[0] = treeArray[current].massCenter[0];
    centerMass[1] = treeArray[current].massCenter[1];
    centerMass[2] = treeArray[current].massCenter[2];

    double d = sqrt( (centerMass[0] - r[idx][0]) * (centerMass[0] - r[idx][0])                    
                   + (centerMass[1] - r[idx][1]) * (centerMass[1] - r[idx][1])
                   + (centerMass[2] - r[idx][2]) * (centerMass[2] - r[idx][2]));                  
    
    if (2.*treeArray[current].w/d <= theta) {                                                               
        std::cout << "NO" << std::endl;
        double rji[3];
        rji[0] = centerMass[0] - r[idx][0];                                                       
        rji[1] = centerMass[1] - r[idx][1];                                                       
        rji[2] = centerMass[2] - r[idx][2];
        double r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];                          
        double denom = r2 * sqrt(r2);
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
            rji[0] = r[idx_i][0] - r[idx][0];                                                     
            rji[1] = r[idx_i][1] - r[idx][1];                                                     
            rji[2] = r[idx_i][2] - r[idx][2];
            double r2 = rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2];                      
            double denom = r2 * sqrt(r2); 
            double a_i = -g * treeArray[current].mass / denom;
            a[idx][0] -= a_i * rji[0];
            printf("%f \n" , a[idx][0]);
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
