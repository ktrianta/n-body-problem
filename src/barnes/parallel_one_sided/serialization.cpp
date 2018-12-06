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

    
void Treenode::set(sim::data_type px, sim::data_type py, sim::data_type pz, sim::data_type pw, sim::data_type ph, sim::data_type pt) {
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

Serialization::Serialization()
    : size(0), treeArray(NULL)
    {}

Serialization::Serialization(sim::data_type px, sim::data_type py, sim::data_type pz, sim::data_type pw, sim::data_type ph, sim::data_type pt)
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

void Serialization::insert(int j, sim::data_type rx, sim::data_type ry, sim::data_type rz, sim::data_type m) {
    queue<size_t> list;
    list.push(0);

    //std::cout << " ---------------- insert " << j << std::endl;
    while (!list.empty()) {
        size_t f = list.front();
        //std::cout << "popping " << f << std::endl;
        list.pop();

        //std::cout << list.size() << std::endl;

        const sim::data_type x = treeArray[f].x;
        const sim::data_type y = treeArray[f].y;
        const sim::data_type z = treeArray[f].z;
        const sim::data_type w = treeArray[f].w;
        const sim::data_type h = treeArray[f].h;
        const sim::data_type t = treeArray[f].t;

        if (x - w > rx) continue;
        if (x + w < rx) continue;
        if (y - h > ry) continue;
        if (y + h < ry) continue;
        if (z - t > rz) continue;
        if (z + t < rz) continue;

        size_t cum_size = ++treeArray[f].cum_size;

        if (treeArray[f].leaf && cum_size == 1) {
            //std::cout << "leaf " << f << " index " << j << std::endl;
            treeArray[f].index = j;
            treeArray[f].massCenter[0] = rx;
            treeArray[f].massCenter[1] = ry;
            treeArray[f].massCenter[2] = rz;
            treeArray[f].mass = m;
            return;
        }

        if (treeArray[f].leaf) {
            //std::cout << f << " sub " << rx << " " << ry << std::endl;
            //std::cout << "cen " << treeArray[f].massCenter[0] << " " << treeArray[f].massCenter[1] << std::endl;
            subdivide(f);
        }

        sim::data_type mult1 = (sim::data_type) (cum_size - 1) / cum_size;
        sim::data_type mult2 = 1.0 / (sim::data_type) cum_size;
        //std::cout << "mult1 " << mult1 << std::endl;
        //std::cout << "mult2 " << mult2 << std::endl;
        for (int d = 0; d < 3; d++) treeArray[f].massCenter[d] *= mult1;
        treeArray[f].massCenter[0] += mult2 * rx;
        treeArray[f].massCenter[1] += mult2 * ry;
        treeArray[f].massCenter[2] += mult2 * rz;
        //std::cout << f << " cen " << treeArray[f].massCenter[0] << " " << treeArray[f].massCenter[1] << std::endl;
        treeArray[f].mass += m;

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

bool Serialization::insertInLeaf(int current, int index, sim::data_type rx, sim::data_type ry, sim::data_type rz, sim::data_type m) {
    //std::cout << "!!!" << std::endl;
    const sim::data_type x = treeArray[current].x;
    const sim::data_type y = treeArray[current].y;
    const sim::data_type z = treeArray[current].z;
    const sim::data_type w = treeArray[current].w;
    const sim::data_type h = treeArray[current].h;
    const sim::data_type t = treeArray[current].t;

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
    //std::cout << "index " << treeArray[current].index << " position " << rx << " " << ry << " " << rz << std::endl;
    //std::cout << x << " " << y << " " << z << " " << w << " " << h << " " << t << std::endl;
    //std::cout << index << " in leaf " << current << std::endl;
    return true;
}

void Serialization::subdivide(int current) {
    if (position + 8 >= size) {
        //std::cout << "size " << size*2 << std::endl;
        size_t old_size = size;
        size *= 2;
        Treenode *temp = new Treenode[size];

        std::copy(treeArray, treeArray + old_size, temp);
        delete[] treeArray;
        treeArray = temp;
    }

    const sim::data_type x = treeArray[current].x;
    const sim::data_type y = treeArray[current].y;
    const sim::data_type z = treeArray[current].z;
    const sim::data_type w = treeArray[current].w;
    const sim::data_type h = treeArray[current].h;
    const sim::data_type t = treeArray[current].t;
    treeArray[position+1].set(x - 0.5*w, y - 0.5*h, z + 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+2].set(x + 0.5*w, y - 0.5*h, z + 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+3].set(x - 0.5*w, y + 0.5*h, z + 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+4].set(x + 0.5*w, y + 0.5*h, z + 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+5].set(x - 0.5*w, y - 0.5*h, z - 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+6].set(x + 0.5*w, y - 0.5*h, z - 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+7].set(x - 0.5*w, y + 0.5*h, z - 0.5*t, 0.5*w, 0.5*h, 0.5*t);  
    treeArray[position+8].set(x + 0.5*w, y + 0.5*h, z - 0.5*t, 0.5*w, 0.5*h, 0.5*t);  


    int index = treeArray[current].index;
    sim::data_type rx = treeArray[current].massCenter[0];
    sim::data_type ry = treeArray[current].massCenter[1];
    sim::data_type rz = treeArray[current].massCenter[2];

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

void Serialization::computeAcceleration(int idx, sim::data_type r[3], sim::data_type a[3], sim::data_type g, sim::data_type theta, size_t& calcs) {
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
            calcs += 1;
            continue;
        }

        if (2 * current.w / d <= theta) {
            sim::data_type denom = (r2+sim::e2) * d;
            sim::data_type a_i = -g * current.mass / denom;
            a[0] -= a_i * rji[0];
            a[1] -= a_i * rji[1];
            a[2] -= a_i * rji[2];
            calcs += 1;
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
