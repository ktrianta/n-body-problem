#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include "types.hpp"

using namespace std;


struct Treenode {
sim::data_type x;
sim::data_type y;
sim::data_type z;
sim::data_type w;
sim::data_type h;
sim::data_type t;
sim::data_type root;
sim::data_type massCenter[3];
sim::data_type mass;
int index;
bool leaf;
size_t cum_size;
size_t child;

Treenode();
void set(sim::data_type px, sim::data_type py, sim::data_type pz, sim::data_type pw, sim::data_type ph, sim::data_type pt);
};



class Serialization{
public: 
    size_t size;
    size_t position;
    Treenode *treeArray;

    Serialization();
    Serialization(sim::data_type px, sim::data_type py, sim::data_type pz, sim::data_type pw, sim::data_type ph, sim::data_type pt);
    ~Serialization();

    void insert(int j, sim::data_type rx, sim::data_type ry, sim::data_type rz, sim::data_type m);
    bool insertInLeaf(int i, int j, sim::data_type rx, sim::data_type ry, sim::data_type rz, sim::data_type m);
    void subdivide(int current);
    void computeAcceleration(int idx, sim::data_type r[3], sim::data_type a[3], sim::data_type g, sim::data_type theta, size_t&);
};

#endif  // SERIALIZATION_H
