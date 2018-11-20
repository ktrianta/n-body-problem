#ifndef OCTREE_H
#define OCTREE_H

#include "types.hpp"

class Rectangle {
    public:
        sim::data_type x; // the x-center of the rectangle
        sim::data_type y; // the y-center of the rectangle
        sim::data_type z; // the z-center of the rectangle
        sim::data_type w; // the width/2  of the rectangle
        sim::data_type h; // the height/2 of the rectangle
        sim::data_type t; // the depth/2  of the rectangle

        Rectangle(sim::data_type px, sim::data_type py, sim::data_type pz,
                  sim::data_type pw, sim::data_type ph, sim::data_type pt);
        bool containPoint(const sim::data_type point[]);
};

class Octree {
public:
    // Static fields
    // Need to be initialized by caller
    static sim::data_type (*r)[3];
    static sim::data_type (*a)[3];
    static sim::data_type *m;
    static sim::data_type g;
    static sim::data_type theta;

    // Methods
    Octree(size_t N, sim::data_type px, sim::data_type py, sim::data_type pz, sim::data_type pw,
           sim::data_type ph, sim::data_type pt);
    Octree(sim::data_type px, sim::data_type py, sim::data_type pz, sim::data_type pw,
           sim::data_type ph, sim::data_type pt);

    ~Octree();

    void print();
    void computeAcceleration(int idx);
private:
    static const size_t dimension = 3;

    Rectangle boundary;
    sim::data_type massCenter[dimension];
    sim::data_type mass;
    int index;
    bool leaf;
    size_t cumSize;

    // Children
    Octree* fnorthWest;
    Octree* fnorthEast;
    Octree* fsouthWest;
    Octree* fsouthEast;
    Octree* bnorthWest;
    Octree* bnorthEast;
    Octree* bsouthWest;
    Octree* bsouthEast;

    // Methods
    void fill(size_t N);
    bool insert(int new_index);
    void subdivide();
};

#endif  // OCTREE_H
