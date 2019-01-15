#ifndef SERIALIZATION_H
#define SERIALIZATION_H

#include "types.hpp"

using namespace std;


struct Treenode {
double x;
double y;
double z;
double w;
double h;
double t;
double massCenter[3];
double mass;
int index;
bool leaf;
size_t cum_size;
size_t child;

Treenode();
void set(double px, double py, double pz, double pw, double ph, double pt);
};



class Serialization{
public: 
    size_t size;
    size_t position;
    Treenode *treeArray;

    Serialization(double px, double py, double pz, double pw, double ph, double pt);
    Serialization();
    ~Serialization();

    void insert(int j, double rx, double ry, double rz, double m);
    bool insertInLeaf(int i, int j, double rx, double ry, double rz, double m);
    void subdivide(int current);
    void computeAcceleration(int current, int idx, double (*r)[3], double (*a)[3], double g, double theta);
};

#endif  // SERIALIZATION_H
