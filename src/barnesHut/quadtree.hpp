#ifndef QUADTREE_H
#define QUADTREE_H

#include"../utils/array2d.hpp"

class Rectangle{
    public:
        double x; // the x-center of the rectangle
        double y; // the y-center of the rectangle
        double h; // the height/2 of the rectangle
        double w; // the width/2  of the rectanglea
        bool containPoint(double point[]);
};


class QuadTree{

    static const int capacity = 1;
    static const int dimension = 2;
    double centerMass[dimension];
   

    Rectangle boundary;
    int index[capacity];
    double (*data)[2];
    
// Parent node
    QuadTree* parent;
    bool leaf;
    int size;
    int cum_size;

// Children
    QuadTree* northWest;
    QuadTree* northEast;
    QuadTree* southWest;
    QuadTree* southEast;
    

// private functions
public:
    QuadTree(double (*r)[2], int N, double inp_x, double inp_y, double inp_hw, double inp_hh);
    QuadTree(QuadTree* inp_parent, double (*r)[2], int N, double inp_x, double inp_y, double inp_hw, double inp_hh);
    QuadTree(QuadTree* inp_parent, double (*r)[2], double inp_x, double inp_y, double inp_hw, double inp_hh);
    ~QuadTree();
    bool insert(int new_index);
    void subdivide();

private:
    void initialization(QuadTree* inp_parent, double (*r)[2], double inp_x, double inp_y, double inp_w, double inp_h);
    void fill(int N);
};



// daclare functions

#endif  // QUADTREE_H
