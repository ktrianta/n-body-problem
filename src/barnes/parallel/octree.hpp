#ifndef OCTREE_H
#define OCTREE_H

#include <vector>

class Rectangle{
    public:
        double x; // the x-center of the rectangle
        double y; // the y-center of the rectangle
        double z; // the z-center of the rectangle
        double h; // the height/2 of the rectangle
        double w; // the width/2  of the rectangle
        double t; // the 3rd dimension/2  of the rectangle
        bool containPoint(double point[]);
};


using namespace std;

class Octree{

    static const int dimension = 3;
    double centerMass[3];
    double cum_Mass;
   

    Rectangle boundary;
    int index;
    double (*data)[3];
    vector<double>  mass;
    
//   Parent node
    Octree* parent;
    bool leaf;
    int size;
    int cum_size;

// Children
    Octree* fnorthWest;
    Octree* fnorthEast;
    Octree* fsouthWest;
    Octree* fsouthEast;
    Octree* bnorthWest;
    Octree* bnorthEast;
    Octree* bsouthWest;
    Octree* bsouthEast;
    

// private functions
public:
    Octree(double (*r)[3], vector<double>& m, int N, double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht);
    Octree(Octree* inp_parent, double (*r)[3], vector<double>& m, int N, double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht);
    Octree(Octree* inp_parent, double (*r)[3], vector<double>& m,  double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht);
    ~Octree();
    bool insert(int new_index);
    void subdivide();
    void print();
    void computeAcceleration(int idxGlobal, int idxLocal, double (*)[3], double (*)[3], double g, double theta);
private:
    void initialization(Octree* inp_parent, double (*r)[3], vector<double>& m, double inp_x, double inp_y, double inp_z, double inp_w, double inp_h, double inp_t);
    void fill(int N);
};




#endif  // OCTREE_H
