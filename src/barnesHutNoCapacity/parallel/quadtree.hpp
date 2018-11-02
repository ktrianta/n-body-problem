#ifndef QUADTREE_H
#define QUADTREE_H

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

class QuadTree{

    static const int dimension = 3;
    double centerMass[3];
    double cum_Mass;
   

    Rectangle boundary;
    int index;
    double (*data)[3];
    vector<double>  mass;
    
//   Parent node
    QuadTree* parent;
    bool leaf;
    int size;
    int cum_size;

// Children
    QuadTree* fnorthWest;
    QuadTree* fnorthEast;
    QuadTree* fsouthWest;
    QuadTree* fsouthEast;
    QuadTree* bnorthWest;
    QuadTree* bnorthEast;
    QuadTree* bsouthWest;
    QuadTree* bsouthEast;
    

// private functions
public:
    QuadTree(double (*r)[3], vector<double>& m, int N, double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht);
    QuadTree(QuadTree* inp_parent, double (*r)[3], vector<double>& m, int N, double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht);
    QuadTree(QuadTree* inp_parent, double (*r)[3], vector<double>& m,  double inp_x, double inp_y, double inp_z, double inp_hw, double inp_hh, double inp_ht);
    ~QuadTree();
    bool insert(int new_index);
    void subdivide();
    void print();
    void computeAcceleration(int idxGlobal, int idxLocal, double (*)[3], double (*)[3], double g, double theta);
private:
    void initialization(QuadTree* inp_parent, double (*r)[3], vector<double>& m, double inp_x, double inp_y, double inp_z, double inp_w, double inp_h, double inp_t);
    void fill(int N);
};




#endif  // QUADTREE_H
