#include <cmath>
#include <ctime>
#include <cstdlib>
#include "serialization.hpp"
#include <stdio.h> //for printf
#include <queue>

Treenode::Treenode()
    : x(0), y(0), z(0), w(0), h(0), t(0), index(-1), leaf(true), cum_size(0), child(0)
    {}

    
void Treenode::set(double px, double py, double pz, double pw, double ph, double pt){
     x = px;
     y = py;
     z = pz;
     w = pw;
     h = ph;
     t = pt;
     index = -1;
     leaf = true;
     cum_size = 0;
     child = 0;
     }


Serialization::Serialization(double px, double py, double pz, double pw, double ph, double pt)
    : treeArray(new Treenode[64])
    {
    treeArray[0].x = px;
    treeArray[0].y = py;
    treeArray[0].z = pz;
    treeArray[0].w = pw;
    treeArray[0].h = ph;
    treeArray[0].t = pt;
    position = 0;
    current = 1;
    }
 
//functions

bool Serialization::insert(int j, double rx, double ry, double rz, double m){

    queue<int> list;
    list.push (0);

    while(!list.empty())
    {
        int f = list.front();
        list.pop();
        const double x = treeArray[f].x;
        const double y = treeArray[f].y;
        const double z = treeArray[f].z;
        const double w = treeArray[f].w;
        const double h = treeArray[f].h;
        const double t = treeArray[f].t;
        


        if (x - w > rx) return false;
        if (x + w < rx) return false;
        if (y - h > ry) return false;
        if (y + h < ry) return false;
        if (z - t > rz) return false;
        if (z + t < rz) return false;

        size_t cum_size = ++treeArray[f].cum_size;
        double mult1 = (double) (cum_size - 1) / (double) treeArray[f].cum_size;                                   
        double mult2 = 1.0 / (double) cum_size;                                                       
        for (int d = 0; d < 3; d++) treeArray[f].massCenter[d] *= mult1;                                   
        treeArray[f].massCenter[0] += mult2 * x;
        treeArray[f].massCenter[1] += mult2 * y;
        treeArray[f].massCenter[2] += mult2 * z;
        treeArray[f].mass += m; 

        if (treeArray[f].leaf && treeArray[f].cum_size == 1)
        {
            treeArray[f].index = j;
            return true;
        }

        if (treeArray[f].leaf)
        { 
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
    return true;
}

void Serialization::subdivide(int current){

    if(position + 8 > size){
      size *= 2;
      Treenode *temp = new Treenode[size];
      std::copy(treeArray, treeArray +size/2, temp);
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
    
    treeArray[position].leaf = false;
    treeArray[position].index = -1;


    position += 8;
}
void Serialization::computeAcceleration(int current, int idx, double (*r)[3], double (*a)[3], double g, double theta) {
    if (treeArray[current].cum_size == 0 || (treeArray[current].leaf == true && treeArray[current].cum_size == 1 && treeArray[current].index == idx)) {                           
        return;                                                                                   
    }                                                                                             
                                                                                                  
    double  centerMass[3];
    centerMass[0] = treeArray[current].massCenter[0];
    centerMass[1] = treeArray[current].massCenter[1];
    centerMass[2] = treeArray[current].massCenter[2];

    double d = sqrt( (centerMass[0] - r[idx][0]) * (centerMass[0] - r[idx][0])                    
                   + (centerMass[1] - r[idx][1]) * (centerMass[1] - r[idx][1])
                   + (centerMass[2] - r[idx][2]) * (centerMass[2] - r[idx][2]));                  
    
    if (2.*treeArray[current].w/d <= theta) {                                                               
        
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
        this->computeAcceleration(current, idx, r, a, g, theta); 
        current ++;
        this->computeAcceleration(current, idx, r, a, g, theta); 
        current ++;
        this->computeAcceleration(current, idx, r, a, g, theta); 
        current ++;
        this->computeAcceleration(current, idx, r, a, g, theta); 
        current ++;
        this->computeAcceleration(current, idx, r, a, g, theta); 
        current ++;
        this->computeAcceleration(current, idx, r, a, g, theta); 
        current ++;
        this->computeAcceleration(current, idx, r, a, g, theta); 
        current ++;
        this->computeAcceleration(current, idx, r, a, g, theta); 
    }                                                                                             
}    




