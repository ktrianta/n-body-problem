#include <math.h> //for sqrt function
#include <stdio.h> //for printf
#include <stdlib.h> //for srand/rand prng
#include <time.h> //for seeding prng
#include <iostream>
#include <fstream>

using namespace std;

const int N=5;      // the number of particles
const float g=1 ;     // gravitational constant
const float epsilon=0.001;

// min ksexaseis na arxikopoihseis mazes kai taxitites

void InitialPositions(float *x,float *y)
{
    srand(time(NULL));
    for (int i=0; i<N; i++)
    {
        x[i] = ((double) rand()) / RAND_MAX;
        y[i] = ((double) rand()) / RAND_MAX;
    }
}
void ComputeF(int j,float *x,float *y, float *m,float *ax)
{
    for (int i=j+1; i<N; i++)
    {
        if(i!=j)
        {
            float d=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
                float denom = sqrt((d*d+epsilon*epsilon)*(d*d+epsilon*epsilon)*(d*d+epsilon*epsilon));
                float force = -g*m[j]*m[i]*(x[j]-x[i])/denom;
                ax[j]+=force/m[j];
                ax[i]-=force/m[i];
              //  sumF+=F;
        }
    }
}


int main(int argc,char** argv)
{
    float x[N],y[N],ux[N],uy[N],m[N],ax[N],ay[N];
    ofstream myfile;
    myfile.open ("output.dat");

    InitialPositions(x,y);
    for (int i; i<N; i++)
    { 
        ux[i]=0;
        uy[i]=0;
        ax[i]=0;
        ay[i]=0;
        m[i]=1;
        myfile <<  x[i] << "   " 
               <<  y[i] << "   " 
               << ux[i] << "   " 
               << uy[i] << "\n";
    }

    const float T=10;
    const float dt=0.00001; 
    const int Ntimesteps=T/dt+1;
    
    for (int k=0; k<N; k++)
    {
          ComputeF(k,x,y,m,ax);
          ComputeF(k,y,x,m,ay);
    }

    for (int t=0; t<Ntimesteps; t++)
    {
        cout << "time=" <<  t*dt << "\n";
        for (int j=0; j<N; j++)
        {
            ux[j]+=0.5*ax[j]*dt;
            uy[j]+=0.5*ay[j]*dt;
//          compute new x position        
            x[j]+=ux[j]*dt;
            y[j]+=uy[j]*dt;
        }
        for (int j=0; j<N; j++)
        {
            ax[j]=0.;
            ay[j]=0.;
        }
        for (int j=0; j<N; j++)
        {
//          compute the new accelaration ~ force which is required for the velocity
            ComputeF(j,x,y,m,ax);
            ComputeF(j,y,x,m,ay);
//          compute the new velocity
            ux[j]+=0.5*ax[j]*dt;
            uy[j]+=0.5*ay[j]*dt;
            if(t%200==0)
            {
                myfile <<  x[j] << "   " 
                       <<  y[j] << "   " 
                       << ux[j] << "   " 
                       << uy[j] << "\n";
            }
        }
    }
    return 0;
}
