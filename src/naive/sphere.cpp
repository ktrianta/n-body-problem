#include<iostream>
#include<cmath>
#include <cstdlib>

/* This codes initiliazes the particules in a spherical shape */

int N = 5; /* Number of particules */
double r[N][2];
double u[N][2];
double m[N];

double randomNum(double a, double b)
{
	return rand()%(b-a + 1) + a;
}

void initiliazeOnSphere(double r[][2],double u[][2],int N,double m[])
{
	for (int k=0;k<N;k++)
	{
		m[k] = 1.0/(double)(N);
		u[k][0] = 0;
		u[k][1] = 0;
	}

	double r = pow(randomNum(0,1),1.0/3.0);
	double theta = randomNum(0,2*M_PI);

	for (int k=0;k<N;k++)
	{
		r[k][0] = r * cos(theta);
		r[k][1] = r * sin(theta);
	}

}