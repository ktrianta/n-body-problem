#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include <vector>
#include "array2d.hpp"
#include "types.hpp"

void initializePositionOnUnitSquare(int N, Array2D<sim_data_type>& r);
void initializePositionOnSphere(int N, Array2D<sim_data_type>& r);

#endif  // INITIALIZATION_H
