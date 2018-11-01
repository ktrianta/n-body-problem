#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include <vector>
#include "array2d.hpp"
#include "types.hpp"

void initializePositionOnUnitSquare(int N, sim_data_type (*r)[3]);
void initializePositionOnSphere(int N, sim_data_type (*r)[3]);

#endif  // INITIALIZATION_H
