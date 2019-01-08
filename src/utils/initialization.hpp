#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "types.hpp"

void initializePositionOnUnitSquare(int N, sim::data_type (*r)[3]);
void initializePositionOnSphere(int N, sim::data_type (*r)[8]);
void initializePositionOnSphere(int N, sim::data_type (*r)[7]);
void initializePositionOnSphere(int N, sim::data_type (*r)[3], sim::data_type *m, sim::data_type (*u)[3]);
void initializePositionOnSphere(int N, sim::data_type *rx, sim::data_type *ry, sim::data_type *rz, sim::data_type *m, sim::data_type (*u)[3]);

#endif  // INITIALIZATION_H
