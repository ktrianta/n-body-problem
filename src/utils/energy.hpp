#ifndef ENERGY_H
#define ENERGY_H

#include "types.hpp"
#include <iostream>
#include <math.h>

sim::data_type energy(const size_t N, sim::data_type (*r)[7]);
sim::data_type energy(const size_t N, sim::data_type (*r)[8]);
sim::data_type energy(const size_t N, sim::data_type (*r)[3], sim::data_type (*u)[3], sim::data_type *m);
sim::data_type energy(const size_t N, sim::data_type *rx, sim::data_type *ry, sim::data_type *rz, sim::data_type (*u)[3], sim::data_type *m);
void printEnergy(sim::data_type a, sim::data_type b);

#endif  // ENERGY_H
