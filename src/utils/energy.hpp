#ifndef ENERGY_H
#define ENERGY_H

#include "types.hpp"
#include <iostream>
#include <math.h>

sim::data_type energy(const size_t N, sim::data_type (*r)[7]);
void printEnergy(sim::data_type a, sim::data_type b);

#endif  // ENERGY_H
