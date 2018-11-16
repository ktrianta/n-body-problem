#ifndef BOXCOMPUTATION_H
#define BOXCOMPUTATION_H

#include "types.hpp"

void boxComputation(const size_t N, sim::data_type (*r)[3], sim::data_type &x, sim::data_type &y, sim::data_type &z, sim::data_type &w, sim::data_type &h, sim::data_type &t);

#endif  // BOXCOMPUTATION_H
