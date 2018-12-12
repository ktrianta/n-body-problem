#ifndef SORT_HPP
#define SORT_HPP

#include <algorithm>
#include "types.hpp"

void p_sort(sim::data_type (*ar)[7], const size_t size, const size_t p);
void p_sort(sim::data_type (*ar)[8], const size_t size, const size_t p, int* local_N);

#endif  // SORT_HPP
