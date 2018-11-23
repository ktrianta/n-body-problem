#include <algorithm>
#include "sort.hpp"

int comp0(const void *arg1, const void *arg2) {
    sim::data_type const *lhs = static_cast<sim::data_type const*>(arg1);
    sim::data_type const *rhs = static_cast<sim::data_type const*>(arg2);
    return (lhs[1] < rhs[1]) ? -1 : ((lhs[1] > rhs[1]) ? 1 : 0);
}

int comp1(const void *arg1, const void *arg2) {
    sim::data_type const *lhs = static_cast<sim::data_type const*>(arg1);
    sim::data_type const *rhs = static_cast<sim::data_type const*>(arg2);
    return (lhs[2] < rhs[2]) ? -1 : ((lhs[2] > rhs[2]) ? 1 : 0);
}

int comp2(const void *arg1, const void *arg2) {
    sim::data_type const *lhs = static_cast<sim::data_type const*>(arg1);
    sim::data_type const *rhs = static_cast<sim::data_type const*>(arg2);
    return (lhs[3] < rhs[3]) ? -1 : ((lhs[3] > rhs[3]) ? 1 : 0);
}

void sort1(sim::data_type (*ar)[7], const size_t size, const size_t p);
void sort2(sim::data_type (*ar)[7], const size_t size, const size_t p);

void p_sort(sim::data_type (*ar)[7], const size_t size, const size_t p) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp0);
    sort1(ar, size/2, p/2);
    sort1(ar + size/2, size/2, p/2);
}

void sort1(sim::data_type (*ar)[7], const size_t size, const size_t p) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp1);
    sort2(ar, size/2, p/2);
    sort2(ar + size/2, size/2, p/2);
}

void sort2(sim::data_type (*ar)[7], const size_t size, const size_t p) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp2);
    p_sort(ar, size/2, p/2);
    p_sort(ar + size/2, size/2, p/2);
}
