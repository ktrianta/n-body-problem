#include <algorithm>
#include "sort.hpp"

int comp0(const void *arg1, const void *arg2) {
    int const *lhs = static_cast<int const*>(arg1);
    int const *rhs = static_cast<int const*>(arg2);
    return (lhs[0] < rhs[0]) ? -1 : ((lhs[0] > rhs[0]) ? 1 : 0);
}

int comp1(const void *arg1, const void *arg2) {
    int const *lhs = static_cast<int const*>(arg1);
    int const *rhs = static_cast<int const*>(arg2);
    return (lhs[1] < rhs[1]) ? -1 : ((lhs[1] > rhs[1]) ? 1 : 0);
}

int comp2(const void *arg1, const void *arg2) {
    int const *lhs = static_cast<int const*>(arg1);
    int const *rhs = static_cast<int const*>(arg2);
    return (lhs[2] < rhs[2]) ? -1 : ((lhs[2] > rhs[2]) ? 1 : 0);
}

void sort1(int (*ar)[3], const size_t size, const size_t p);
void sort2(int (*ar)[3], const size_t size, const size_t p);

void p_sort(int (*ar)[3], const size_t size, const size_t p) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp0);
    sort1(ar, size/2, p/2);
    sort1(ar + size/2, size/2, p/2);
}

void sort1(int (*ar)[3], const size_t size, const size_t p) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp1);
    sort2(ar, size/2, p/2);
    sort2(ar + size/2, size/2, p/2);
}

void sort2(int (*ar)[3], const size_t size, const size_t p) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp2);
    p_sort(ar, size/2, p/2);
    p_sort(ar + size/2, size/2, p/2);
}
