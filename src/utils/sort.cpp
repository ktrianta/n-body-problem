#include <algorithm>
#include <iostream>
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
    
    if (p % 2 == 0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort1(ar, size/2, p/2);
        sort1(ar + size/2, size/2, p/2);
    } else if ( p % 3 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort1(ar, size/3, p/3);
        sort1(ar + size/3, size/3, p/3);
        sort1(ar + 2*size/3, size/3, p/3);
    } else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort1(ar, size/5, p/5);
        sort1(ar + size/5, size/5, p/5);
        sort1(ar + 2*size/5, size/5, p/5);
        sort1(ar + 3*size/5, size/5, p/5);
        sort1(ar + 4*size/5, size/5, p/5);
    }
}

void sort1(sim::data_type (*ar)[7], const size_t size, const size_t p) {
    if (p == 1) return;

    if (p % 2 == 0)
    {
        std::qsort(ar, size, sizeof(*ar), comp1);
        sort2(ar, size/2, p/2);
        sort2(ar + size/2, size/2, p/2);
    } else if ( p % 3 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort2(ar, size/3, p/3);
        sort2(ar + size/3, size/3, p/3);
        sort2(ar + 2*size/3, size/3, p/3);
    } else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort2(ar, size/5, p/5);
        sort2(ar + size/5, size/5, p/5);
        sort2(ar + 2*size/5, size/5, p/5);
        sort2(ar + 3*size/5, size/5, p/5);
        sort2(ar + 4*size/5, size/5, p/5);
    }
    
}

void sort2(sim::data_type (*ar)[7], const size_t size, const size_t p) {
    if (p == 1) return;

    if (p % 2 == 0)
    {
        std::qsort(ar, size, sizeof(*ar), comp2);
        p_sort(ar, size/2, p/2);
        p_sort(ar + size/2, size/2, p/2);
    } else if ( p % 3 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        p_sort(ar, size/3, p/3);
        p_sort(ar + size/3, size/3, p/3);
        p_sort(ar + 2*size/3, size/3, p/3);
    } else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        p_sort(ar, size/5, p/5);
        p_sort(ar + size/5, size/5, p/5);
        p_sort(ar + 2*size/5, size/5, p/5);
        p_sort(ar + 3*size/5, size/5, p/5);
        p_sort(ar + 4*size/5, size/5, p/5);
    }
    
}


void sort0(sim::data_type (*ar)[8], const size_t size, const size_t p, int idx, int* local_N);
void sort1(sim::data_type (*ar)[8], const size_t size, const size_t p, int idx, int* local_N);
void sort2(sim::data_type (*ar)[8], const size_t size, const size_t p, int idx, int* local_N);

void p_sort(sim::data_type (*ar)[8], const size_t size, const size_t p, int* local_N) {
    sort0(ar, size, p, 0, local_N);
}

void sort0(sim::data_type (*ar)[8], const size_t size, const size_t p, int idx, int* local_N) {
    if (p == 1) {
        local_N[idx] = size;
        return;
    }
  
    if (p % 2 == 0) {
        std::qsort(ar, size, sizeof(*ar), comp0);

        long* costs = new long[size]; costs[0] = (long) ar[0][7];
        for (size_t i = 1; i < size; i++) {
            costs[i] = costs[i-1] + ar[i][7];
        }

        long w2 = costs[size-1] / 2;
        size_t cut = std::lower_bound(costs, costs+size, w2) - costs;
        if (w2 == 0) cut = size/2;

        sort1(ar, cut, p/2, idx*2, local_N);
        sort1(ar + cut, size - cut, p/2, idx*2+1, local_N);
    }/* else if ( p % 3 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort1(ar, size/3, p/3);
        sort1(ar + size/3, size/3, p/3);
        sort1(ar + 2*size/3, size/3, p/3);
    } else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort1(ar, size/5, p/5);
        sort1(ar + size/5, size/5, p/5);
        sort1(ar + 2*size/5, size/5, p/5);
        sort1(ar + 3*size/5, size/5, p/5);
        sort1(ar + 4*size/5, size/5, p/5);
    }*/
}

void sort1(sim::data_type (*ar)[8], const size_t size, const size_t p, int idx, int* local_N) {
    if (p == 1) {
        local_N[idx] = size;
        return;
    }

    if (p % 2 == 0) {
        std::qsort(ar, size, sizeof(*ar), comp0);

        double* costs = new double[size]; costs[0] = ar[0][7];
        for (size_t i = 1; i < size; i++) {
            costs[i] = costs[i-1] + ar[i][7];
        }

        double w2 = costs[size-1] / 2;
        size_t cut = std::lower_bound(costs, costs+size, w2) - costs;
        if (w2 == 0) cut = size/2;

        sort2(ar, cut, p/2, idx*2, local_N);
        sort2(ar + cut, size - cut, p/2, idx*2+1, local_N);
    } /*else if ( p % 3 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort2(ar, size/3, p/3);
        sort2(ar + size/3, size/3, p/3);
        sort2(ar + 2*size/3, size/3, p/3);
    } else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort2(ar, size/5, p/5);
        sort2(ar + size/5, size/5, p/5);
        sort2(ar + 2*size/5, size/5, p/5);
        sort2(ar + 3*size/5, size/5, p/5);
        sort2(ar + 4*size/5, size/5, p/5);
    }*/
    
}

void sort2(sim::data_type (*ar)[8], const size_t size, const size_t p, int idx, int* local_N) {
    if (p == 1) {
        local_N[idx] = size;
        return;
    }

    if (p % 2 == 0) {
        std::qsort(ar, size, sizeof(*ar), comp0);

        long* costs = new long[size]; costs[0] = (long) ar[0][7];
        for (size_t i = 1; i < size; i++) {
            costs[i] = costs[i-1] + ar[i][7];
        }

        long w2 = costs[size-1] / 2;
        size_t cut = std::lower_bound(costs, costs+size, w2) - costs;
        if (w2 == 0) cut = size/2;

        sort0(ar, cut, p/2, idx*2, local_N);
        sort0(ar + cut, size - cut, p/2, idx*2+1, local_N);
    } /*else if ( p % 3 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        p_sort0(ar, size/3, p/3);
        p_sort0(ar + size/3, size/3, p/3);
        p_sort0(ar + 2*size/3, size/3, p/3);
    } else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp0);
        p_sort0(ar, size/5, p/5);
        p_sort0(ar + size/5, size/5, p/5);
        p_sort0(ar + 2*size/5, size/5, p/5);
        p_sort0(ar + 3*size/5, size/5, p/5);
        p_sort0(ar + 4*size/5, size/5, p/5);
    }*/
    
}
