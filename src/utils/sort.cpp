#include <algorithm>
#include <iostream>
#include <vector>
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
        std::qsort(ar, size, sizeof(*ar), comp1);
        sort2(ar, size/3, p/3);
        sort2(ar + size/3, size/3, p/3);
        sort2(ar + 2*size/3, size/3, p/3);
    } else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp1);
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
        std::qsort(ar, size, sizeof(*ar), comp2);
        p_sort(ar, size/3, p/3);
        p_sort(ar + size/3, size/3, p/3);
        p_sort(ar + 2*size/3, size/3, p/3);
    } else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp2);
        p_sort(ar, size/5, p/5);
        p_sort(ar + size/5, size/5, p/5);
        p_sort(ar + 2*size/5, size/5, p/5);
        p_sort(ar + 3*size/5, size/5, p/5);
        p_sort(ar + 4*size/5, size/5, p/5);
    }
    
}

void sort0(sim::data_type (*ar)[8], const size_t size, const size_t p, int start, long* costs, std::vector<int>& cuts);
void sort1(sim::data_type (*ar)[8], const size_t size, const size_t p, int start, long* costs, std::vector<int>& cuts);
void sort2(sim::data_type (*ar)[8], const size_t size, const size_t p, int start, long* costs, std::vector<int>& cuts);

void p_sort(sim::data_type (*ar)[8], const size_t size, const size_t p, int* local_N) {
    long* costs = new long[size];
    std::vector<int> cuts; cuts.reserve(p);
    sort0(ar, size, p, 0, costs, cuts);
    
    cuts.push_back(size);
    std::sort(cuts.begin(), cuts.end());
    local_N[0] = cuts[0];
    for (size_t i = 1; i < p; i++) {
        local_N[i] = cuts[i] - cuts[i-1];
    }
    delete[] costs;
}

void sort0(sim::data_type (*ar)[8], const size_t size, const size_t p, int start, long* costs, std::vector<int>& cuts) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp0);
    costs[0] = ar[0][7];
    for (size_t i = 1; i < size; i++) {
        costs[i] = costs[i-1] + ar[i][7];
    }

    if (p % 2 == 0) {
        int cut = size / 2;
        long w = costs[size-1] / 2;
        if (w != 0) cut = std::lower_bound(costs, costs+size, w) - costs;
        cuts.push_back(start+cut);

        sort1(ar, cut, p/2, start, costs, cuts);
        sort1(ar + cut, size - cut, p/2, start+cut, costs+cut, cuts);
    } else if (p % 3 == 0) {
        int cut1 = size / 3;
        int cut2 = 2*size / 3;
        long w = costs[size-1] / 3;
        if (w != 0) {
            cut1 = std::lower_bound(costs, costs+size, w) - costs;
            cut2 = std::lower_bound(costs, costs+size, 2*w) - costs;
        }
        cuts.push_back(start+cut1);
        cuts.push_back(start+cut2);

        sort1(ar, cut1, p/3, start, costs, cuts);
        sort1(ar + cut1, cut2 - cut1, p/3, start+cut1, costs+cut1, cuts);
        sort1(ar + cut2, size - cut2, p/3, start+cut2, costs+cut2, cuts);
    } /*else if (p % 5 == 0) {
        std::qsort(ar, size, sizeof(*ar), comp0);
        sort1(ar, size/5, p/5);
        sort1(ar + size/5, size/5, p/5);
        sort1(ar + 2*size/5, size/5, p/5);
        sort1(ar + 3*size/5, size/5, p/5);
        sort1(ar + 4*size/5, size/5, p/5);
    }*/
}

void sort1(sim::data_type (*ar)[8], const size_t size, const size_t p, int start, long* costs, std::vector<int>& cuts) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp1);
    costs[0] = ar[0][7];
    for (size_t i = 1; i < size; i++) {
        costs[i] = costs[i-1] + ar[i][7];
    }

    if (p % 2 == 0) {
        int cut = size / 2;
        long w = costs[size-1] / 2;
        if (w != 0) cut = std::lower_bound(costs, costs+size, w) - costs;
        cuts.push_back(start+cut);

        sort2(ar, cut, p/2, start, costs, cuts);
        sort2(ar + cut, size - cut, p/2, start+cut, costs+cut, cuts);
    } else if (p % 3 == 0) {
        int cut1 = size / 3;
        int cut2 = 2*size / 3;
        long w = costs[size-1] / 3;
        if (w != 0) {
            cut1 = std::lower_bound(costs, costs+size, w) - costs;
            cut2 = std::lower_bound(costs, costs+size, 2*w) - costs;
        }
        cuts.push_back(start+cut1);
        cuts.push_back(start+cut2);

        sort2(ar, cut1, p/3, start, costs, cuts);
        sort2(ar + cut1, cut2 - cut1, p/3, start+cut1, costs+cut1, cuts);
        sort2(ar + cut2, size - cut2, p/3, start+cut2, costs+cut2, cuts);
    } /*else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp1);
        sort2(ar, size/5, p/5);
        sort2(ar + size/5, size/5, p/5);
        sort2(ar + 2*size/5, size/5, p/5);
        sort2(ar + 3*size/5, size/5, p/5);
        sort2(ar + 4*size/5, size/5, p/5);
    }*/
    
}

void sort2(sim::data_type (*ar)[8], const size_t size, const size_t p, int start, long* costs, std::vector<int>& cuts) {
    if (p == 1) return;

    std::qsort(ar, size, sizeof(*ar), comp1);
    costs[0] = ar[0][7];
    for (size_t i = 1; i < size; i++) {
        costs[i] = costs[i-1] + ar[i][7];
    }

    if (p % 2 == 0) {
        int cut = size / 2;
        long w = costs[size-1] / 2;
        if (w != 0) cut = std::lower_bound(costs, costs+size, w) - costs;
        cuts.push_back(start+cut);

        sort2(ar, cut, p/2, start, costs, cuts);
        sort2(ar + cut, size - cut, p/2, start+cut, costs+cut, cuts);
    } else if (p % 3 == 0) {
        int cut1 = size / 3;
        int cut2 = 2*size / 3;
        long w = costs[size-1] / 3;
        if (w != 0) {
            cut1 = std::lower_bound(costs, costs+size, w) - costs;
            cut2 = std::lower_bound(costs, costs+size, 2*w) - costs;
        }
        cuts.push_back(start+cut1);
        cuts.push_back(start+cut2);

        sort2(ar, cut1, p/3, start, costs, cuts);
        sort2(ar + cut1, cut2 - cut1, p/3, start+cut1, costs+cut1, cuts);
        sort2(ar + cut2, size - cut2, p/3, start+cut2, costs+cut2, cuts);
    } /*else if ( p % 5 ==0){
        std::qsort(ar, size, sizeof(*ar), comp2);
        p_sort0(ar, size/5, p/5);
        p_sort0(ar + size/5, size/5, p/5);
        p_sort0(ar + 2*size/5, size/5, p/5);
        p_sort0(ar + 3*size/5, size/5, p/5);
        p_sort0(ar + 4*size/5, size/5, p/5);
    }*/
    
}
