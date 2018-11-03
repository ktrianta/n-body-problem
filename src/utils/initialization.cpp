#include <cmath>
#include <ctime>
#include <cstdlib>
#include "initialization.hpp"

void initializePositionOnUnitSquare(int N, sim::data_type (*r)[3])
{
    std::srand(std::time(NULL));

    for (int i = 0; i < N; i++)
    {
        r[i][0] = ((sim::data_type) std::rand()) / RAND_MAX;
        r[i][1] = ((sim::data_type) std::rand()) / RAND_MAX;
        r[i][2] = ((sim::data_type) std::rand()) / RAND_MAX;
    }
}

sim::data_type randomNum(float a, float b)
{
    sim::data_type num = ((float) std::rand()) / RAND_MAX;
    return a + (b-a) * num;
}

/* This codes initiliazes the particules in a circular shape */
void initializePositionOnSphere(int N, sim::data_type (*r)[3])
{
    std::srand(std::time(NULL));

    for (int i = 0; i < N ; i++)
    {
        sim::data_type p = pow(randomNum(0, 1), 1.0 / 3.0);
        sim::data_type theta = randomNum(0, M_PI);
        sim::data_type phi = randomNum(0, 2 * M_PI);
        r[i][0] = p * sin(theta) * cos(phi);
        r[i][1] = p * sin(theta) * sin(phi);
        r[i][2] = p * cos(theta);
    }
}
