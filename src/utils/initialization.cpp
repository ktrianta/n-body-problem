#include <cmath>
#include <ctime>
#include <cstdlib>
#include "initialization.hpp"

void initializePositionOnUnitSquare(int N, sim_data_type (*r)[2])
{
    std::srand(std::time(NULL));

    for (int i = 0; i < N; i++)
    {
        r[i][0] = ((sim_data_type) std::rand()) / RAND_MAX;
        r[i][1] = ((sim_data_type) std::rand()) / RAND_MAX;
    }
}

sim_data_type randomNum(float a, float b)
{
    sim_data_type num = ((float) std::rand()) / RAND_MAX;
    return a + (b-a) * num;
}

/* This codes initiliazes the particules in a circular shape */
void initializePositionOnSphere(int N, sim_data_type (*r)[2])
{
    std::srand(std::time(NULL));

    for (int i = 0; i < N ; i++)
    {
        sim_data_type p = pow(randomNum(0, 1), 1.0 / 3.0);
        sim_data_type theta = randomNum(0,2 * M_PI);
        r[i][0] = p * cos(theta);
        r[i][1] = p * sin(theta);
    }
}
