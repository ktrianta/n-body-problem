#include <cmath>
#include <ctime>
#include <cstdlib>

void initializePositionOnUnitSquare(int N, float r[][2])
{
    std::srand(std::time(NULL));

    for (int i = 0; i < N; i++)
    {
        r[i][0] = ((float) std::rand()) / RAND_MAX;
        r[i][1] = ((float) std::rand()) / RAND_MAX;
    }
}

float randomNum(float a, float b)
{
    float num = ((float) std::rand()) / RAND_MAX;
    return a + (b-a) * num;
}

/* This codes initiliazes the particules in a circular shape */
void initializePositionOnSphere(int N, float r[][2])
{
    std::srand(std::time(NULL));

    for (int k = 0; k < N ; k++)
    {
        float p = pow(randomNum(0, 1), 1.0 / 3.0);
        float theta = randomNum(0,2 * M_PI);
        r[k][0] = p * cos(theta);
        r[k][1] = p * sin(theta);
    }
}
