#include <iostream>
#include <math.h>

#define PI acos(-1)
#define GAMMA 0.5772156649;

double first_bessel_0(double x, int n)
{
    double sum;
    for(int m = 1; m <= pow(2, n) - 1; sum += cos(x * sin(m * PI * pow(2, -n))), m++);
    return pow(2, -n) * sum;
}

// Blows up to inf and -inf
double second_bessel_0(double x, int n)
{
    if(x == 0)
        return 0;
    double sum = 0;
    for(int m = 1; m <= pow(2, n) - 1; m++)
    {
        double p1 = cos(x * cos(m * PI * pow(2, (-n) - 1)));
        double p2 = log(2.0 * x * pow(sin(m * PI * pow(2, (-n) - 1)), 2));
        p2 += GAMMA;
        sum += p1 * p2;
    }
    return (2 / PI) * pow(2, -n) * sum;
}