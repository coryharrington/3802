#include <iostream>
#include <math.h>

#define PI acos(-1)
#define GAMMA 0.5772156649;

const double outer = 0.075;
const double inner = 0.325;

double first_bessel_0(double x, int n)
{
    double sum = 0;
    for(int m = 1; m <= pow(2, n) - 1; sum += cos(x * sin(m * PI * pow(2, -n))), m++);
    return pow(2, -n) * sum;
}

double first_bessel_prime_0(double x, int n)
{
    double sum = 0;
    for(int m = 1; m <= pow(2, n) - 1; m++)
    {
        double term = sin(m * PI * pow(2, -n));
        sum += (-sin(x * term)) * term;
    }
    return pow(2, -n) * sum;
}

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

double second_bessel_prime_0(double x, int n)
{
    if(x == 0)
        return 0;
    double sum = 0;
    for(int m = 1; m <= pow(2, n) - 1; m++)
    {
        double t1 = cos(m * PI * pow(2, (-n) - 1));
        double t2 = sin(m * PI * pow(2, (-n) - 1));
        double h = log(2.0 * x * pow(sin(m * PI * pow(2, (-n) - 1)), 2));
        h += GAMMA;
        double p1 = -sin(x * t1) * t1 * h;
        double p2 = cos(x * t1) * pow(x * pow(t2 , 2), -1) * pow(t2, 2);
        sum += p1 + p2;
    }
    return (2 / PI) * pow(2, -n) * sum;
}

// J0(ka)Y0(kb) - J0(kb)Y0(ka) = 0
double zero_function(double k, int n)
{
    return
        (first_bessel_0(k*outer, n)*second_bessel_0(k*inner, n)) - \
        (first_bessel_0(k*inner, n)*second_bessel_0(k*outer, n));
}

double zero_function_prime(double k, int n)
{
    return
        (first_bessel_prime_0(k*outer, n)*second_bessel_0(k*inner, n)) + \
        (first_bessel_0(k*outer, n)*(second_bessel_prime_0(k*inner, n))) - \
        (
            (first_bessel_prime_0(k*inner, n)*second_bessel_0(k*outer, n)) + \
            (first_bessel_0(k*inner, n)*second_bessel_prime_0(k*outer, n))
        );
}