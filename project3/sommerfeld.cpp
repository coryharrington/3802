#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>
#include <stdio.h>

#define PI acos(-1)

using namespace std;
using namespace std::complex_literals;

const double light_spd = 2.99792458e8;

// Besel's from method of exhastion (triangles)
double first_bessel_0(double x, int n)
{
    double sum = 0;
    for(int m = 1; m <= pow(2, n) - 1; sum += cos(x * sin(m * PI * pow(2, -n))), m++);
    return pow(2, -n) * sum;
}

// Rhomberg iteration
complex<double> rhomberg_iteration(int n, int m, double a, double b, double k, double r, complex<double>(*function)(double, double, double))
{
    auto h = [a, b](int n){return (b - a) / pow(2, n);};
    double hn = h(n);
    if(n == 0 && m == 0)
        return h(1) * (function(a, k, r) + function(b, k, r));
    else if(m == 0)
    {
        complex<double> sum = 0;
        for(int k = 1; k <= pow(2, n - 1); k++) sum += function(a + hn * (2*k-1), k, r);
        return 0.5 * rhomberg_iteration(n-1, 0, a, b, k, r, function) + hn * sum;
    } else 
    {
        complex<double> r1 = rhomberg_iteration(n, m-1, a, b, k, r, function);
        complex<double> r2 = rhomberg_iteration(n-1, m-1, a, b, k, r, function);
        return r1 + (1 / (pow(4, m) + 1)) * (r1 - r2);
    }  
}

complex<double> rhomberg_iteration_testing(int n, int m, double a, double b, complex<double>(*function)(double))
{
    auto h = [a, b](int n){return (b - a) / pow(2, n);};
    double hn = h(n);
    if(n == 0 && m == 0)
        return h(1) * (function(a) + function(b));
    else if(m == 0)
    {
        complex<double> sum = 0;
        for(int k = 1; k <= pow(2, n - 1); k++) sum += function(a + hn * (2*k-1));
        return 0.5 * rhomberg_iteration_testing(n-1, 0, a, b, function) + hn * sum;
    } else 
    {
        complex<double> r1 = rhomberg_iteration_testing(n, m-1, a, b, function);
        complex<double> r2 = rhomberg_iteration_testing(n-1, m-1, a, b, function);
        return r1 + (1 / (pow(4, m) + 1)) * (r1 - r2);
    }  
}

// Rectangle exhastion
complex<double> method_of_exhastion(int N, double a, double b, double k, double r, complex<double> (*function)(double, double, double))
{
    complex<double> sum = (0, 0);
    for(int n = 1; n <= pow(2, N - 1); n++)
    {
        double point = ((2*n - 1) * b + (pow(2, N) - (2*n - 1)) * a) / pow(2, N);
        sum += function(point, k, r);
    }
    return (b - a) * sum / (pow(2, N-1));
}

// Rectangle exhastion
complex<double> method_of_exhastion_testing(int N, double a, double b, complex<double> (*function)(double))
{
    complex<double> sum = (0, 0);
    for(int n = 1; n <= pow(2, N - 1); n++) sum += function(((2*n - 1) * b + (pow(2, N) - (2*n - 1)) * a) / pow(2, N));
    return (b - a) * sum / (pow(2, N-1));
}

// f => k
double wave_number(double f)
{
    return 2 * PI * f / (light_spd);
}

// Testing
complex<double> testing(double x)
{
    return sin(x);
}

// Transform of bessel => cylindrical coordiantes.
complex<double> sommerfeld_identity_analytic(double k, double r)
{
    return (1 / r) * (cos(k*r) - 1i*sin(k*r));
}

// Analytic solution 
complex<double> sommerfeld_identity_integrand(double a, double k, double r)
{
    if(a == k) return 0i;
    return a < k ? 1i * -a * first_bessel_0(a*r, 10) / sqrt(pow(k, 2) - pow(a, 2)) : a * first_bessel_0(a*r, 10) / sqrt(pow(a, 2) - pow(k, 2));
}

void help_msg()
{
    cout << "Please use ./sommerfeld dx f0 f r0\n";
    exit(0);
}

int main(int argc, char **argv)
{
    // Invalid args
    if(argc != 5 || argv[1] == "help" || argv[1] == "--help")
        help_msg();
    // Go through arguments as double's
    double params[4];
    std::string::size_type sz;
    for(int i = 1; i < argc; i++)
    {
        const char *str = argv[i];
        params[i - 1] = std::stod(str, &sz);
    }
    if(params[0] >= params[2] - params[1])
        cout << "Give a correct step size.\n";
    if(params[1] > params[2])
    {
        cout << "Give a correct range of frequencies, f > 0, f > f0\n";
        return 0;
    }
    complex<double> (*testing_integrand)(double) = testing;
    cout << "Integrating sin(x) from 0 to pi\n";
    cout << method_of_exhastion_testing(20, 0, PI, testing_integrand) << "\n";
    cout << rhomberg_iteration_testing(10, 3, 0, PI, testing_integrand) << "\n";
    complex<double> (*integrand)(double, double, double) = sommerfeld_identity_integrand;
    double terror1, terror2 = 0;
    int win1 = 0;
    int win2 = 0;
    ofstream moe_out, rom_out;
    moe_out.open("moe_out.csv");
    rom_out.open("rom_out.csv");
    for(double f = params[1]; f <= params[2]; f += params[0])
    {
        double k = wave_number(f);
        complex<double> a_result = sommerfeld_identity_analytic(k, params[3]);
        complex<double> moe_result = method_of_exhastion(12, 0, k, k, params[3], integrand) + method_of_exhastion(12, k, 50 * k, k, params[3], integrand);
        complex<double> rom_result = rhomberg_iteration(8, 4, 0, 50 * k, k, params[3], integrand);
        double error1 = std::abs(a_result - moe_result);
        double error2 = std::abs(a_result - rom_result);
        // Write results to file
        moe_out << f << "," << moe_result << "," << a_result << "," << k << "," << params[3] << "," << error1 / std::abs(a_result) << "\n";
        rom_out << f << "," << rom_result << "," << a_result << "," << k << "," << params[3] << "," << error2 / std::abs(a_result) << "\n";
        string win;
        if(error1 > error2)
            win1++;
        else
            win2++;
        terror1 += error1;
        terror2 += error2;
    }
    moe_out.close();
    rom_out.close();
    cout << "total abs_error over frequencies, " << terror1 << ", " << terror2 << "\n" << "with respective " << win1 << ", " << win2 << "\n";
    return 0;
}
