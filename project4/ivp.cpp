#include <iostream>
#include <fstream>
#include <functional>
#include <math.h>
#include <vector>

using namespace std;

void eulers(double h, double t0, double y0, function<double(double)> a, function<double(double, double)> f, vector<int> points)
{
    // map, n -> t0, y0, f(t0, y0), y_n+1
    // assume points are sorted
    double maximum = points.back();
    int maximum_n = ceil((maximum - t0) / h);
    tuple<double, double, double, double> table[maximum_n];
    // open file for writing
    ofstream e1_out, e2_out;
    e1_out.open("e1_out.csv");
    e2_out.open("e2_out.csv");
    // set tile of table
    e1_out << "tn" << ',' << "yn" << ',' << "f(tn yn)" << ',' << "y_n+1\n";
    e2_out << "tn" << ',' << "yn" << ',' << "f(tn yn)" << ',' << "y_n+1\n";
    double t_n, y_n, f_n;
    double sum = 0;
    int i = 0;
    for(int n = 0; n < maximum_n; n++)
    {
        // set table row = tuple of t_n, y_n, f_n, y_n+1
        t_n = t0 + n*h;
        if(n > 0)
            y_n = get<3>(table[n - 1]);
        else
            y_n = y0;
        f_n = f(t_n, y_n);
        table[n] = make_tuple(t_n, y_n, f_n, y_n + h * f_n);
        double ase = pow(a(get<0>(table[n])) - get<1>(table[n]), 2);
        sum += ase;
        // write to csv file, please move entries to save multiple runs of say different h's
        e2_out << get<0>(table[n]) << ',' << get<1>(table[n]) << ',' << get<2>(table[n]) << ',' << get<3>(table[n]) << '\n';
        if((int) t_n == points[i])
        {
            // write to csv file, please move entries to save multiple runs of say different h's, only at the points specified
            e1_out << get<0>(table[n]) << ',' << get<1>(table[n]) << ',' << get<2>(table[n]) << ',' << get<3>(table[n]) << '\n';
            // a is the analytical result, in other words the answer
            cout << get<0>(table[n]) << ", y_n : " << get<1>(table[n]) << " => RMS : " << sqrt(sum/n) << '\n';
            i++;
        }
    }
    cout << "maximum_n is " << maximum_n << '\n';
    cout << "Totall RMS over " << maximum_n << " points, RMS : " << sqrt(sum/maximum_n) << "\n";
    e1_out.close();
    e2_out.close();
}

void modified_eulers(double h, double t0, double y0, function<double(double)> a, function<double(double, double)> f, vector<int> points)
{
    // map, n -> t0, y0, f(t0, y0), y_n+1
    // assume points are sorted
    double maximum = points.back();
    int maximum_n = ceil((maximum - t0) / h);
    tuple<double, double, double, double> table[maximum_n];
    // open file for writing
    ofstream e1_out, e2_out;
    e1_out.open("me1_out.csv");
    e2_out.open("me2_out.csv");
    // set tile of table
    e1_out << "tn" << ',' << "yn" << ',' << "f(tn yn)" << ',' << "y_n+1\n";
    e2_out << "tn" << ',' << "yn" << ',' << "f(tn yn)" << ',' << "y_n+1\n";
    double t_n, y_n, f_n;
    double sum = 0;
    int i = 0;
    for(int n = 0; n < maximum_n; n++)
    {
        // set table row = tuple of t_n, y_n, f_n, y_n+1
        t_n = t0 + n*h;
        if(n > 0)
            y_n = get<3>(table[n - 1]);
        else
            y_n = y0;
        f_n = f(t_n, y_n);
        table[n] = make_tuple(t_n, y_n, f_n, y_n + (h/2) * (f_n + f(t_n + h, y_n + h*f_n)));
        double ase = pow(a(get<0>(table[n])) - get<1>(table[n]), 2);
        sum += ase;
        // write to csv file, please move entries to save multiple runs of say different h's
        e2_out << get<0>(table[n]) << ',' << get<1>(table[n]) << ',' << get<2>(table[n]) << ',' << get<3>(table[n]) << '\n';
        if((int) t_n == points[i])
        {
            // write to csv file, please move entries to save multiple runs of say different h's, only at the points specified
            e1_out << get<0>(table[n]) << ',' << get<1>(table[n]) << ',' << get<2>(table[n]) << ',' << get<3>(table[n]) << '\n';
            // a is the analytical result, in other words the answer
            cout << get<0>(table[n]) << ", y_n : " << get<1>(table[n]) << " => RMS : " << sqrt(sum/n) << '\n';
            i++;
        }
    }
    cout << "Totall RMS over " << maximum_n << " points, RMS : " << sqrt(sum/maximum_n) << "\n";
    e1_out.close();
    e2_out.close();
}

void rk4(double h, double t0, double y0, function<double(double)> a, function<double(double, double)> f, vector<int> points)
{
    // map, n -> t0, y0, f(t0, y0), y_n+1
    // assume points are sorted
    double maximum = points.back();
    int maximum_n = ceil((maximum - t0) / h);
    tuple<double, double, double, double> table[maximum_n];
    // open file for writing
    ofstream r1_out, r2_out;
    r1_out.open("r1_out.csv");
    r2_out.open("r2_out.csv");
    // set tile of table
    r1_out << "tn" << ',' << "yn" << ',' << "f(tn yn)" << ',' << "y_n+1\n";
    r2_out << "tn" << ',' << "yn" << ',' << "f(tn yn)" << ',' << "y_n+1\n";
    double t_n, y_n, f_n, k1, k2, k3, k4;
    double sum = 0;
    int i = 0;
    for(int n = 0; n < maximum_n; n++)
    {
        // set table row = tuple of t_n, y_n, f_n, y_n+1
        t_n = t0 + n*h;
        if(n > 0)
            y_n = get<3>(table[n - 1]);
        else
            y_n = y0;
        f_n = f(t_n, y_n);
        k1 = h * f_n;
        k2 = h * f(t_n + (h/2), y_n + (k1/2));
        k3 = h * f(t_n + (h/2), y_n + (k2/2));
        k4 = h * f(t_n + h, y_n + k3);
        table[n] = make_tuple(t_n, y_n, f_n, y_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4));
        double ase = pow(a(get<0>(table[n])) - get<1>(table[n]), 2);
        sum += ase;
        // write to csv file, please move entries to save multiple runs of say different h's
        r2_out << get<0>(table[n]) << ',' << get<1>(table[n]) << ',' << get<2>(table[n]) << ',' << get<3>(table[n]) << '\n';
        if((int) t_n == points[i])
        {
            // write to csv file, please move entries to save multiple runs of say different h's, only at the points specified
            r1_out << get<0>(table[n]) << ',' << get<1>(table[n]) << ',' << get<2>(table[n]) << ',' << get<3>(table[n]) << '\n';
            // a is the analytical result, in other words the answer
            cout << get<0>(table[n]) << ", y_n : " << get<1>(table[n]) << " => RMS : " << sqrt(sum/n) << '\n';
            i++;
        }
    }
    cout << "Totall RMS over " << maximum_n << " points, RMS : " << sqrt(sum/maximum_n) << "\n";
    r1_out.close();
    r2_out.close();
}

int main(int argc, char **argv)
{
    // Problem given is f(x, y) = y' = x(y^2)(e^-(x^4)), y(0) = 2 => analytically y = (-4 / sqrt(M_1_PI)) * (1 / (erf(sqrt(pow(x, 4))) - erf(0)))
    // Test a discrete subset of the family of solutions y'=y, and y''+
    if(argc != 1)
    {
        cout << "No paramaters are allowed\n";
        return 0;
    }
    // Analytic result f(x, y) = y'
    auto a = [](double x) { return (-4 / sqrt(M_1_PI)) * (1 / (erf(sqrt(pow(x, 4))) - erf(0))); };
    auto f = [](double x, double y) { return x * pow(y, 2) * exp(-1 * pow(x, 4)); };
    // Testing f(x, y) = y' = y
    auto b = [](double x) { return exp(x); };
    auto ft = [](double x, double y) { return y; };
    // Step size can be played with, but watch how quickly the RMS goes down with the step size in [0, 1] <-
    double h = 0.0001;
    // IV
    double t0 = 0;
    double y0 = 2;
    // Wide range, writen to csv and displayed to stdout
    vector<int> points;
    for(int i = 0; i < 11; i++) points.push_back(i);
    // Supply all information to our mehtod's, now with problem data, testing for y'= y commented out below
    eulers(h, 0, 1, b, ft, points);
    eulers(h, t0, y0, a, f, points);
    modified_eulers(h, 0, 1, b, ft, points);
    modified_eulers(h, t0, y0, a, f, points);
    rk4(h, 0, 1, b, ft, points);
    rk4(h, t0, y0, a, f, points);
    return 0;
}