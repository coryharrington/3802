#include <fstream>
#include <functional>
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;

double frequency[16] = {9.0, 9.6, 10.2, 10.8, 11.4, 12.0, 12.6, 13.2, 13.8, 14.40, 15.0, 15.6, 16.2, 16.8, 17.4, 18.0};
double resistance[16] = {6.58, 8.05, 9.9, 12.25, 15.34, 19.55, 25.54, 34.62, 49.54, 76.91, 135.3, 284.09, 618.48, 547.86, 252.51, 136.87};
double reactance[16] = {229.49, 190.98, 154, 117.74, 81.36, 43.84, 3.9, -40.33, -91.81, -155.33, -237.61, -330.25, -231.01, 229.8, 247.89, 152.36};

// Need to make dynamic1
int length(double * arr) { return 16; }

function<double(double)> lagrangian_interpolation(double * data_x, double * data_y)
{
    function<double(double)> interpolation = [data_x, data_y](double x)
    {
        double sum = 0;
        for(int i = 0; i < length(data_x); i++)
        {
            double num = 1;
            double denom = 1;
            for(int n = 0; n < length(data_x); n++)
            {
                if(n != i)
                {
                    num *= (x - data_x[n]);
                    denom *= (data_x[i] - data_x[n]);
                }
            }
            sum += (num * data_y[i]) / denom;
        }
        return sum;
    };
    return interpolation;
}

function<double(double)> quadratic_splines(double * data_x, double * data_y)
{
    double * a = data_y;
    double b[16], c[16];
    // Create coefficients, b0 is from secant line.
    b[0] = (data_y[1] - data_y[0]) / (data_x[1] - data_x[0]);
    for(int i = 1; i < length(data_y); i++)
    {
        b[i] = 2 * (a[i]-a[i-1]) / (data_x[i] - data_x[i-1]) - b[i-1];
        c[i] = 0.5 * (b[i] - b[i-1]) / data_x[i] - data_x[i-1];
    }
    function<double(double)> interpolation = [a, b, c, data_x](double x)
    {
        // Find interval
        int i = 0;
        while(x > data_x[i]) i++;
        // Return spline evaluated
        return a[i] + b[i] * (x - data_x[i]) + c[i] * pow(x - data_x[i], 2);
    };
    return interpolation;
}

double modified_regula_falsi(double tol, int iter, double guess, function<double(double)> f)
{
    double a, b, c_old, c_new, fa, fb, fc, result;
    c_old = guess;
    cout << "\nc is starting at : " << c_old;
    a = c_old - 0.1;
    b = c_old + 0.1;
    while(f(a)*f(b) > 0)
    {
        cout << "\nNo guarantee there is a zero in the given interval, increasing bracket by magnitude 0.2";
        a -= 0.1; b += 0.1;
    }
    for(int n = 1; n <= iter; n++)
    {
        fa = f(a);
        fb = f(b);
        c_new = (a * fb - b * fa) / (fb - fa);
        if(c_new == 0)
            cout << "\nConverging on the TRIVIAL solution!";
        fc = f(c_new);
        // Use IVT with midpoint to see where to walk
        cout << "\nc_n : " << c_new << ", c_n+1 : " << c_old << ", on " << a << " to " << b << "c_new : " << c_new << ", iter " << n;
        if(fa*fc > 0)
            a = c_new;
        else
            b = c_old;
        // Found zero within tolerance?
        if(fabs(fc) < tol)
            return c_new;
        c_old = c_new;
    }
    cout << "Reached maximum number of iterations, returning initial guess.";
    return guess;
}

int main(const int argc, const char** args)
{
    int n, maxiter;
    double tol;
    // Method returns an interpolation, or function<double(double)>
    function<double(double)> interpolation = NULL;
    function<function<double(double)>(double*, double*)> method;
    do
    {
        cout << \
            "\nPlease choose a method, { 0 : Lagrangian, 1 : Quadratic Splines, 2: Zero Find, 3 : Interpolate(x), 4 : Chebyshev nodes, _ : Help } : ";
        cin >> n;
        switch(n)
        {
            case 0: method = lagrangian_interpolation; break;
            case 1: method = quadratic_splines; break;
            case 2:
                if(interpolation != NULL)
                {
                    double guess;
                    cout << "\nPlease input a zero guess : ";
                    cin >> guess;
                    cout << "\nPlease input the tolerance on the solution : ";
                    cin >> tol;
                    cout << "\nPlease input the maximum number of iterations : ";
                    cin >> maxiter;
                    double result = modified_regula_falsi(tol, maxiter, guess, interpolation);
                    cout << "\nResult : " << result;
                }
                else
                    cout << "\nCannot zero find without an interpolation function.";
                break;
            case 3:
                if(interpolation != NULL)
                {
                    cout << "\nPlease enter a value to interpolate : ";
                    double temp;
                    cin >> temp;
                    cout << "\nf(" << temp << ") : " << interpolation(temp);
                }
                else
                    cout << "\nCannot interpolate without an interpolation function.";
                break;
            case 4:
                int option;
                double a, b;
                double * data;
                cout << "\nPlease enter the data_y in question, {0 : R, 1 : X} : ";
                cin >> option;
                // Data in question and get interval, try to find chebyshev nodes
                data = option ? reactance : resistance;
                a = data[0];
                b = data[length(data) - 1];
                for(int i = 0; i < length(data); i++)
                {
                    int k = 1;
                    double term;
                    do
                    {
                        term = 0.5 * (a + b) + 0.5 * (a - b) * cos((2*k - 1) * M_PI / (length(data) - 1));
                        // We will consider nodes within a half of a chebyshev node.
                        if(fabs(term - data[i]) < 0.5)
                            cout << "\nNode found at : " << data[i] << ", as cos(2*" << k << " - 1) * pi / 2n = " << term; break;
                        k++;
                    } while (k < length(data)); 
                }
                break;
            default:
                cout << "\nIf you do not have a function interpolation, you cannot run the zero finder. Modified Regula falsi.";
        }
        if(n == 0 || n == 1)
        {
            cout << "\nPlease enter the data_y in question, {0 : R, 1 : X} : ";
            int option;
            cin >> option;
            interpolation = option ? method(frequency, reactance) : method(frequency, resistance);
            // Testing
            for(double i = frequency[0]; i <= frequency[length(frequency) - 1]; i += 0.05) cout << i << "," << interpolation(i) << "\n";
            cout << '\n';
        }
    } while(!feof(stdin));
    return 0;
}
