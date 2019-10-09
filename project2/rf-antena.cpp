#include <functional>
#include <iostream>
// #include <armadillo>

using namespace std;
// using namespace arma;

function<double(double)> lagrangian_interpolation(double tol)
{
    return NULL;
}

function<double(double)> cubic_splines(double tol)
{
    return NULL;
}

int main()
{
    int n;
    double tol;
    // Method returns an interpolation, or function<double(double)>
    function<double(double)> results;
    function<function<double(double)>(double)> method;
    do
    {
        do {
            cout << "\nPlease choose a method, { 0 : Lagrangian, 1 : Cubic Splines, _ : Help } : ";
            cin >> n;
            switch(n)
            {
                case 0: method = lagrangian_interpolation; break;
                case 1: method = cubic_splines; break;
                default:
                    cout << "\n If you do not have a function interpolation, you cannot run the zero finder.";
            }
        } while(n == 0 || n == 1);
        cout << "\nPlease enter the tolerance on the solution : ";
        cin >> tol;
        method(tol);
        cout << '\n';
    } while (!feof(stdin));
    return 0;
}
