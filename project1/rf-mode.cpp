#include <iostream>
#include <vector>
#include "bessel.h"

using namespace std;

const double a = 0.075;
const double b = 0.325;
static double tol;

int main(int argc, const char* argv[])
{
    int n;
    double k;
    vector<double> results;
    cout << "Please enter the maximum number of iterations : ";
    cin >> n;
    cout << "\nPlease enter the tolerance on the solution : ";
    cin >> tol;
    cout << "\nPlease enter the an initial estimate of the solution : ";
    cin >> k;
    return 0;
}
