#include <iostream>
#include <string>
#include <vector>
#include "bessel.h"

using namespace std;

int main(int argc, const char* argv[])
{
    double range[argc];
    cout << "Please specify { 0 : JO(X), 1 : Y0(X) }\n";
    int offset, option, i;
    cin >> option;
    function<double(double, int)> bessel;
    if(option == 0) {
        bessel = first_bessel_0;
        offset = 0;
    } else {
        bessel = second_bessel_0;
        offset = 302;
    }
    vector<string> testing;
    for(i = 0; i < argc && i <= 302; i++)
    {
        double input = strtod(argv[i], NULL);
        char temp[32];
        sprintf(temp, "X: %8.1f, F(X): %8.6f", input, bessel(input, 10));
        testing.push_back(string(temp));
    }
    // Start at the correct column from the AWK'ed input
    for(i += offset; i < argc && i < 604 + offset; i++) 
    {
        char temp[32];
        sprintf(temp, " = %8.6f\n", strtod(argv[i], NULL));
        cout << testing[i - offset - 302] << string(temp);
    }
    return 0;
}
