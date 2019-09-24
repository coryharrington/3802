#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include "bessel.h"

using namespace std;

int main(int argc, const char* argv[])
{
    double range[argc];
    cout << "Please specify { 0 : JO(X), 1 : Y0(X), 2 : Zero(X) }\n";
    int offset, option, i;
    cin >> option;
    function<double(double, int)> bessel;
    switch(option)
    {
        case 0:
            bessel = first_bessel_0;
            offset = 0;
            break;
        case 1:
            bessel = second_bessel_0;
            offset = 302;
            break;
        case 2:
            bessel = zero_function;
            break;
        default: break;
    }
    vector<string> testing;
    for(i = 0; i < argc && i <= 302; i++)
    {
        double input = strtod(argv[i], NULL);
        char temp[32];
        sprintf(temp, "X: %8.1f, F(X): %8.6f", input, bessel(input, 10));
        testing.push_back(string(temp));
    }
    // Zero function testing
    if(option == 2)
    {
        for(i = 0; i < testing.size(); cout << '\n' << testing[i], i++);
        return 0;
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
