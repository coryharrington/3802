#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

double tol;

int main()
{
    int n;
    function<void(double)> method;
    do
    {
        do {
            cout << "\nPlease choose a method, { 0 : NR iter, 1 : MCB, 2: RF, 3 : SS, _ : Help } : ";
            cin >> n;
            switch(n)
            {
                case 0: 
                    method = newton_raphson;
                    cout << "\nNR iter skip's most roots by itself, so I implimented a decreasing learning rate if it overshoots.";
                    cout << "\nYou must choose an inital value less than it to gaurantee the learning algorithm.";
                    break;
                case 1: method = monte_carlo_bisection; break;
                case 2: method = regula_falsi; break;
                case 3: 
                    method = stephansons_secant;
                    cout << "\nSS iter converges like NR iter, and so I kept the learning rate approach.";
                    break;
                default:
                    cout << "\n0 : Newton Raphson Iteration\n1 : Monte Carlo Bisection\n2 : Regula Falsi (Illinois)\n3 : Stephansons Secant";
            }
        } while(n == 0 || n == 1);
        cout << "\nPlease enter the tolerance on the solution : ";
        cin >> tol;
        cout << "\nPlease enter the an initial estimate of the solution : ";
        cin >> k;
        cout << "\nFinding wavenumber solutions (rad / in)";
        // Initial conditions, run proceedure (void), inform results to user
        results = vector<double>();
        knew = kold = k;
        method(k);
        for(int i = 0; i < results.size(); cout << "\nSol #" << i << " : " << results[i], i++);
        if(results.empty())
            cout << "\nNo solutions found... Divergent";
        cout << '\n';
    } while (!feof(stdin));
    return 0;
    // cout << A*B.t() << endl;
    return 0;
}