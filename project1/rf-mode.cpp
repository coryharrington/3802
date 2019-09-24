#include <functional>
#include <iostream>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <vector>
#include "bessel.h"

using namespace std;

static double k, knew, kold, maxiter, tol;
vector<double> results;

double random_in(double a, double b)
{
    uniform_real_distribution<double> range(a, b);
    default_random_engine ng;
    return range(ng);
}

// Cannot give initial value > first root, or it will skip such root, finding non-consecutive roots for the rest.
// Quadratic convergence
void newton_raphson(double k0)
{
    double alpha = 1;
    for(int n = 1; n <= maxiter; n++)
    {
compute_0:
        double fold = zero_function(kold, 10);
        double div = zero_function_prime(kold, 10);
        // Move zero_function_prime off zero
        if(div == 0)
            knew += 0.1;
        else
            knew = kold - alpha * (fold / zero_function_prime(kold, 10));
        double fnew = zero_function(knew, 10);
        cout << "\nknew : " << knew << ", kold : " << kold << ", f(knew) : " << fnew << ", f(kold) : " << fold << ", " << n;
        // Overshoot ? recompute with alpha_n+1 << alpha_n
        //&& alpha >= 0.1
        if(fnew * fold < 0 && alpha >= 0.1)
        {
            alpha *= 0.1;
            cout << "\nAdjusting learning rate because : " << fnew * fold;
            goto compute_0;
        }
        kold = knew;
        if(fabs(fnew) < tol)
        {
            results.push_back(knew);
            cout << "\nSol #" << results.size() << " found!";
            goto again_0;
        }
    }
    cout << "\nMaximum iterations reached";
    return;
again_0:
    if(results.size() < 4)
    {
        // Try the next integer multiple, as we are finding eigen's in a sense
        kold += (results[0] - 1);
        newton_raphson(kold);
    } else return;
}

// Gaurantee convergence from IVT. Expectation running time <= O(log(n))
void monte_carlo_bisection(double k0)
{
    double a = k0 - 1;
    double b = k0 + 1;
    double fa, fb, fm;
bracket_0:
    if(a <= 0)
    {
        cout << "\nZero function undefined at the domain 0";
        cout << "\nPlease enter the a different wave number : ";
        cin >> k;
        results = vector<double>();
        monte_carlo_bisection(k);
    }
    if(zero_function(a, 10)*zero_function(b, 10) > 0)
    {
        cout << "\nNo guarantee there is a zero in the given interval, increasing bracket";
        a--;
        b++;
        goto bracket_0;
    }
    for(int n = 1; n <= maxiter; n++)
    {
        cout << "\na is " << a << ", b is " << b;
        knew = random_in(a, b);
        if(knew == 0)
            cout << "\nConverging on the TRIVIAL solution! No such eigen";
        fa = zero_function(a, 10);
        fb = zero_function(b, 10);
        fm = zero_function(knew, 10);

        // Use IVT with random "midpoint" to see where to walk
        cout << "\nknew : " << knew << ", kold : " << kold << ", on " << a << " to " << b << "f(knew) : " << fm << ", iter " << n;
        if(fa*fm > 0)
            a = knew;
        else
            b = knew;
        if(fabs(fm) < tol)
        {
            results.push_back(knew);
            cout << "\nSol #" << results.size() << " found!";
            break;
        }
        // Get out of converging to the trivial solution or getting stuck
        if(fabs(knew - kold) < tol)
            b++;
        kold = knew;
    }
    if(results.empty())
        cout << "\nMaximum iterations reached";
    else
        if(results.size() < 4)
        {
            // Try the next integer multiple, as we are finding eigen's in a sense
            kold += results[0];
            monte_carlo_bisection(kold);
        }
}

// Gaurantee of convergence, at a superlinear rate, THIS THIS BLOWS EVERYTHING ELSE OUT THE WATER
// (Just kidding, SS method is great)
void regula_falsi(double k0)
{
    double a = k0 - 1;
    double b = k0 + 1;
    double fa, fb, fc;
bracket_1:
    if(a <= 0)
    {
        cout << "\nZero function undefined at the domain 0";
        cout << "\nPlease enter the a different wave number : ";
        cin >> k;
        results.pop_back();
        regula_falsi(k);
    }
    if(zero_function(a, 10)*zero_function(b, 10) > 0)
    {
        cout << "\nNo guarantee there is a zero in the given interval, increasing bracket";
        a--;
        b++;
        goto bracket_1;
    }
    for(int n = 1; n <= maxiter; n++)
    {
        cout << "\na is " << a << ", b is " << b;
        fa = zero_function(a, 10);
        fb = zero_function(b, 10);
        // Move dx from vanishing
        if(fb - fa == 0)
            a += 0.1;
        knew = (a * fb - b * fa) / (fb - fa);
        if(knew == 0)
            cout << "\nConverging on the TRIVIAL solution! No such eigen";
        fc = zero_function(knew, 10);

        // Use IVT with random "midpoint" to see where to walk
        cout << "\nknew : " << knew << ", kold : " << kold << ", on " << a << " to " << b << "f(knew) : " << fc << ", iter " << n;
        if(fa*fc > 0)
            a = knew;
        else
            b = knew;
        if(fabs(fc) < tol)
        {
            results.push_back(knew);
            cout << "\nSol #" << results.size() << " found!";
            break;
        }
        // Get out of converging to the trivial solution or getting stuck
        if(fabs(knew -kold) < tol)
            b++;
        kold = knew;
    }
    if(results.empty())
        cout << "\nMaximum iterations reached";
    else
        if(results.size() < 4)
        {
            // Try the next integer multiple, as we are finding eigen's in a sense.
            knew = kold = knew + results[0];
            regula_falsi(knew);
        }
}

// Quadratic convergence, like Newton Raphson
void stephansons_secant(double k0)
{
        double alpha = 1;
    for(int n = 1; n <= maxiter; n++)
    {
compute_1:
        double fold = zero_function(kold, 10);
        double div = zero_function((kold + zero_function(kold, 10)), 10) - fold;
        // Move zero_function_prime off zero
        if(div == 0)
            knew += 0.1;
        else
            knew = kold - (pow(fold, 2) / div);
        double fnew = zero_function(knew, 10);
        cout << "\nknew : " << knew << ", kold : " << kold << ", f(knew) : " << fnew << ", f(kold) : " << fold << ", " << n;
        // Overshoot ? recompute with alpha_n+1 << alpha_n
        // && alpha >= 0.1
        if(fnew * fold < 0 && alpha >= 0.1)
        {
            alpha *= 0.1;
            cout << "\nAdjusting learning rate because : " << fnew * fold;
            goto compute_1;
        }
        kold = knew;
        if(fabs(fnew) < tol)
        {
            results.push_back(knew);
            cout << "\nSol #" << results.size() << " found!";
            goto again_1;
        }
    }
    cout << "\nMaximum iterations reached";
    return;
again_1:
    if(results.size() < 4)
    {
        // Try the next integer multiple, as we are finding eigen's in a sense
        kold += (results[0] - 1);
        stephansons_secant(kold);
    } else return;
}

// ENTRYPOINT (Main function)
int main(int argc, const char* argv[])
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
        } while(n < 0 || n > 3);
        cout << "\nPlease enter the maximum number of iterations : ";
        cin >> maxiter;
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
}
