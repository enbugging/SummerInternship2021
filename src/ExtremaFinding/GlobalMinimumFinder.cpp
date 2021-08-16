#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "GlobalMinimumFinder.h"
#ifndef _NR_
#define _NR_
#include "../numerical_recipes/nr.h" // numerical recipes
#endif

vector<double> standard_cauchy_vector(
    int D, 
    mt19937& rng)
{
    cauchy_distribution<double> cauchy_dist (0.0,1.0);
    vector<double> v(D-1), w(D-1), x(D);
    for (int i = 0; i < D-1; i++)
    {
        v[i] = cauchy_dist(rng);
    }
    for (int i = 0; i < D-1; i++)
    {
        student_t_distribution<double> student_t_dist (i + 2);
        w[i] = student_t_dist(rng);
    }
    x[0] = v[0];
    for (int i = 0; i < D-1; i++)
    {
        v[i] *= v[i];
        if (i)
        {
            v[i] += v[i-1];
        }
    }
    for (int i = 1; i < D; i++)
    {
        x[i] = w[i-1] * sqrt((1 + v[i-1]) / (1 + i));
    }
    return x;
}

/* 
Implementation follows Fast simulated annealing with a multivariate 
Cauchy distribution and the configuration’s initial temperature
May 2015Journal- Korean Physical Society 66(10):1457-1466
DOI:10.3938/jkps.66.1457
Retrieved from https://www.researchgate.net/publication/279241787
 */
function<double(vector<double>&)> target_function;

DP objective_function_wrapper(Vec_I_DP &x) {
    vector<double> y;
    y.resize(x.size());
    for (int i = 0; i < (int) x.size(); i++)
    {
        y[i] = x[i];
    }
    return target_function(y);
}

vector<double> simulated_annealing(
    function<double(vector<double>&)> objective_function, 
    int D, 
    double r, 
    int number_of_steps)
{
    // Power's algorithm initialization
    DP FTOL = 1.0e-12;
    const int IMAXSTEP = 1000000;
    int iter;
    DP fret;
    Vec_DP p(0.0, D); // variables
    Mat_DP xi(D, D);
    target_function = objective_function;

    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(seed);
    uniform_real_distribution<double> uniform_init (-r, r), uniform (0, 1);

    vector<double> x(D), dx(D), x_dx(D);
    for (int i = 0; i < D; i++)
    {
        x[i] = uniform_init(rng);
    }

    double 
        beta = 100, 
        T_0g = r / (2 * tan(M_PI/(2 * (1 + beta)))), 
        former_value = objective_function(x);
    for (int i = 2; i <= number_of_steps; i++)
    {
        dx = standard_cauchy_vector(D, rng);
        for (int j = 0; j < D; j++)
        {
            x_dx[j] = x[j] + T_0g / i * dx[j];
        }
        double 
            new_value = objective_function(x_dx), 
            delta_f = new_value - former_value;
        if (uniform(rng) <= exp(-delta_f * i / T_0g))
        {
            former_value = new_value;
            x = x_dx;
        }
    }

    // further polishing using Powell's algorithm, provided by Professor Aleksandrov
    for (int i = 0; i < D; i++)
    {
        p[i] = x[i];
    }
    for (int i = 0; i < D; i++)
    {
        for (int j = 0; j < D; j++)
        {
            xi[i][j] = (i == j ? x[i] : 0.0);
        }
    }
    NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, objective_function_wrapper);
    for (int i = 0; i < D; i++)
    {
        x[i] = p[i];
    }

    return x;
}