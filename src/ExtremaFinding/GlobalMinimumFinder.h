#include <bits/stdc++.h>
using namespace std;

static vector<double> dummy;

/* 
Simulated annealing optimization.
objective_function: function to be optimized
D: the dimension of the input
r: the initial radius
number_of_steps: number of steps, default to be
1e6.
 */
vector<double> simulated_annealing(
    function<double(vector<double>&)> objective_function, 
    int D, 
    double r, 
    int number_of_steps = 1000000);