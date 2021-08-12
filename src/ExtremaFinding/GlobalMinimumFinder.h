#include <bits/stdc++.h>
using namespace std;

static vector<double> dummy;

vector<double> simulated_annealing(
    function<double(vector<double>&)> objective_function, 
    int D, 
    double r, 
    int number_of_steps = 4000);