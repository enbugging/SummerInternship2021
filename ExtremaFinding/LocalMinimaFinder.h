#include <bits/stdc++.h>
using namespace std;

#include "../nomad/interfaces/CInterface/NomadStdCInterface.h"

struct Point
{
    vector<double> x;
    double value;
    bool clustered;
    bool operator<( const Point& other) const {
        return value < other.value;
    }
};

vector<Point> MLSL_MADS(
    int n, 
    double l, 
    function<double(vector<double>&)> objective_function);