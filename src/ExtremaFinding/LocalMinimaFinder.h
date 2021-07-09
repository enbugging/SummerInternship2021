#include <bits/stdc++.h>
using namespace std;

struct Point
{
    vector<double> x;
    double value;
    bool clustered;
    bool operator<( const Point& other) const {
        return value < other.value;
    }
};

vector<Point> MLSL(
    int n, 
    double l, 
    function<double(vector<double>&)> objective_function);