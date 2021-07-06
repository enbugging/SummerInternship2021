#include <bits/stdc++.h>
using namespace std;

class LocalMinimaFinder
{
    public:
        struct Point;

    private:
        struct MADS;

        unsigned int seed;
        mt19937 rng;
        double l;
        int n;
        bool calculated;
        vector<Point> X_star;

        double random_step(
            double min,
            double max);
        
        bool min_finder_starting_point_check(
            vector<double>& x,
            vector<vector<double> >& X_star,
            vector<vector<double> >& V, 
            double min_norm_squared_between_minima, 
            double r
        );

        function <double(vector<double>)> objective_function;

    public:
        LocalMinimaFinder(int dimension, double length, std::function <double(vector<double>)> target_function)
        {
            // random seed, can be fixed for rerun of experiments
            unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
            mt19937 rng(seed);
            l = length;
            calculated = false;
            n = dimension;
            objective_function = target_function;
        }

        vector<Point> MLSL_MADS();
};