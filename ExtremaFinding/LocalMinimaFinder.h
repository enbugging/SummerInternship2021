#include <bits/stdc++.h>
using namespace std;


#include "../nomad/interfaces/CInterface/NomadStdCInterface.h"

class LocalMinimaFinder
{
    public:
        function <double(vector<double>&)> objective_function;

    private:
        struct Point;
        
        unsigned int seed;
        mt19937 rng;
        double l;
        int n;
        bool calculated;
        vector<Point> X_star;

        double random_step(
            double min,
            double max);
        
        vector<double> MADS(
            vector<double>& x0
        );

        double norm(
            Point a, 
            Point b);
        
        void clustering(
            vector<Point>& C_twiddle, 
            vector<Point>& clustered,
            double r);

    public:
        LocalMinimaFinder(int dimension, double length, std::function <double(vector<double>&)> target_function)
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

/*
#include <bits/stdc++.h>
using namespace std;


#include "../nomad/interfaces/CInterface/NomadStdCInterface.h"

class LocalMinimaFinder
{
    public:
        struct Point;

    private:
        double l;
        int n;
        bool calculated;
        function <double(array<double, n>)> objective_function;
        unsigned int seed;
        mt19937 rng;
        vector<Point> X_star;

        double random_step(
            double min,
            double max);


        bool LocalMinimaFinder::objective_function_for_Nomad(
            int nb_inputs, 
            double *x, 
            int nb_outputs, 
            double *bb_outputs, 
            bool *count_eval, 
            NomadUserDataPtr data);
        
        vector<double> LocalMinimaFinder::MADS(
    array<int, 6> x0
);

    public:
        LocalMinimaFinder(int dimension, double length, std::function <double(array<double, n>)> target_function)
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
*/