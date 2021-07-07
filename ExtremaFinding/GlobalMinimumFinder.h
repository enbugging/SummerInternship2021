#include <bits/stdc++.h>
using namespace std;

class GlobalMinimumFinder
{
    private:
        vector<int> multiplicities;
        mt19937 rng;

        double random_step(
            double min,
            double max);
        
        double coefficient(
            double t, 
            double threshold);

        double rmse_with_threshold(
            double threshold, 
            vector<double>& force_constants,
            vector<double>& angles,
            vector<double>& quantum_mechanics_data_points);
    public:
        GlobalMinimumFinder()
        {
            multiplicities = {1, 2, 3, 4, 6};
            // random seed, can be fixed for rerun of experiments
            unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
            mt19937 rng1(seed);
            rng = rng1;
        }

        double force_field_calculate(
            vector<double>& force_constants,
            double angle
        );

        double rmse(
            vector<double>& force_constant,
            vector<double>& angles,
            vector<double>& quantum_mechanics_data_points);

        vector<double> simulated_annealing(
            vector<double>& angles,
            vector<double>& quantum_mechanics_data_points, 
            int number_of_terms,
            int number_of_steps, 
            double initial_temperature, 
            double initial_radius, 
            double threshold = 0
        );

        vector<double> threshold_accepting(
            vector<double>& angles,
            vector<double>& quantum_mechanics_data_points, 
            int number_of_terms,
            int number_of_steps, 
            double initial_temperature, 
            double initial_radius, 
            double threshold = 0
            );

        vector<double> simulated_annealing_with_threshold(
            vector<double>& angles,
            vector<double>& quantum_mechanics_data_points, 
            int number_of_terms,
            int number_of_steps, 
            double initial_temperature, 
            double initial_radius, 
            double threshold = 0);
};