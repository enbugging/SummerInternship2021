#include <bits/stdc++.h>
using namespace std;

#include "ExtremaFinding.h"

vector<int>
	multiplicities = {1, 2, 3, 4, 6};
// random seed, can be fixed for rerun of experiments
unsigned int
    seed = chrono::steady_clock::now().time_since_epoch().count();
mt19937 rng(seed);

int number_of_steps = 50000;

double force_field_calculate(
    vector<double> force_constants,
    double angle)
{
    double E = 0;
    for (int i = 0; i < (int)force_constants.size(); i++)
    {
        E += force_constants[i] * (1 + cos(multiplicities[i] * angle));
    }
    return E;
}

double square_error(
    vector<double> force_constants,
    vector<double> angles,
    vector<double> quantum_mechanics_data_points)
{
    double sum_of_squares_of_error = 0;
    for (int i = 0; i < (int)angles.size(); i++)
    {
        double error = force_field_calculate(force_constants, angles[i]) - quantum_mechanics_data_points[i];
        sum_of_squares_of_error += error * error;
    }
    return sum_of_squares_of_error;
}

double random_step(
    double min,
    double max)
{
    double f = (double)rng() / ULONG_MAX;
    return min + f * (max - min);
}

vector<double> simulatedAnnealing(
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double radius, 
    int number_of_terms,
    vector<double> angles,
    vector<double> quantum_mechanics_data_points)
{
    // initialization
    double best_square_error = numeric_limits<double>::max();
    vector<double> force_constants;
    force_constants.resize(number_of_terms);

    // random initial set of parameters
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = random_step(-initial_radius, initial_radius);
    }

    for (int i = 0; i < number_of_steps; i++)
    {
        double T = initial_temperature * exp(-4.0 * i / number_of_steps);
        vector<double> new_force_constants = force_constants;
        
        // random distortion to the parameters
        for (int i = 0; i < number_of_terms; i++)
        {
            new_force_constants[i] += random_step(-radius, radius);
            // random_step(-radius / (1 + i), radius / (1 + i));
            // random_step(-radius, radius);
        }
        double new_square_error = square_error(new_force_constants, angles, quantum_mechanics_data_points);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        // printf("%lf %lf %lf\n", T, new_square_error, best_square_error);
        if (new_square_error <= best_square_error || 
            random_step(0, 1) <= exp(-(new_square_error - best_square_error) / T))
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }

        vector<double> molecular_mechanics_data_points;
        molecular_mechanics_data_points.resize((int) angles.size());
            
        if(0)
        {
            for (int j = 0; j < (int) angles.size(); j++)
            {
                molecular_mechanics_data_points[j] = force_field_calculate(new_force_constants, angles[j]);
            }
            for (int j = 0; j < (int) angles.size(); j++)
            {
                printf("%lf %lf\n", molecular_mechanics_data_points[j], quantum_mechanics_data_points[j]);
            }
            printf("SUMMARY\n");
            printf("---------------------------------------------\n");
            printf("Temperature: %lf\n", T);
            double sum_of_squares_of_error = square_error(new_force_constants, angles, quantum_mechanics_data_points);
            printf("\nSquare error: %lf\n", sqrt(sum_of_squares_of_error/(int) angles.size()));
            sum_of_squares_of_error = square_error(force_constants, angles, quantum_mechanics_data_points);
            printf("\nBest square error: %lf\n", sqrt(sum_of_squares_of_error/(int) angles.size()));
        }
    }

    return force_constants;
}