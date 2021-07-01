#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "ExtremaFinding.h"

vector<int>
	multiplicities = {1, 2, 3, 4, 6};
// random seed, can be fixed for rerun of experiments
unsigned int
    seed = chrono::steady_clock::now().time_since_epoch().count();
mt19937 rng(seed);

normal_distribution<double> gaussian(0.0, 1.0);

double force_field_calculate(
    vector<double> force_constants,
    double angle)
{
    double E = 0;
    for (int i = 0; i < (int)force_constants.size(); i++)
    {
        E += force_constants[i] * (1 + cos(multiplicities[i] * angle / 180.0 * M_PI));
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

double coefficient(
    double t, 
    double threshold
)
{
    double border_width = 0.002;
    if (abs(t) >= threshold + border_width) return 1;
    if (abs(t) < threshold) return 0;
    return (1.0 - cos(M_PI / border_width * (abs(t) - threshold)))/2.0;
}

double square_error_with_threshold(
    double threshold, 
    vector<double> force_constants,
    vector<double> angles,
    vector<double> quantum_mechanics_data_points)
{
    vector<double> new_force_constants = force_constants;
    for (int i = 0; i < (int) force_constants.size(); i++)
    {
        new_force_constants[i] *= coefficient(force_constants[i], threshold);
    }
    double sum_of_squares_of_error = 0;
    for (int i = 0; i < (int)angles.size(); i++)
    {
        double error = force_field_calculate(new_force_constants, angles[i]) - quantum_mechanics_data_points[i];
        sum_of_squares_of_error += error * error;
    }
    return sum_of_squares_of_error;
}

double random_step(
    double min,
    double max)
{
    double f = (double)rng() / UINT_MAX;
    return min + f * (max - min);
}


vector<double> simulated_annealing(
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
        double T = initial_temperature * pow(((double) (number_of_steps - i) / number_of_steps), 2);
        vector<double> new_force_constants = force_constants;
        
        // random distortion to the parameters
        for (int j = 0; j < number_of_terms; j++)
        {
            new_force_constants[j] += 
            T * tan(random_step(-M_PI/2, M_PI/2));
        }
        double new_square_error = square_error(new_force_constants, angles, quantum_mechanics_data_points);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        if (new_square_error <= best_square_error || 
            random_step(0, 1) <= exp(-(new_square_error - best_square_error) / T))
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }
    }

    return force_constants;
}

vector<double> threshold_accepting(
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
        double T = initial_temperature * pow(((double) (number_of_steps - i) / number_of_steps), 2);
        vector<double> new_force_constants = force_constants;
        
        // random distortion to the parameters
        for (int j = 0; j < number_of_terms; j++)
        {
            new_force_constants[j] += 
            T * tan(random_step(-M_PI/2, M_PI/2));
        }
        double new_square_error = square_error(new_force_constants, angles, quantum_mechanics_data_points);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        if (new_square_error <= best_square_error + T)
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }
    }

    return force_constants;
}

vector<double> simulated_annealing_with_threshold(
    double threshold, 
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double radius, 
    int number_of_terms,
    vector<double> angles,
    vector<double> quantum_mechanics_data_points)
{
    vector<double> best_force_constants;
    best_force_constants.resize(number_of_terms);
    double optimal_square_error = numeric_limits<double>::max();

    for (int mask = 1; mask < (1<<number_of_terms); mask++)
    {
        // initialization
        vector<double> force_constants;
        force_constants.resize(number_of_terms);
        double best_square_error = numeric_limits<double>::max();

        // random initial set of parameters
        for (int i = 0; i < number_of_terms; i++)
        {
            if((mask >> i) & 1)
            {
                force_constants[i] = random_step(-initial_radius, initial_radius);
            }
            else 
            {
                force_constants[i] = 0.0;
            }
        }

        for (int i = 0; i < number_of_steps; i++)
        {
            double T = initial_temperature * pow(((double) (number_of_steps - i) / number_of_steps), 2);
            vector<double> new_force_constants = force_constants;
            
            // random distortion to the parameters
            for (int j = 0; j < number_of_terms; j++)
            {
                if((mask >> j) & 1) new_force_constants[j] += T * tan(random_step(-M_PI/2, M_PI/2));
            }
            double new_square_error = square_error_with_threshold(threshold, new_force_constants, angles, quantum_mechanics_data_points);

            // if the new set of parameters is better, then we accept
            // else, we accept, with a probability corresponding to the temperature
            if (new_square_error <= best_square_error || 
                random_step(0, 1) <= exp(-(new_square_error - best_square_error) / T))
            {
                best_square_error = new_square_error,
                force_constants = new_force_constants;
            }
        }

        if (optimal_square_error > best_square_error)
        // if (mask == (1<<number_of_terms) - 1)
        {
            optimal_square_error = best_square_error;
            best_force_constants = force_constants;
        }

        cerr << "Result for choosing harmonic(s) ";
        for(int i = 0; i < number_of_terms; i++)
        {
            if((mask >> i) & 1) cerr << i << " ";
            else cerr << "- ";
        }
        cerr << ": ";
        for (int i = 0; i < number_of_terms; i++)
        {
            cerr << force_constants[i] << ' ';
        }
        cerr << "\nRMSE: " << sqrt(square_error(force_constants, angles, quantum_mechanics_data_points)/(double) angles.size()) << '\n';
        force_constants.clear();
    }
    
    // setting trapped coefficients to be zero
    for (int i = 0; i < number_of_terms; i++)
    {
        if (abs(best_force_constants[i]) < threshold)
        {
            best_force_constants[i] = 0;
        }
    }

    return best_force_constants;
}