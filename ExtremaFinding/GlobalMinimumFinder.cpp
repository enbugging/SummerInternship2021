#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "GlobalMinimumFinder.h"

#define random_step(rng,min,max) min + (double)rng() / UINT_MAX * (max - min)

/*	-----------------------------------------------------------------------------------
 *	Utility functions, supporthing simulated-annealing-related functions
 * */
vector<int> multiplicities = {1, 2, 3, 4, 6};

double force_field_calculate(
    vector<double>& force_constants,
    double angle)
{
    double E = 0;
    for (int i = 0; i < (int)force_constants.size(); i++)
    {
        E += force_constants[i] * (1 + cos(multiplicities[i] * angle / 180.0 * M_PI));
    }
    return E;
}

double rmse(
    vector<double>& force_constants,
    vector<double>& angles,
    vector<double>& quantum_mechanics_data_points)
{
    if (angles.empty()) return 0;
    double sum_of_squares_of_error = 0;
    for (int i = 0; i < (int)angles.size(); i++)
    {
        double error = force_field_calculate(force_constants, angles[i]) - quantum_mechanics_data_points[i];
        sum_of_squares_of_error += error * error;
    }
    sum_of_squares_of_error = sqrt(sum_of_squares_of_error/angles.size());
    return sum_of_squares_of_error;
}

double coefficient(
    double t, 
    double threshold)
{
    double
        border_width = 0.001, 
        A, 
        B, 
        extra_wall = 0, 
        center = border_width + threshold;
    if (abs(t) >= center) A = 1;
    else if (abs(t) < threshold) A = 0;
    else A = (1.0 - cos(M_PI / border_width * (abs(t) - threshold)))/2.0;

    if (abs(abs(t) - center) < 2 * border_width) B = extra_wall * (1 + cos(M_PI / (2 * border_width) * (abs(t) - center)))/2.0;
    else B = 0;
    return A + B;
}

double rmse_with_threshold(
    double threshold, 
    vector<double>& force_constants,
    vector<double>& angles,
    vector<double>& quantum_mechanics_data_points)
{
    if (angles.empty()) return 0;
    vector<double> new_force_constants = force_constants;
    double c = 0.02;
    for (int i = 0; i < (int) force_constants.size(); i++)
    {
        new_force_constants[i] *= coefficient(force_constants[i], threshold);
    }
    double sum_of_squares_of_error = 0;
    for (int i = 0; i < (int) angles.size(); i++)
    {
        double error = force_field_calculate(new_force_constants, angles[i]) - quantum_mechanics_data_points[i];
        sum_of_squares_of_error += error * error;
    }
    sum_of_squares_of_error = sum_of_squares_of_error/angles.size();
    for (int i = 0; i < (int) force_constants.size(); i++)
    {
        sum_of_squares_of_error += c * coefficient(force_constants[i], threshold) *  coefficient(force_constants[i], threshold);
    }
    return sum_of_squares_of_error;
}

/*	-----------------------------------------------------------------------------------
 *	Simulated annealing and its variants. The implementations used Cauchy mutation scheme
 *  for best stability and power.
 * */
vector<double> simulated_annealing(
    vector<double>& angles,
    vector<double>& quantum_mechanics_data_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double threshold)
{
    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // initialization
    vector<double> force_constants;
    force_constants.resize(number_of_terms);
    double best_square_error = numeric_limits<double>::max();

    // random initial set of parameters
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = random_step(rng, -initial_radius, initial_radius);
    }

    for (int i = 0; i < number_of_steps; i++)
    {
        double T = initial_temperature * pow(((double) (number_of_steps - i) / number_of_steps), 2);
        vector<double> new_force_constants = force_constants;
        
        // random distortion to the parameters
        for (int j = 0; j < number_of_terms; j++)
        {
            new_force_constants[j] += T * tan(random_step(rng, -M_PI/2, M_PI/2));
        }
        double new_square_error = rmse_with_threshold(threshold, new_force_constants, angles, quantum_mechanics_data_points);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        if (new_square_error <= best_square_error || 
            random_step(rng, 0, 1) <= exp(-(new_square_error - best_square_error) / T))
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }
    }

    // setting trapped coefficients to be zero
    for (int i = 0; i < number_of_terms; i++)
    {
        if (abs(force_constants[i]) < threshold)
        {
            force_constants[i] = 0;
        }
    }

    return force_constants;
}

vector<double> threshold_accepting(
    vector<double>& angles,
    vector<double>& quantum_mechanics_data_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double threshold)
{
    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // initialization
    vector<double> force_constants;
    force_constants.resize(number_of_terms);
    double best_square_error = numeric_limits<double>::max();

    // random initial set of parameters
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = random_step(rng, -initial_radius, initial_radius);
    }

    for (int i = 0; i < number_of_steps; i++)
    {
        double T = initial_temperature * pow(((double) (number_of_steps - i) / number_of_steps), 2);
        vector<double> new_force_constants = force_constants;
        
        // random distortion to the parameters
        for (int j = 0; j < number_of_terms; j++)
        {
            new_force_constants[j] += T * tan(random_step(rng, -M_PI/2, M_PI/2));
        }
        double new_square_error = rmse_with_threshold(threshold, new_force_constants, angles, quantum_mechanics_data_points);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        if (new_square_error <= best_square_error + T)
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }
    }

    // setting trapped coefficients to be zero
    for (int i = 0; i < number_of_terms; i++)
    {
        if (abs(force_constants[i]) < threshold)
        {
            force_constants[i] = 0;
        }
    }

    return force_constants;
}

vector<double> simulated_annealing_with_threshold(
    vector<double>& angles,
    vector<double>& quantum_mechanics_data_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double threshold)
{
    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(seed);
    
    vector<double> best_force_constants;
    best_force_constants.resize(number_of_terms);
    // double optimal_square_error = numeric_limits<double>::max();

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
                force_constants[i] = random_step(rng, -initial_radius, initial_radius);
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
                if((mask >> j) & 1) new_force_constants[j] += T * tan(random_step(rng, -M_PI/2, M_PI/2));
            }
            double new_square_error = rmse_with_threshold(threshold, new_force_constants, angles, quantum_mechanics_data_points);

            // if the new set of parameters is better, then we accept
            // else, we accept, with a probability corresponding to the temperature
            if (new_square_error <= best_square_error || 
                random_step(rng, 0, 1) <= exp(-(new_square_error - best_square_error) / T))
            {
                best_square_error = new_square_error,
                force_constants = new_force_constants;
            }
        }

        //if (optimal_square_error > best_square_error)
        if (mask == (1<<number_of_terms) - 1)
        {
            // optimal_square_error = best_square_error;
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
        cerr << "\nRMSE: " << best_square_error << ' ' << rmse(force_constants, angles, quantum_mechanics_data_points) << '\n';
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