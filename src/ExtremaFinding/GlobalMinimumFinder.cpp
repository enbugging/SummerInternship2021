#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "GlobalMinimumFinder.h"
#ifndef _NR_
#define _NR_
#include "../numerical_recipes/nr.h" // numerical recipes
#endif
#define random_step(rng,min,max) min + (double)rng() / UINT_MAX * (max - min)

/*	-----------------------------------------------------------------------------------
 *	Utility functions, supporthing simulated-annealing-related functions
 * */
vector<int> multiplicities = {0, 1, 2, 3, 4, 6};
vector<double> 
    dummy_angles,
    dummy_quantum_mechanics_points, 
    dummy_quantum_mechanics_weights;
double dummy_threshold;

double force_field_calculate(
    vector<double>& force_constants,
    double angle)
{
    double E = force_constants[0];
    for (int i = 1; i < (int)force_constants.size(); i++)
    {
        E += force_constants[i] * (cos(multiplicities[i] * angle / 180.0 * M_PI));
    }
    return E;
}

double rmse(
    vector<double>& force_constants,
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    vector<double>& quantum_mechanics_weights)
{
    if (angles.empty()) return 0;
    // if the function is unweighted
    if (quantum_mechanics_weights.empty())
    {
        quantum_mechanics_weights.assign(angles.size(), 1.0);
    }
    double sum_of_squares_of_error = 0, sum_of_weights = 0;
    for (int i = 0; i < (int)angles.size(); i++)
    {
        double error = force_field_calculate(force_constants, angles[i]) - quantum_mechanics_points[i];
        sum_of_squares_of_error += quantum_mechanics_weights[i] * error * error;
        sum_of_weights += quantum_mechanics_weights[i];
    }

    sum_of_squares_of_error = sqrt(sum_of_squares_of_error / sum_of_weights);
    return sum_of_squares_of_error;
}

double coefficient(
    double t, 
    double threshold)
{
    double
        border_width = 0.00001, 
        A, 
        B, 
        extra_wall = 0.001, 
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
    vector<double>& quantum_mechanics_points, 
    vector<double>& quantum_mechanics_weights)
{
    if (angles.empty()) return 0;
    // if the function is unweighted
    if (quantum_mechanics_weights.empty())
    {
        quantum_mechanics_weights.assign(angles.size(), 1.0);
    }

    vector<double> new_force_constants = force_constants;
    double c = 0.01;
    for (int i = 0; i < (int) force_constants.size(); i++)
    {
        new_force_constants[i] *= coefficient(force_constants[i], threshold);
    }
    double sum_of_squares_of_error = 0, sum_of_weights = 0;
    for (int i = 0; i < (int) angles.size(); i++)
    {
        double error = force_field_calculate(new_force_constants, angles[i]) - quantum_mechanics_points[i];
        sum_of_squares_of_error += quantum_mechanics_weights[i] * error * error;
        sum_of_weights += quantum_mechanics_weights[i];
    }
    sum_of_squares_of_error = sqrt(sum_of_squares_of_error/sum_of_weights);
    for (int i = 0; i < (int) force_constants.size(); i++)
    {
        sum_of_squares_of_error += c * coefficient(force_constants[i], threshold);
    }
    return sum_of_squares_of_error;
}

DP rmse_with_threshold_wrapper(Vec_I_DP &x) {
    vector<double> y;
    y.resize(x.size());
    for (int i = 0; i < (int) x.size(); i++)
    {
        y[i] = x[i];
    }
    return rmse_with_threshold(dummy_threshold, y, dummy_angles, dummy_quantum_mechanics_points, dummy_quantum_mechanics_weights);
}

/*	-----------------------------------------------------------------------------------
 *	Simulated annealing and its variants. The implementations used Cauchy mutation scheme
 *  for best stability and power.
 * */
vector<double> simulated_annealing(
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double threshold, 
    vector<double>& quantum_mechanics_weights)
{
    // if the function is unweighted
    if (quantum_mechanics_weights.empty())
    {
        quantum_mechanics_weights.assign(angles.size(), 1.0);
    }

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
        double new_square_error = rmse_with_threshold(threshold, new_force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        if (new_square_error <= best_square_error || 
            random_step(rng, 0, 1) <= exp(-(new_square_error - best_square_error) / T))
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }
    }

    // further polishing
    // Powell's algorithm, using Numerical Recipes, provided by Prof. Alexandrov
    // setup minimizer
    /*
    dummy_angles = angles;
    dummy_quantum_mechanics_points = quantum_mechanics_points;
    dummy_quantum_mechanics_weights = quantum_mechanics_weights;
    dummy_threshold = threshold;
    DP FTOL = 1.0e-5;
    const int IMAXSTEP = 10000;
    int iter;
    DP fret;
    Vec_DP p(0.0, number_of_terms); // variables
    Mat_DP xi(number_of_terms, number_of_terms);
    for (int i = 0; i < number_of_terms; i++)
    {
        p[i] = force_constants[i];
    }
    for (int i = 0; i < number_of_terms; i++)
    {
        for (int j = 0; j < number_of_terms; j++)
        {
            xi[i][j] = (i == j ? force_constants[i] : 0.0);
        }
    }
    
    NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, rmse_with_threshold_wrapper);
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = p[i];
    }
    */

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
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double threshold, 
    vector<double>& quantum_mechanics_weights)
{
    // if the function is unweighted
    if (quantum_mechanics_weights.empty())
    {
        quantum_mechanics_weights.assign(angles.size(), 1.0);
    }

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
        double new_square_error = rmse_with_threshold(threshold, new_force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        if (new_square_error <= best_square_error + T)
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }
    }

    // further polishing
    // Powell's algorithm, using Numerical Recipes, provided by Prof. Alexandrov
    // setup minimizer
    /*
    dummy_angles = angles;
    dummy_quantum_mechanics_points = quantum_mechanics_points;
    dummy_quantum_mechanics_weights = quantum_mechanics_weights;
    dummy_threshold = threshold;
    DP FTOL = 1.0e-5;
    const int IMAXSTEP = 10000;
    int iter;
    DP fret;
    Vec_DP p(0.0, number_of_terms); // variables
    Mat_DP xi(number_of_terms, number_of_terms);
    for (int i = 0; i < number_of_terms; i++)
    {
        p[i] = force_constants[i];
    }
    for (int i = 0; i < number_of_terms; i++)
    {
        for (int j = 0; j < number_of_terms; j++)
        {
            xi[i][j] = (i == j ? force_constants[i] : 0.0);
        }
    }
    
    NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, rmse_with_threshold_wrapper);
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = p[i];
    }
    */

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
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double threshold, 
    vector<double>& quantum_mechanics_weights)
{
    // if the function is unweighted
    if (quantum_mechanics_weights.empty())
    {
        quantum_mechanics_weights.assign(angles.size(), 1.0);
    }

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
            double new_square_error = rmse_with_threshold(threshold, new_force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights);

            // if the new set of parameters is better, then we accept
            // else, we accept, with a probability corresponding to the temperature
            if (new_square_error <= best_square_error || 
                random_step(rng, 0, 1) <= exp(-(new_square_error - best_square_error) / T))
            {
                best_square_error = new_square_error,
                force_constants = new_force_constants;
            }
        }

        // further polishing
        // Powell's algorithm, using Numerical Recipes, provided by Prof. Alexandrov
        // setup minimizer
        /*
        dummy_angles = angles;
        dummy_quantum_mechanics_points = quantum_mechanics_points;
        dummy_quantum_mechanics_weights = quantum_mechanics_weights;
        dummy_threshold = threshold;
        DP FTOL = 1.0e-5;
        const int IMAXSTEP = 10000;
        int iter;
        DP fret;
        Vec_DP p(0.0, number_of_terms); // variables
        Mat_DP xi(number_of_terms, number_of_terms);
        for (int i = 0; i < number_of_terms; i++)
        {
            p[i] = force_constants[i];
        }
        for (int i = 0; i < number_of_terms; i++)
        {
            for (int j = 0; j < number_of_terms; j++)
            {
                xi[i][j] = (i == j ? force_constants[i] : 0.0);
            }
        }
        
        NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, rmse_with_threshold_wrapper);
        for (int i = 0; i < number_of_terms; i++)
        {
            force_constants[i] = p[i];
        }
        */

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
        cerr << "\nRMSE: " << best_square_error << ' ' << rmse(force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights) << '\n';
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