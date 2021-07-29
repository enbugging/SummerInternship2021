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
vector<double> 
    dummy_angles,
    dummy_quantum_mechanics_points, 
    dummy_quantum_mechanics_weights;
int 
    multiplicities[6] = {0, 1, 2, 3, 4, 6},
    dummy_mask;
double 
    dummy_cutoff, 
    simplicity_accuracy_trading[6] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.35};

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

double cubical(double x)
{
    return -2 * x * x * x + 3 * x * x;
}

double non_analytic_smooth_function(double x)
{
    return x > 0 ? exp(-1/x) : 0.0;
}

double smooth_transition_function(double x)
{
    return non_analytic_smooth_function(x) / (non_analytic_smooth_function(x) + non_analytic_smooth_function(1 - x));
}

double coefficient(
    double t, 
    double cutoff)
{
    double 
        border_width = min(1e-4, cutoff/2), 
        center = cutoff - border_width;
    if (abs(t) >= cutoff) return 1;
    else if (abs(t) <= center) return 0;
    //else return (1.0 - cos(M_PI / border_width * (abs(t) - center)))/2.0;
    else return cubical((abs(t) - cutoff)/border_width + 1);
    //else return (abs(t) - cutoff)/border_width + 1;
    //else return smooth_transition_function((abs(t) - center) / border_width);
    //return 0.9999;
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
    for (int i = 0; i < (int) angles.size(); i++)
    {
        double error = force_field_calculate(force_constants, angles[i]) - quantum_mechanics_points[i];
        sum_of_squares_of_error += quantum_mechanics_weights[i] * error * error;
        sum_of_weights += quantum_mechanics_weights[i];
    }
    return sqrt(sum_of_squares_of_error/sum_of_weights);
}

double rmse_with_cutoff_and_simplicity_accuracy_trading(
    vector<double>& force_constants,
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    double cutoff, 
    vector<double>& quantum_mechanics_weights)
{
    if (angles.empty()) return 0;
    // if the function is unweighted
    if (quantum_mechanics_weights.empty())
    {
        quantum_mechanics_weights.assign(angles.size(), 1.0);
    }

    double sum_of_squares_of_error = 0, sum_of_weights = 0;
    for (int i = 0; i < (int) angles.size(); i++)
    {
        double error = force_field_calculate(force_constants, angles[i]) - quantum_mechanics_points[i];
        sum_of_squares_of_error += quantum_mechanics_weights[i] * error * error;
        sum_of_weights += quantum_mechanics_weights[i];
    }
    sum_of_squares_of_error = sqrt(sum_of_squares_of_error/sum_of_weights);
    for (int i = 0; i < (int) force_constants.size(); i++)
    {
        sum_of_squares_of_error += simplicity_accuracy_trading[i] * coefficient(force_constants[i], cutoff);
    }
    return sum_of_squares_of_error;
}

DP rmse_with_cutoff_and_simplicity_accuracy_trading_wrapper(Vec_I_DP &x) {
    vector<double> y;
    y.resize(x.size());
    for (int i = 0; i < (int) x.size(); i++)
    {
        y[i] = x[i];
    }
    return rmse_with_cutoff_and_simplicity_accuracy_trading(y, dummy_angles, dummy_quantum_mechanics_points, dummy_cutoff, dummy_quantum_mechanics_weights);
}

DP rmse_wrapper(Vec_I_DP &x) {
    vector<double> y;
    y.resize(x.size());
    for (int i = 0; i < (int) x.size(); i++)
    {
        if((dummy_mask >> i) & 1)
        {
            y[i] = x[i];
        }
    }
    return rmse(y, dummy_angles, dummy_quantum_mechanics_points, dummy_quantum_mechanics_weights);
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
    double cutoff, 
    vector<double>& quantum_mechanics_weights)
{
    // if the function is unweighted
    if (quantum_mechanics_weights.empty())
    {
        quantum_mechanics_weights.assign(angles.size(), 1.0);
    }

    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count(); // 17 for testing
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
        double new_square_error = rmse(new_force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        if (new_square_error <= best_square_error || 
            random_step(rng, 0, 1) <= exp(-(new_square_error - best_square_error) / T))
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }
        /*
        for (int i = 0; i < number_of_terms; i++)
        {
            cerr << force_constants[i] << ' ';
        }
        cerr << " - Error: " << best_square_error << '\n';
        */
    }

    // further polishing
    // Powell's algorithm, using Numerical Recipes, provided by Prof. Alexandrov
    // setup minimizer
    ///*
    dummy_mask = (1<<number_of_terms) - 1;
    dummy_angles = angles;
    dummy_quantum_mechanics_points = quantum_mechanics_points;
    dummy_quantum_mechanics_weights = quantum_mechanics_weights;
    dummy_cutoff = cutoff;
    DP FTOL = 1.0e-6;
    const int IMAXSTEP = 100000;
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
    
    NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, rmse_wrapper);
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = p[i];
    }
    //*/
    
    /*
    for (int i = 0; i < number_of_terms; i++)
    {
        cerr << force_constants[i] << ' ';
    }
    cerr << " - Error: " << best_square_error << '\n';
    */
    // setting trapped coefficients to be zero
    for (int i = 0; i < number_of_terms; i++)
    {
        //cerr << abs(force_constants[i]) - cutoff << '\n';
        if (abs(force_constants[i]) <= cutoff)
        {
            force_constants[i] = 0;
        }
    }

    /*
    for (int i = 0; i < number_of_terms; i++)
    {
        cerr << force_constants[i] << ' ';
    }
    cerr << " - Error: " << best_square_error << '\n';
    */
    return force_constants;
}

vector<double> simulated_annealing_with_simplicity_accuracy_trading(
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double cutoff, 
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
        double new_square_error = rmse_with_cutoff_and_simplicity_accuracy_trading(new_force_constants, angles, quantum_mechanics_points, cutoff, quantum_mechanics_weights);

        // if the new set of parameters is better, then we accept
        // else, we accept, with a probability corresponding to the temperature
        if (new_square_error <= best_square_error || 
            random_step(rng, 0, 1) <= exp(-(new_square_error - best_square_error) / T))
        {
            best_square_error = new_square_error,
            force_constants = new_force_constants;
        }
            
        for (int i = 0; i < number_of_terms; i++)
        {
            cerr << force_constants[i] << ' ';
        }
        cerr << " - Error: " << best_square_error << ' ' << rmse(force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights) << '\n';
    }

    // further polishing
    // Powell's algorithm, using Numerical Recipes, provided by Prof. Alexandrov
    // setup minimizer
    ///*
    dummy_angles = angles;
    dummy_quantum_mechanics_points = quantum_mechanics_points;
    dummy_quantum_mechanics_weights = quantum_mechanics_weights;
    dummy_cutoff = cutoff;
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
    
    NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, rmse_with_cutoff_and_simplicity_accuracy_trading_wrapper);
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = p[i];
    }
    //*/
    
    for (int i = 0; i < number_of_terms; i++)
    {
        cerr << force_constants[i] << ' ';
    }
    cerr << " - Error: " << best_square_error << ' ' << rmse(force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights) << '\n';
    // setting trapped coefficients to be zero
    for (int i = 0; i < number_of_terms; i++)
    {
        if (abs(force_constants[i]) <= cutoff)
        {
            force_constants[i] = 0;
        }
    }
    
    for (int i = 0; i < number_of_terms; i++)
    {
        cerr << force_constants[i] << ' ';
    }
    cerr << " - Error: " << best_square_error << ' ' << rmse(force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights) << '\n';
    
    return force_constants;
}

vector<double> threshold_accepting(
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double cutoff, 
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
        double new_square_error = rmse(new_force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights);

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
    ///*
    dummy_mask = (1<<number_of_terms) - 1;
    dummy_angles = angles;
    dummy_quantum_mechanics_points = quantum_mechanics_points;
    dummy_quantum_mechanics_weights = quantum_mechanics_weights;
    dummy_cutoff = cutoff;
    DP FTOL = 1.0e-6;
    const int IMAXSTEP = 100000;
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
    
    NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, rmse_wrapper);
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = p[i];
    }
    //*/

    // setting trapped coefficients to be zero
    for (int i = 0; i < number_of_terms; i++)
    {
        if (abs(force_constants[i]) <= cutoff)
        {
            force_constants[i] = 0;
        }
    }

    return force_constants;
}

vector<double> simulated_annealing_brute_force(
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double cutoff, 
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
    double optimal_square_error = numeric_limits<double>::max();

    for (int mask = 0; mask < (1<<number_of_terms); mask++)
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
                if((mask >> j) & 1)
                {
                    new_force_constants[j] += T * tan(random_step(rng, -M_PI/2, M_PI/2));
                }
            }
            double new_square_error = rmse(new_force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights);

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
        ///*
        dummy_mask = mask;
        dummy_angles = angles;
        dummy_quantum_mechanics_points = quantum_mechanics_points;
        dummy_quantum_mechanics_weights = quantum_mechanics_weights;
        dummy_cutoff = cutoff;
        DP FTOL = 1.0e-6;
        const int IMAXSTEP = 100000;
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
        
        NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, rmse_wrapper);
        for (int i = 0; i < number_of_terms; i++)
        {
            force_constants[i] = p[i];
        }
        //*/
        
        best_square_error = rmse(force_constants, angles, quantum_mechanics_points, quantum_mechanics_weights);
        for (int i = 0; i < number_of_terms; i++)
        {
            if (abs(force_constants[i]) <= cutoff)
            {
                force_constants[i] = 0;
            }
            else
            {
                best_square_error += simplicity_accuracy_trading[i];
            }
        }

        if (optimal_square_error > best_square_error)
        //if (mask == (1<<number_of_terms) - 1)
        {
            optimal_square_error = best_square_error;
            best_force_constants = force_constants;
        }

        ///*
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
        //*/
    }

    return best_force_constants;
}