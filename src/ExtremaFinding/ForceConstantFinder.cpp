#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "GlobalMinimumFinder.h"
#ifndef _NR_
#define _NR_
#include "../numerical_recipes/nr.h" // numerical recipes
#endif


double optimal_force_constant(
    vector<double>& angles, 
    vector<double>& quantum_mechanics_points, 
    int harmonic)
{
    double ans = 0;
    for (int i = 0; i < (int) angles.size(); i++)
    {
        ans += quantum_mechanics_points[i] * cos(angles[i] * M_PI / 180.0 * harmonic);
    }
    return (harmonic ? 2 * ans : ans) / (int) angles.size();
}

vector<double> find_with_simplicity_accuracy_trading(
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    int number_of_terms, 
    double cutoff = 0)
{
    // initialization
    int multiplicities[6] = {0, 1, 2, 3, 4, 6};
    double simplicity_accuracy_trading[6] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.35};
    vector<double> force_constants, threshold;
    vector<pair<double, int> > auxiliary;
    force_constants.resize(number_of_terms);
    threshold.resize(number_of_terms);
    auxiliary.resize(number_of_terms);

    double sum_of_squares_of_energy = 0;
    for (int i = 0; i < (int) angles.size(); i++)
    {
        sum_of_squares_of_energy += quantum_mechanics_points[i] * quantum_mechanics_points[i];
    }

    // calculate force constants directly as Cosine Transform coefficients
    for (int i = 0; i < number_of_terms; i++)
    {
        force_constants[i] = optimal_force_constant(angles, quantum_mechanics_points, multiplicities[i]);
        double error_contributor = (i ? angles.size() / 2.0 : angles.size()) * force_constants[i] * force_constants[i];
        double scaled_sat = simplicity_accuracy_trading[i] * sqrt(angles.size());
        if (error_contributor < scaled_sat * scaled_sat)
        {
            threshold[i] = 0;
        }
        else
        {
            threshold[i] = error_contributor + 
                        (error_contributor - 
                            scaled_sat * 
                            scaled_sat)/
                        (2.0 * scaled_sat) * 
                        (error_contributor - 
                            scaled_sat * 
                            scaled_sat)/
                        (2.0 * scaled_sat);
        }
        auxiliary[i] = {threshold[i], i};
    }
    sort(auxiliary.begin(), auxiliary.end(), greater<pair<double, int> >());
    for (int i = 0; i < number_of_terms; i++)
    {
        int idx = auxiliary[i].second;
        if (sum_of_squares_of_energy <= auxiliary[i].first)
        {
            sum_of_squares_of_energy -= (idx ? angles.size() / 2.0 : angles.size()) * force_constants[idx] * force_constants[idx];
        }
        else
        {
            force_constants[idx] = 0;
        }
    }

    return force_constants;
}