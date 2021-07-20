#include <bits/stdc++.h>
using namespace std;

double force_field_calculate(
    vector<double>& force_constants,
    double angle
);

static vector<double> dummy;

double rmse(
    vector<double>& force_constants,
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    vector<double>& quantum_mechanics_weights = dummy);

double rmse_with_cutoff(
    double cutoff, 
    vector<double>& force_constants,
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    vector<double>& quantum_mechanics_weights = dummy);

vector<double> simulated_annealing(
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double cutoff = 0, 
    vector<double>& quantum_mechanics_weights = dummy);

vector<double> threshold_accepting(
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double cutoff = 0, 
    vector<double>& quantum_mechanics_weights = dummy);

vector<double> simulated_annealing_with_cutoff(
    vector<double>& angles,
    vector<double>& quantum_mechanics_points, 
    int number_of_terms,
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double cutoff = 0, 
    vector<double>& quantum_mechanics_weights = dummy);