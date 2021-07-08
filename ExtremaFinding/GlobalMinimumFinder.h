#include <bits/stdc++.h>
using namespace std;

double force_field_calculate(
    vector<double>& force_constants,
    double angle
);

double rmse(
    vector<double>& force_constants,
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