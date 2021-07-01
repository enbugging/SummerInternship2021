#include <bits/stdc++.h>
using namespace std;

double force_field_calculate(
    vector<double> force_constants,
    double angle);

double square_error(
    vector<double> force_constant,
    vector<double> angles,
    vector<double> quantum_mechanics_data_points);

vector<double> simulated_annealing(
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double radius, 
    int number_of_terms,
    vector<double> angles,
    vector<double> quantum_mechanics_data_points);

vector<double> threshold_accepting(
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double radius, 
    int number_of_terms,
    vector<double> angles,
    vector<double> quantum_mechanics_data_points);

vector<double> simulated_annealing_with_threshold(
    double threshold, 
    int number_of_steps, 
    double initial_temperature, 
    double initial_radius, 
    double radius, 
    int number_of_terms,
    vector<double> angles,
    vector<double> quantum_mechanics_data_points);