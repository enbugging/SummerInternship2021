#include <bits/stdc++.h>
using namespace std;

double rmse(
	vector<double>& set_of_force_constants, 
    vector<double>& energy, 
    vector<vector<double> >& angles, 
	int multiplicities[], 
	int number_of_data_points, 
    int number_of_angles, 
    int number_of_terms);

double sign_of_main_multiplicity(
	vector<double>& set_of_force_constants, 
    int multiplicities[], 
    int number_of_angles, 
    int number_of_terms, 
    int main_mult);