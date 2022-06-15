#include <bits/stdc++.h>
using namespace std;

/* 
Calculate the root mean square error of reference data
and the prediction made by force constants and
offset constant contained in set_of_force_constant.
 */
double rmse(
	vector<double>& set_of_force_constants, 
    vector<double>& energy, 
    vector<vector<double> >& angles, 
	int multiplicities[], 
	int number_of_data_points, 
    int number_of_angles, 
    int number_of_terms);

/* 
Calculate the root mean square (normalized) of the
magnitude of difference between force constants of 
main multiplicites.
 */
double magnitude_main_multiplicity(
	vector<double>& set_of_force_constants, 
    int multiplicities[], 
    int number_of_angles, 
    int number_of_terms, 
    int main_mult);

/* 
Calculate the root mean square (normalized) of the
magnitude of difference in soft sign of force constants 
of main multiplicites.
 */
double sign_of_main_multiplicity(
	vector<double>& set_of_force_constants, 
    int multiplicities[], 
    int number_of_angles, 
    int number_of_terms, 
    int main_mult);