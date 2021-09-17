#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "ObjectiveFunctions.h"

double rmse(
	vector<double>& set_of_force_constants, 
    vector<double>& energy, 
    vector<vector<double> >& angles, 
	int multiplicities[], 
	int number_of_data_points, 
    int number_of_angles, 
    int number_of_terms)
{
	double sum_of_square_error = 0;
	for (int i = 0; i < number_of_data_points; i++)
	{
		double error = energy[i] - set_of_force_constants[set_of_force_constants.size()-1];
		for (int j = 0; j < number_of_angles; j++)
		{
			double angle = angles[i][j];
			for (int k = 0; k < number_of_terms; k++)
			{
				error -= 
				set_of_force_constants[j * number_of_terms + k] * 
				cos(multiplicities[k] * angle / 180.0 * M_PI);
			}
		}
		sum_of_square_error += error * error;
	}
	return sqrt(sum_of_square_error/number_of_data_points);
}

double soft_sign(
	double x)
{
	double epsilon = 1e-4;
	if(abs(x) >= epsilon)
	{
		return (x > 0 ? 1 : -1);
	}
	else
	{
		return x/epsilon;
	}
}

double sign_of_main_multiplicity(
	vector<double>& set_of_force_constants, 
    int multiplicities[], 
    int number_of_angles, 
    int number_of_terms, 
    int main_mult)
{
	double res = 0;
	int idx_main_mult = -1;
	for (int i = 0; i < number_of_terms; i++)
	{
		if (multiplicities[i] == main_mult)
		{
			idx_main_mult = i;
			break;
		}
	}
	if (idx_main_mult != -1)
	{
		for (int i = 0; i < number_of_angles; i++)
		{
			for (int j = 0; j < i; j++)
			{
				double d = 
					soft_sign(
						set_of_force_constants[
							i * number_of_terms + 
							idx_main_mult]
						) - 
					soft_sign(
						set_of_force_constants[
							j * number_of_terms + 
							idx_main_mult]
						);
				res = max(res, abs(d));
			}
		}
	}
	return res;
}