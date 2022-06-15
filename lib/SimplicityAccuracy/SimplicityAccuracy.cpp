#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "SimplicityAccuracy.h"

double 
    simplicity_accuracy_trading[5] = {0.01, 0.015, 0.02, 0.025, 0.035};

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
    return 
		non_analytic_smooth_function(x) / 
		(
			non_analytic_smooth_function(x) + 
			non_analytic_smooth_function(1 - x)
		);
}

double coefficient(
    double t, 
    double cutoff)
{
    double 
        border_width = 1e-2, 
        center = cutoff + border_width;
    if (abs(t) >= center) return 1;
    else if (abs(t) <= cutoff) return 0;
    //else return (1.0 - cos(M_PI / border_width * (abs(t) - cutoff)))/2.0;
    else return cubical((abs(t) - cutoff)/border_width);
    //else return (abs(t) - cutoff)/border_width;
    //else return smooth_transition_function((abs(t) - cutoff) / border_width);
    //return 0.9999;
}

double simplicity_vs_accuracy(
    vector<double>& set_of_force_constants, 
    int multiplicities[], 
    int number_of_angles, 
    int number_of_terms, 
    int main_mult, 
	double cutoff)
{
    double error = 0;
	for (int i = 0; i < number_of_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] == main_mult)
			{
				error += 
				simplicity_accuracy_trading[j] * 
				coefficient(
					set_of_force_constants[i * number_of_terms + j], 
					cutoff);
			}
		}
	}
    return error;
}