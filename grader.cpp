#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "src/ExtremaFinding/GlobalMinimumFinder.h"

string
	quantum_mechanics_url_points = "ref.dat", 
	quantum_mechanics_url_weights = "weight.dat";
int
	number_of_terms = 6,
	number_of_data_points = 36;
vector<double>
	angles,
	quantum_mechanics_points,
	quantum_mechanics_weights,
	force_constants;
vector<Point> result;

void preprocess()
{
	// initialization
	angles.resize(number_of_data_points);
	quantum_mechanics_points.resize(number_of_data_points);
	quantum_mechanics_weights.resize(number_of_data_points);

	// read input
	ifstream quantum_mechanics_file;

	quantum_mechanics_file.open(quantum_mechanics_url_points);
	for (int i = 0; i < number_of_data_points; i++)
	{
		quantum_mechanics_file >> angles[i] >> quantum_mechanics_points[i];
	}
	quantum_mechanics_file.close();
	quantum_mechanics_file.open(quantum_mechanics_url_weights);
	for (int i = 0; i < number_of_data_points; i++)
	{
		quantum_mechanics_file >> angles[i] >> quantum_mechanics_weights[i];
	}
	quantum_mechanics_file.close();
}


double rmse(
    vector<double>& force_constants)
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

void summary()
{
	// summary
	printf("SUMMARY\n");
	printf("---------------------------------------------\n");
	printf("Force constants: ");
	for (int i = 0; i < number_of_terms; i++)
	{
		printf("%lf ", force_constants[i]);
	}
	double r = rmse(force_constants);
	printf("\nSquare error: %lf\n", r);
}

int main()
{
	// preprocessing
	preprocess();

	// GLOBAL MINIMUM FINDING
	// simulated annealing
    force_constants = simulated_annealing(rmse, number_of_terms, 3.0);
	//force_constants = threshold_accepting(angles, quantum_mechanics_points, number_of_terms, 5000, 50.0, 3.0);
		
	// grading
	summary();
}
