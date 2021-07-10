#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "src/ExtremaFinding/GlobalMinimumFinder.h"
#include "src/ExtremaFinding/LocalMinimaFinder.h"

string
	quantum_mechanics_url_points = "ref.dat", 
	quantum_mechanics_url_weights = "weight.dat";
int
	number_of_terms = 4,
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
	double r = rmse(force_constants, angles, quantum_mechanics_points);
	printf("\nSquare error: %lf\n", r);
}

double wrapper_function_for_finding_local_minima(
	vector<double> force_constants)
{
	return rmse(force_constants, angles, quantum_mechanics_points);
}

int main()
{
	// preprocessing
	preprocess();

	// GLOBAL MINIMUM FINDING
	// simulated annealing
	//force_constants = simulated_annealing(angles, quantum_mechanics_points, number_of_terms, 5000, 1.0, 3.0);
	//force_constants = threshold_accepting(angles, quantum_mechanics_points, number_of_terms, 5000, 50.0, 3.0);
	//force_constants = simulated_annealing_with_threshold(angles, quantum_mechanics_points, number_of_terms, 10000, 1.5, 3.0, 0.03);
		
	// grading
	//summary();

	/*
	double average = 0;
	vector<double> run;
	for (int i = 0; i < 10; i++)
	{
		force_constants = simulated_annealing_with_threshold(angles, quantum_mechanics_points, number_of_terms, 5000, 1.0, 3.0, 0.0);
		double sum_of_squares_of_error = rmse(force_constants, angles, quantum_mechanics_points);
		average += sum_of_squares_of_error;
		run.push_back(sum_of_squares_of_error);
	}
	double sd = 0;
	for (int i = 0; i < 10; i++)
	{
		sd += (run[i] - average) * (run[i] - average);
	}
	printf("Average: %lf\nStandard deviation: %lf", average/10, sqrt(sd/10));
	*/

	// LOCAL MINIMA FINDING
	result = MLSL(number_of_terms, 3.0, wrapper_function_for_finding_local_minima);\
	cerr << result.size() << '\n';
	for (int i = 0; i < (int) result.size(); i++)
	{
		cerr << "Minima " << i << " - RMSE: " << result[i].value << '\n';
		/*
		cerr << "Force constants: ";
		for (int j = 0; j < number_of_terms; j++)
		{
			cerr << result[i].x[j] << " ";
		}
		cerr << '\n';
		*/
	}
}
