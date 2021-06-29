#include <bits/stdc++.h>
using namespace std;

#include "ExtremaFinding.h"

string
	quantum_mechanics_data_url = "qm.dat";
int
	number_of_terms = 3,
	number_of_data_points = 36;
vector<double>
	angles,
	quantum_mechanics_data_points,
	molecular_mechanics_data_points,
	force_constants;

void preprocess()
{
	// initialization
	angles.resize(number_of_data_points);
	quantum_mechanics_data_points.resize(number_of_data_points);

	// read input
	ifstream quantum_mechanics_data_file;
	quantum_mechanics_data_file.open(quantum_mechanics_data_url);
	for (int i = 0; i < number_of_data_points; i++)
	{
		quantum_mechanics_data_file >> angles[i] >> quantum_mechanics_data_points[i];
	}
	quantum_mechanics_data_file.close();
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
	double sum_of_squares_of_error = square_error(force_constants, angles, quantum_mechanics_data_points);
	printf("\nSquare error: %lf", sqrt(sum_of_squares_of_error/number_of_data_points));
}

int main()
{
	// preprocessing
	preprocess();

	/*
	// simulated annealing
	force_constants = simulatedAnnealing(50000, 0.5, 1, 0.3, number_of_terms, angles, quantum_mechanics_data_points);

	// grading
	summary();
	*/
	double average = 0;
	vector<double> run;
	for (int i = 0; i < 10; i++)
	{
		force_constants = simulatedAnnealing(5000, 1000, 1, 0.5, number_of_terms, angles, quantum_mechanics_data_points);
		double sum_of_squares_of_error = square_error(force_constants, angles, quantum_mechanics_data_points);
		average += sqrt(sum_of_squares_of_error/number_of_data_points);
		run.push_back(sqrt(sum_of_squares_of_error/number_of_data_points));
	}
	double sd = 0;
	for (int i = 0; i < 10; i++)
	{
		sd += (run[i] - average) * (run[i] - average);
	}
	printf("Average: %lf\nStandard deviation: %lf", average/10, sqrt(sd/10));
}
