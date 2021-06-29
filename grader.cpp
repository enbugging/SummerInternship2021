#include <bits/stdc++.h>
using namespace std;

#include "ExtremaFinding.h"

string 
	quantum_mechanics_data_url = "abc";
int 
	number_of_terms = 3, 
	number_of_data_points = 36;
vector<int> 
	multiplicities = { 1, 2, 3, 4, 6};
vector<double> 
	angles, 
	quantum_mechanics_data_points,
	molecular_mechanics_data_points,
	force_constants;
double square_error;

void preprocess()
{
	// initialization
	angles.resize(number_of_data_points);
	quantum_mechanics_data_points.resize(number_of_data_points);
	square_error = 0;
	
	// read input
	ifstream quantum_mechanics_data_file;
	quantum_mechanics_data_file.open(quantum_mechanics_data_url);
	for(int i = 0; i < number_of_data_points; i++)
	{
		quantum_mechanics_data_file >> angles[i] >> quantum_mechanics_data_points[i];
	}
	quantum_mechanics_data_file.close();
}

double force_field_calculate(double angle)
{
	double E = 0;
	for(int i = 0; i < number_of_terms; i++)
	{
		E += force_constants[i] * (1 + cos(multiplicities[i] * angle));
	}
	return E;
}

void grading()
{
	// square error calculation
	for (int i = 0; i < number_of_data_points; i++)
	{
		double error = (force_field_calculate(angles[i]) - quantum_mechanics_data_points[i]);
		square_error += error * error;
	}

	// summary
	printf("SUMMARY\n");
	printf("----------------------------------\n");
	printf("Force constants: ");
	for(int i = 0; i < number_of_terms; i++)
	{
		printf("%lf ", force_constants[i]);
	}
	printf("\nSquare error: %lf", square_error);
}

int main()
{
	// preprocessing
	preprocess();

	// simulated annealing
	force_constants = simulatedAnnealing(number_of_terms, angles, quantum_mechanics_data_points);

	// grading
	grading();
}
