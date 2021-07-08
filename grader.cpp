#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "ExtremaFinding/GlobalMinimumFinder.h"
//#include "ExtremaFinding/LocalMinimaFinder.h"

string
	quantum_mechanics_data_url = "ref.dat";
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
	molecular_mechanics_data_points.resize(number_of_data_points);

	// read input
	ifstream quantum_mechanics_data_file;

	vector<double> charmm;
	quantum_mechanics_data_file.open(quantum_mechanics_data_url);
	for (int i = 0; i < number_of_data_points; i++)
	{
		quantum_mechanics_data_file >> angles[i] >> quantum_mechanics_data_points[i];
	}
	quantum_mechanics_data_file.close();
}





/*
double test_rmse(
    vector<double>& force_constants)
{
    if (angles.empty()) return 0;
    double sum_of_squares_of_error = 0;
    for (int i = 0; i < (int)angles.size(); i++)
    {
        double error = force_field_calculate(force_constants, angles[i]) - quantum_mechanics_data_points[i];
        sum_of_squares_of_error += error * error;
    }
    sum_of_squares_of_error = sqrt(sum_of_squares_of_error/angles.size());
    return sum_of_squares_of_error;
}
vector<Point> result = MLSL_MADS(number_of_terms, 3.0, test_rmse);
*/











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
	double r = rmse(force_constants, angles, quantum_mechanics_data_points);
	printf("\nSquare error: %lf\n", r);
}

int main()
{
	// preprocessing
	preprocess();

	
	// simulated annealing
	//force_constants = global_minimum_finder.simulated_annealing(angles, quantum_mechanics_data_points, number_of_terms, 5000, 1.0, 3.0);
	//force_constants = global_minimum_finder.threshold_accepting(angles, quantum_mechanics_data_points, number_of_terms, 5000, 50.0, 3.0);
	//force_constants = global_minimum_finder.simulated_annealing_with_threshold(angles, quantum_mechanics_data_points, number_of_terms, 5000, 1.0, 3.0, 0.1);
		
	// grading
	//summary();

	///*
	double average = 0;
	vector<double> run;
	for (int i = 0; i < 10; i++)
	{
		force_constants = simulated_annealing_with_threshold(angles, quantum_mechanics_data_points, number_of_terms, 5000, 1.0, 3.0, 0.0);
		double sum_of_squares_of_error = rmse(force_constants, angles, quantum_mechanics_data_points);
		average += sum_of_squares_of_error;
		run.push_back(sum_of_squares_of_error);
	}
	double sd = 0;
	for (int i = 0; i < 10; i++)
	{
		sd += (run[i] - average) * (run[i] - average);
	}
	printf("Average: %lf\nStandard deviation: %lf", average/10, sqrt(sd/10));
	//*/
}
