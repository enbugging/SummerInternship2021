#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "../../src/ExtremaFinding/GlobalMinimumFinder.h"
#include "../../src/ObjectiveFunctions/ObjectiveFunctions.h"
#include "../../src/SimplicityAccuracy/SimplicityAccuracy.h"

//string quantum_mechanic_data = "ethane_dihe_c1_c2.dat";
//string quantum_mechanic_data = "propane_dihe_c1_c2.dat";
//string quantum_mechanic_data = "butane_dihe_c2_c3.dat";
const int MAX_LENGTH = 10;
int
	number_of_angles, 
	number_of_terms = 3,
	number_of_data_points = 36, 
	number_of_distinct_angles, 
	nbrs_of_central_atom_1,
	nbrs_of_central_atom_2, 
	main_mult, 
	angles_id[9], 
	multiplicities[5] = {1, 2, 3, 4, 6};
double
	epsilon_main, // error by optimization of principle force constants only, named by Prof. Alexandrov
	epsilon_0 = 0.001, 
	w_total;
vector<vector<double> > 
	angles, 
	interaction, 
	pearson_interact_vs_angles, 
	weights;
vector<double> 
	energy,
	magnitude_of_interaction, 
	set_of_force_constants_without_offset_constant;
map<string, int> 
	angles_id_dict;
map<int, string>
	distinct_angles_represent;
ofstream log_file;
	
/* 
--------------------------------------------------------------------------------------
Helper functions
 */
/* 
Principle multiplicity, taken from 
QuickFF: toward a generally applicable methodology
to quickly derive accurate force fields for
Metal-Organic Frameworks from ab initio input
 */
int main_multiplicity(
	int a, 
	int b)
{
	if (a < b) swap(a, b);
	if (a == 4 && b == 4) return 3;
	if (a == 4 && b == 3) return 6;
	if (a == 4 && b == 2) return 3;
	if (a == 3 && b == 3) return 2;
	if (a == 3 && b == 2) return 2;
	if (a == 2 && b == 2) return 2;
	return -1;
}

/* 
Cheking if a string is a number
 */
bool isNumber(const string& str)
{
	if (str.length() == 0) return false;
    for (char const &c : str) {
        if (std::isdigit(c) == 0) return false;
    }
    return true;
}

/* 
--------------------------------------------------------------------------------------
Preprocess functions
 */
void input(string quantum_mechanic_data)
{
	// initialization
	angles.resize(number_of_data_points);
	interaction.resize(number_of_data_points);
	energy.resize(number_of_data_points);

	// initialize input stream
	ifstream quantum_mechanics_file;
	string dummy;
	quantum_mechanics_file.open(quantum_mechanic_data);
	
	// read number of neighbors
	getline(quantum_mechanics_file, dummy);
	quantum_mechanics_file >> nbrs_of_central_atom_1 >> nbrs_of_central_atom_2;
	main_mult = 
		main_multiplicity(
			nbrs_of_central_atom_1, 
			nbrs_of_central_atom_2);
	
	// read number of angles
	while(not isNumber(dummy))
	{
		quantum_mechanics_file >> dummy;
	}
	number_of_angles = 0;
	for (int i = 0; i < (int) dummy.length(); i++)
	{
		number_of_angles *= 10;
		number_of_angles += (dummy[i] - '0');
	}
	getline(quantum_mechanics_file, dummy);
	for (int i = 0; i < number_of_data_points; i++)
	{
		angles[i].resize(number_of_angles);
		interaction[i].resize(number_of_angles);
	}

	// read angle descrition	
	for (int i = 0; i < number_of_angles; i++)
	{
		// angle representation
		vector<string> s(4);
		for (int j = 0; j < 6; j++)
		{
			quantum_mechanics_file >> dummy;
		}
		quantum_mechanics_file >> s[0] >> s[1] >> s[2] >> s[3];
		getline(quantum_mechanics_file, dummy);
		string 
			t1 = s[0] + ' ' + s[1] + ' ' + s[2] + ' ' + s[3],
			t2 = s[3] + ' ' + s[2] + ' ' + s[1] + ' ' + s[0];
		
		cerr << t1 << ' ' << t2 << '\n';
		// get id of the angle
		if (angles_id_dict.find(t1) != angles_id_dict.end())
		{
			angles_id[i] = angles_id_dict[t1];
		}
		else if (angles_id_dict.find(t2) != angles_id_dict.end())
		{
			angles_id[i] = angles_id_dict[t2];
		}
		else
		{
			angles_id_dict[t1] = number_of_distinct_angles;
			angles_id_dict[t2] = number_of_distinct_angles;
			angles_id[i] = number_of_distinct_angles;
			distinct_angles_represent[number_of_distinct_angles] = t1;
			number_of_distinct_angles++;
		}
	}

	// read actual data
	getline(quantum_mechanics_file, dummy);
	for (int i = 0; i < number_of_data_points; i++)
	{
		for (int j = 0; j < number_of_angles; j++)
		{
			quantum_mechanics_file >> angles[i][j] >> interaction[i][j];
		}
		quantum_mechanics_file >> energy[i];
	}
	quantum_mechanics_file.close();
}

void correlation_preparation()
{
	pearson_interact_vs_angles.resize(number_of_angles);
	magnitude_of_interaction.resize(number_of_angles);
	for (int i = 0; i < number_of_angles; i++)
	{
		pearson_interact_vs_angles[i].resize(number_of_terms);
		vector<double> a(number_of_data_points);
		double 
			sumA = 0, 
			sumB = 0, 
			sumAB = 0, 
			sumA2 = 0, 
			sumB2 = 0, 
			max_a = numeric_limits<double>::min(), 
			min_a = numeric_limits<double>::max();
		for (int j = 0; j < number_of_data_points; j++)
		{
			a[j] = interaction[j][i];
			sumA += a[j];
			sumA2 += a[j] * a[j];
			min_a = min(min_a, a[j]);
			max_a = max(max_a, a[j]);
		}
		magnitude_of_interaction[i] = max_a - min_a;
		
		// 2 sequences of interaction and cos(multiplicity * angle)
		for (int j = 0; j < number_of_terms; j++)
		{
			sumB = 0, sumB2 = 0, sumAB = 0;
			for (int k = 0; k < number_of_data_points; k++)
			{
				double b = cos(multiplicities[j] * angles[k][i] * M_PI / 180.0);
				sumB += b;
				sumB2 += b * b;
				sumAB += a[k] * b;
			}

			double p = 
				(number_of_data_points * sumAB - sumA * sumB) / 
				sqrt(
					(number_of_data_points * sumA2 - sumA * sumA) * 
					(number_of_data_points * sumB2 - sumB * sumB)
				);
			
			pearson_interact_vs_angles[i][j] = p;
		}
	}
}

double pre_objective_function(
	vector<double>& main_force_constants)
{
	vector<double> set_of_force_constants(number_of_angles * number_of_terms + 1);
	set_of_force_constants[set_of_force_constants.size() - 1] = 
	main_force_constants[main_force_constants.size() - 1];
	double result = 0;
	for (int i = 0; i < number_of_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] == main_mult)
			{	
				set_of_force_constants[number_of_terms * i + j] = 
				main_force_constants[i];
			}
		}
	}
	double r = rmse(
		set_of_force_constants, 
		energy, 
		angles, 
		multiplicities, 
		number_of_data_points, 
		number_of_angles, 
		number_of_terms);
	result += r;
	result += 0.1 * 
		magnitude_main_multiplicity(
			set_of_force_constants, 
			multiplicities, 
			number_of_angles,
			number_of_terms,
			main_mult);
	result += 1000 * exp(-r) * sign_of_main_multiplicity(
			set_of_force_constants, 
			multiplicities, 
			number_of_angles,
			number_of_terms,
			main_mult);
	return result;
}

double objective_function(
	vector<double>& set_of_force_constants_reduced)
{
	// expand the set of force constants
	vector<double> set_of_force_constants(number_of_angles * number_of_terms + 1);
	set_of_force_constants[set_of_force_constants.size() - 1] = 
	set_of_force_constants_reduced[set_of_force_constants_reduced.size() - 1];
	double result = 0;
	for (int i = 0; i < number_of_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			set_of_force_constants[number_of_terms * i + j] = 
			set_of_force_constants_reduced[number_of_terms * angles_id[i] + j];
		}
	}
	double r = rmse(
		set_of_force_constants, 
		energy, 
		angles, 
		multiplicities, 
		number_of_data_points, 
		number_of_angles, 
		number_of_terms);
	result += r;
	result += 0.1 * magnitude_main_multiplicity(
			set_of_force_constants, 
			multiplicities, 
			number_of_angles,
			number_of_terms,
			main_mult);
	result += 1000 * exp(-r) * sign_of_main_multiplicity(
			set_of_force_constants, 
			multiplicities, 
			number_of_angles,
			number_of_terms,
			main_mult);
	r = 0;

	for (int i = 0; i < number_of_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] != main_mult)
			{
				r += 
				weights[i][j]
				 * set_of_force_constants[number_of_terms * i + j] 
				 * set_of_force_constants[number_of_terms * i + j];
			}
		}
	}
	result += w_total * sqrt(r);
	//result += simplicity_vs_accuracy(set_of_force_constants);
	return result;
}

double RMSE_with_respect_to_offset_constant(
	vector<double>& offset_constant)
{
	vector<double> set_of_force_constants(number_of_angles * 6 + 1);
	for (int i = 0; i < (int) set_of_force_constants_without_offset_constant.size(); i++)
	{
		set_of_force_constants[i] = set_of_force_constants_without_offset_constant[i];
	}
	set_of_force_constants[set_of_force_constants.size() - 1] = offset_constant[0];
	return rmse(set_of_force_constants, 
				energy, 
				angles, 
				multiplicities, 
				number_of_data_points, 
				number_of_angles, 
				6);
}


/* 
Main function
 */
vector<double> optimization_scheme()
{
	//cerr << "1\n";
	vector<double> main_force_constants = simulated_annealing(pre_objective_function, number_of_angles + 1, 1.0, 100000);
	vector<double> set_of_force_constants(number_of_angles * number_of_terms + 1);
	set_of_force_constants[set_of_force_constants.size() - 1] = 
	main_force_constants[main_force_constants.size() - 1];
	for (int j = 0; j < number_of_angles; j++)
	{
		for (int k = 0; k < number_of_terms; k++)
		{
			if (multiplicities[k] == main_mult)
			{	
				set_of_force_constants[number_of_terms * j + k] = 
				main_force_constants[j];
				//cerr << "Main force constant " << j << ": " << main_force_constants[j] << '\n';
			}
		}
	}
	epsilon_main = rmse(
		set_of_force_constants, 
		energy, 
		angles, 
		multiplicities, 
		number_of_data_points, 
		number_of_angles, 
		number_of_terms);
	
	//cerr << "2\n";
	double sum_of_energy_square = 0;
	for (int j = 0; j < number_of_data_points; j++)
	{
		sum_of_energy_square += energy[j] * energy[j];
	}
	w_total = 0.01 * sqrt(sum_of_energy_square)/(epsilon_main);
	
	weights.resize(number_of_angles);
	double sum_of_weights = 0;
	for (int j = 0; j < number_of_angles; j++)
	{
		weights[j].resize(number_of_terms);
		for (int k = 0; k < number_of_terms; k++)
		{
			weights[j][k] = 
				1.0 / (
					abs(pearson_interact_vs_angles[j][k])
					* magnitude_of_interaction[j] 
					+ epsilon_0);
			sum_of_weights += weights[j][k];
		}
	}
	// normalization of weights
	for (int j = 0; j < number_of_angles; j++)
	{
		for (int k = 0; k < number_of_terms; k++)
		{
			weights[j][k] /= sum_of_weights;
		}
	}

	//cerr << "3\n";
	vector<double> set_of_force_constants_reduced = simulated_annealing(objective_function, number_of_distinct_angles * number_of_terms + 1, 1.0);
	return set_of_force_constants_reduced;
}

void output(
	ostream& outstream, 
	int& i, 
	double& sum_error, 
	double& prev_mean,
	double& current_mean, 
	double& M, 
	vector<double>& set_of_force_constants_reduced)
{
	//cerr << "4\n";
	// expand the set of force constants
	vector<double> set_of_force_constants(number_of_angles * number_of_terms + 1);
	set_of_force_constants[set_of_force_constants.size() - 1] = 
	set_of_force_constants_reduced[set_of_force_constants_reduced.size() - 1];
	for (int j = 0; j < number_of_angles; j++)
	{
		for (int k = 0; k < number_of_terms; k++)
		{
			set_of_force_constants[number_of_terms * j + k] = 
			set_of_force_constants_reduced[number_of_terms * angles_id[j] + k];
		}
	}
	
	//cerr << "5\n";
	double error = rmse(
		set_of_force_constants, 
		energy, 
		angles, 
		multiplicities, 
		number_of_data_points, 
		number_of_angles, 
		number_of_terms);
	sum_error += error;
	prev_mean = current_mean;
	current_mean = (prev_mean * (i - 1) + error)/i;
	M += (error - current_mean) * (error - prev_mean);
	
	//cerr << "6\n";
	outstream << "Error: " 
			<< error << '\n';
	outstream << "Offset constant: " << set_of_force_constants_reduced[set_of_force_constants_reduced.size() - 1] << '\n';
	outstream << "Force constants of type:\n";
	for (int j = 0; j < number_of_distinct_angles; j++)
	{
		outstream << distinct_angles_represent[j] << ": ";
		for (int k = 0; k < number_of_terms; k++)
		{
			string str = to_string(set_of_force_constants_reduced[j * number_of_terms + k]);
			for (int l = 0; l < MAX_LENGTH - (int) str.length(); l++)
			{
				outstream << ' ';
			}
			outstream << str;
		}
		outstream << '\n';
	}

	//cerr << "7\n";
	outstream << "Energy profile comparison\n";
	outstream << " Reference|   Model|\n";

	for (int i = 0; i < number_of_data_points; i++)
	{
		double energy_prediction = set_of_force_constants[set_of_force_constants.size()-1];
		for (int j = 0; j < number_of_angles; j++)
		{
			double angle = angles[i][j];
			for (int k = 0; k < number_of_terms; k++)
			{
				energy_prediction += 
				set_of_force_constants[j * number_of_terms + k] * 
				cos(multiplicities[k] * angle / 180.0 * M_PI);
			}
		}
		string str = to_string(energy[i]);
		for (int l = 0; l < MAX_LENGTH - (int) str.length(); l++)
		{
			outstream << ' ';
		}
		outstream << str;

		str = to_string(energy_prediction);
		for (int l = 0; l < MAX_LENGTH - (int) str.length(); l++)
		{
			outstream << ' ';
		}
		outstream << str << '\n';
	}

	outstream << "----------------------------------------------------\n";
}

void conclude(
	ostream& outstream, 
	double sum_error,
	int trial,
	double M, 
	chrono::high_resolution_clock::time_point t1, 
	chrono::high_resolution_clock::time_point t2)
{
	outstream << "Error average: " << sum_error/trial << '\n';
	outstream << "Standard deviation: " << M/(trial-1) << '\n';
    outstream << "95% confidence interval:" << 1.960 * M/(trial-1)/sqrt(trial) << '\n';
	outstream << "Total runtime: " << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0 << '\n';
}

int main(int argc, char* argv[])
{
	if (not strcmp(argv[1], "-h") || not strcmp(argv[1], "--help"))
	{
		cerr << "Usage: <BIN_FILE_PATH>\n"
			 << "[-h|--help]\n"
			 << "[-t|--test <DATA_PATH> <FORCE_CONSTANTS_PATH>]\n"
			 << "[-r|--run <DATA_PATH>]\n"
			 << "[-r|--run <DATA_PATH> <LOG_FILE_PATH>]\n";
	}
	else if (not strcmp(argv[1], "-t") || not strcmp(argv[1], "--test"))
	{
		if (argc < 4)
		{
			cerr << "Bad arguments: At least 4 arguments expected. Type -h or --help for complete guide.\n";
			return 1;
		}
		input(argv[2]);

		// initialize input stream
		ifstream force_constants_file;
		force_constants_file.open(argv[3]);
		
		vector<double> set_of_force_constants_reduced(number_of_distinct_angles * 6);
		for (int i = 0; i < number_of_distinct_angles; i++)
		{
			// angle representation
			int idx;
			vector<string> s(4);
			force_constants_file >> s[0] >> s[1] >> s[2] >> s[3];
			string 
				t1 = s[0] + ' ' + s[1] + ' ' + s[2] + ' ' + s[3],
				t2 = s[3] + ' ' + s[2] + ' ' + s[1] + ' ' + s[0];
			
			// get id of the angle
			if (angles_id_dict.find(t1) != angles_id_dict.end())
			{
				idx = angles_id_dict[t1];
			}
			else if (angles_id_dict.find(t2) != angles_id_dict.end())
			{
				idx = angles_id_dict[t2];
			}
			else
			{
				cerr << "Unknown angle types. Please check the force constants file and try again.\n";
				return 1;
			}

			for (int j = 0; j < 6; j++)
			{
				force_constants_file >> set_of_force_constants_reduced[idx * 6 + j];
			}
		}

		// expand the set of force constants
		set_of_force_constants_without_offset_constant.resize(number_of_angles * 6);
		for (int j = 0; j < number_of_angles; j++)
		{
			for (int k = 0; k < 6; k++)
			{
				set_of_force_constants_without_offset_constant[6 * j + k] = 
				set_of_force_constants_reduced[6 * angles_id[j] + k];
			}
		}
		
		vector<double> set_of_force_constants(number_of_angles * 6 + 1);
		for (int i = 0; i < (int) set_of_force_constants_without_offset_constant.size(); i++)
		{
			set_of_force_constants[i] = set_of_force_constants_without_offset_constant[i];
		}
		// find best offset constant
		set_of_force_constants[set_of_force_constants.size() - 1]
		 = simulated_annealing(RMSE_with_respect_to_offset_constant, 
		 					   1, 
							   1.0, 
							   10000)[0];

		cerr << "RMSE: "
		     << rmse(set_of_force_constants, 
					 energy, 
					 angles, 
					 multiplicities, 
					 number_of_data_points, 
					 number_of_angles, 
					 6);
	}
	else if (not strcmp(argv[1], "-r") || not strcmp(argv[1], "--run"))
	{	
		if (argc <= 2)
		{
			cerr << "Bad arguments: At least 3 arguments expected. Type -h or --help for complete guide.\n";
			return 1;
		}
		if (argc > 3)
		{
			log_file.open(argv[3]);
		}
		cerr << "START\n";
		input(argv[2]);
		correlation_preparation();
		int trial = 2;
		double sum_error = 0;
		chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
		double prev_mean = 0, current_mean = 0, M = 0;
		for (int i = 1; i <= trial; i++)
		{
			vector<double> set_of_force_constants_reduced = optimization_scheme();
			if (argc > 3)
			{
				output(log_file, 
					   i, 
					   sum_error, 
					   prev_mean, 
					   current_mean, 
					   M, 
					   set_of_force_constants_reduced);
			}
			else
			{
				output(cout, 
					   i, 
					   sum_error, 
					   prev_mean, 
					   current_mean, 
					   M, 
					   set_of_force_constants_reduced);
			}
		}
		chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
		if (argc > 3)
		{
			conclude(log_file, 
					 sum_error, 
					 trial, 
					 M, 
					 t1, 
					 t2);
			log_file.close();
		}
		else
		{
			conclude(cout, 
					 sum_error, 
					 trial, 
					 M, 
					 t1, 
					 t2);
		}
		cerr << "FINISH\n";
	}
	else
	{
		cerr << "Unknown command. Type -h or --help for complete guide.\n";
	}
}