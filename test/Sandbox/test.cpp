#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "../../src/ExtremaFinding/GlobalMinimumFinder.h"
#include "../../src/ObjectiveFunctions/ObjectiveFunctions.h"
#include "../../src/SimplicityAccuracy/SimplicityAccuracy.h"

const int number_of_angles = 9;
//string quantum_mechanic_data = "ethane_dihe_c1_c2.dat";
//string quantum_mechanic_data = "propane_dihe_c1_c2.dat";
//string quantum_mechanic_data = "butane_dihe_c2_c3.dat";
int
	number_of_terms = 3,
	number_of_data_points = 36, 
	number_of_distinct_angles, 
	nbrs_of_central_atom_1,
	nbrs_of_central_atom_2, 
	main_mult, 
	angles_id[number_of_angles], 
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
	magnitude_of_interaction;
map<string, int> angles_id_dict;

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
--------------------------------------------------------------------------------------
Preprocess functions
 */
void input(string quantum_mechanic_data)
{
	// initialization
	angles.resize(number_of_data_points);
	interaction.resize(number_of_data_points);
	energy.resize(number_of_data_points);
	for (int i = 0; i < number_of_data_points; i++)
	{
		angles[i].resize(number_of_angles);
		interaction[i].resize(number_of_angles);
	}

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
	getline(quantum_mechanics_file, dummy);
	
	// read angle descrition
	getline(quantum_mechanics_file, dummy);
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

int main(int argc, char* argv[])
{
	log_file.open("log.txt");
	cerr << "START\n";
	input(argv[1]);
	correlation_preparation();
	int trial = 2;
	double sum_error = 0;
	auto t1 = chrono::high_resolution_clock::now();
	double prev_mean = 0, current_mean = 0, M = 0;
	for (int i = 1; i <= trial; i++)
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
		w_total = sqrt(sum_of_energy_square)/(epsilon_main);
		
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
		
		//cerr << "4\n";
		// expand the set of force constants
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
		log_file << "Error: " 
				<< error << '\n';
		log_file << "Offset constant: " << set_of_force_constants[set_of_force_constants.size() - 1] << '\n';
		for (int j = 0; j < number_of_angles; j++)
		{
			log_file << "Force constants of angle " << j << ": ";
			for (int k = 0; k < number_of_terms; k++)
			{
				log_file << set_of_force_constants[j * number_of_terms + k] << " ";
			}
			log_file << '\n';
		}
	}
	log_file << "Error average: " << sum_error/trial << '\n';
	log_file << "Standard deviation: " << M/(trial-1) << '\n';
    log_file << "95% confidence interval:" << 1.960 * M/(trial-1)/sqrt(trial) << '\n';
	auto t2 = chrono::high_resolution_clock::now();
	log_file << "Total runtime: " << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0 << '\n';
	//log_file << "Average runtime: " << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0/trial << '\n'; 
	
	log_file.close();
}