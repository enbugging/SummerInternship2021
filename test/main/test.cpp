#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "../../src/ExtremaFinding/GlobalMinimumFinder.h"

const int number_of_angles = 9;
//string quantum_mechanic_data = "ethane_dihe_c1_c2.dat";
//string quantum_mechanic_data = "propane_dihe_c1_c2.dat";
string quantum_mechanic_data = "../../butane_dihe_c2_c3.dat";
int
	number_of_terms = 3,
	number_of_data_points = 36, 
	number_of_distinct_angles, 
	nbrs_of_central_atom_1,
	nbrs_of_central_atom_2, 
	angles_id[number_of_angles], 
	multiplicities[5] = {1, 2, 3, 4, 6};
double 
    simplicity_accuracy_trading[5] = {0.01, 0.015, 0.02, 0.025, 0.035};
vector<vector<double> > 
	angles, 
	interaction, 
	pearson_interact_vs_angles;
vector<double> 
	energy,
	interact_mag;
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
void input()
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
			min_a = numeric_limits<double>::max(), 
			delta;
		for (int j = 0; j < number_of_data_points; j++)
		{
			a[j] = interaction[j][i];
			sumA += a[j];
			sumA2 += a[j] * a[j];
			min_a = min(min_a, a[j]);
			max_a = max(max_a, a[j]);
		}
		delta = max_a - min_a;
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
			
			pearson_interact_vs_angles[i][j] = delta * p;
		}
	}
}

/* 
--------------------------------------------------------------------------------------
Objective functions
 */
double rmse(
	vector<double>& set_of_force_constants)
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

double interaction_correlation(
	vector<double>& set_of_force_constants)
{
	double 
		sumA = 0, 
		sumB = 0, 
		sumAB = 0, 
		sumA2 = 0, 
		sumB2 = 0;
	int cnt = number_of_angles * (number_of_terms - 1) + 1;
	for (int i = 0; i < number_of_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] == main_multiplicity(nbrs_of_central_atom_1, nbrs_of_central_atom_2)) continue;
			double s = set_of_force_constants[number_of_terms * i + j];
			s *= abs(s);
			sumA += s;
			sumA2 += s * s;
			sumB += pearson_interact_vs_angles[i][j];
			sumB2 += pearson_interact_vs_angles[i][j] * pearson_interact_vs_angles[i][j];
			sumAB += s * pearson_interact_vs_angles[i][j];
		}
	}
	if ((cnt * sumA2 - sumA * sumA) == 0 && (cnt * sumB2 - sumB * sumB) == 0) return 1;
	if ((cnt * sumA2 - sumA * sumA) == 0 || (cnt * sumB2 - sumB * sumB) == 0) return 0;
	return (1 -
				(cnt * sumAB - sumA * sumB) / 
				sqrt(
					(cnt * sumA2 - sumA * sumA) * 
					(cnt * sumB2 - sumB * sumB)
				)
			)/2;
}

double variance_of_main_multiplicity(
	vector<double>& set_of_force_constants)
{
	double sum = 0;
	int idx_main_mult = -1;
	for (int i = 0; i < number_of_terms; i++)
	{
		if (multiplicities[i] == main_multiplicity(nbrs_of_central_atom_1, nbrs_of_central_atom_2))
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
					set_of_force_constants[i * number_of_terms + idx_main_mult] - 
					set_of_force_constants[j * number_of_terms + idx_main_mult];
				sum += d * d;
			}
		}
	}
	return sqrt(sum/(number_of_angles * (number_of_angles - 1)/2));
}

double magnitude_of_force_constants(
	vector<double>& set_of_force_constants_reduced)
{
	double sum = 0;
	int no_of_dist_angles = 
		(set_of_force_constants_reduced.size() - 1)/number_of_terms;
	for (int i = 0; i < no_of_dist_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] == main_multiplicity(nbrs_of_central_atom_1, nbrs_of_central_atom_2)) continue;
			sum += 
				set_of_force_constants_reduced[i * number_of_terms + j] * 
				set_of_force_constants_reduced[i * number_of_terms + j] / 
				multiplicities[j];
		}
	}
	return sqrt(sum/(no_of_dist_angles * (number_of_terms - 1)));
}
/* 
--------------------------------------------------------------------------------------
Simplicity-accuracy trading functions
 */
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
	double cutoff = 0.0)
{
    double error = 0;
	for (int i = 0; i < number_of_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			if (
				multiplicities[j] == 
				main_multiplicity(
					nbrs_of_central_atom_1, 
					nbrs_of_central_atom_2)
				)
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
	//*
	double r = rmse(set_of_force_constants);
	result += r;
	result += simplicity_vs_accuracy(set_of_force_constants);
	result += exp(-r) * variance_of_main_multiplicity(set_of_force_constants);
	//result += magnitude_of_force_constants(set_of_force_constants_reduced);
	result += interaction_correlation(set_of_force_constants);
	return result;
}
/* 
--------------------------------------------------------------------------------------
Brute-force method
 */
int dummy_mask, N;

double brute_force_objective_function(
	vector<double>& reduced_input)
{
	// construct reduced set of force constants
	vector<double> set_of_force_constants_reduced(N + 1);
	int cnt = 0;
	for (int i = 0; i < N; i++)
	{
		if ((dummy_mask>>i)&1)
		{
			set_of_force_constants_reduced[i] = reduced_input[cnt++];
		}
	}
	set_of_force_constants_reduced[set_of_force_constants_reduced.size() - 1] 
	= reduced_input[cnt++];
	
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
	//*
	double r = rmse(set_of_force_constants);
	result += r;
	result += exp(-r) * variance_of_main_multiplicity(set_of_force_constants);
	//result += magnitude_of_force_constants(set_of_force_constants_reduced);
	result += interaction_correlation(set_of_force_constants);
	return result;
}

vector<double> brute_force()
{
	double minimal_r = numeric_limits<double>::max();
	vector<double> answer;
	N = number_of_distinct_angles * number_of_terms;
	for (dummy_mask = 0; dummy_mask < (1<<N); dummy_mask++)
	{
		int cnt = 0;
		for (int i = 0; i < N; i++)
		{
			if ((dummy_mask>>i)&1) cnt++;
		}

		vector<double> reduced_input = simulated_annealing(brute_force_objective_function, cnt + 1, 1.0, 1000);
		double r = brute_force_objective_function(reduced_input);
		for (int i = 0; i < number_of_angles; i++)
		{
			for (int j = 0; j < number_of_terms; j++)
			{
				int d = number_of_terms * angles_id[i] + j;
				if ((dummy_mask>>d)&1)
				{
					r += simplicity_accuracy_trading[j];
				}
			}
		}

		// construct reduced set of force constants
		vector<double> set_of_force_constants_reduced(N + 1);
		cnt = 0;
		for (int i = 0; i < N; i++)
		{
			if ((dummy_mask>>i)&1)
			{
				set_of_force_constants_reduced[i] = reduced_input[cnt++];
			}
		}
		set_of_force_constants_reduced[set_of_force_constants_reduced.size() - 1] 
		= reduced_input[cnt++];
		
		// expand the set of force constants
		vector<double> set_of_force_constants(number_of_angles * number_of_terms + 1);
		set_of_force_constants[set_of_force_constants.size() - 1] = 
		set_of_force_constants_reduced[set_of_force_constants_reduced.size() - 1];
		for (int i = 0; i < number_of_angles; i++)
		{
			for (int j = 0; j < number_of_terms; j++)
			{
				set_of_force_constants[number_of_terms * i + j] = 
				set_of_force_constants_reduced[number_of_terms * angles_id[i] + j];
			}
		}
		
		string s = "";
		for (int i = 0; i < N; i++)
		{
			if ((dummy_mask>>i)&1) s += "1";
			else s += "0";
		}

		if (r < minimal_r)
		{	
			minimal_r = r;
			answer = set_of_force_constants;
		}
	}
	return answer;
}


int main()
{
	log_file.open("log.txt");
	cerr << "START\n";
	input();
	correlation_preparation();
	int trial = 2;
	double sum_error = 0;
	auto t1 = chrono::high_resolution_clock::now();
	double prev_mean = 0, current_mean = 0, M = 0;
	//*/
	for (int i = 1; i <= trial; i++)
	{
		vector<double> set_of_force_constants_reduced = simulated_annealing(objective_function, number_of_distinct_angles * number_of_terms + 1, 1.0);
		// expand the set of force constants
		vector<double> set_of_force_constants(number_of_angles * number_of_terms + 1);
		set_of_force_constants[set_of_force_constants.size() - 1] = 
		set_of_force_constants_reduced[set_of_force_constants_reduced.size() - 1];
		for (int i = 0; i < number_of_angles; i++)
		{
			for (int j = 0; j < number_of_terms; j++)
			{
				set_of_force_constants[number_of_terms * i + j] = 
				set_of_force_constants_reduced[number_of_terms * angles_id[i] + j];
			}
		}
		
		double error = rmse(set_of_force_constants);
		sum_error += error;
		prev_mean = current_mean;
		current_mean = (prev_mean * (i - 1) + error)/i;
		M += (error - current_mean) * (error - prev_mean);
		
		log_file << "Error: " 
				<< error << ' ' 
				<< simplicity_vs_accuracy(set_of_force_constants) << ' ' 
				<< interaction_correlation(set_of_force_constants) << ' ' 
				<< exp(-error) * variance_of_main_multiplicity(set_of_force_constants) << '\n';
		log_file << "Offset constant: " << set_of_force_constants[set_of_force_constants.size() - 1] << '\n';
		for (int k = 0; k < number_of_angles; k++)
		{
			log_file << "Force constants of angle " << k << ": ";
			for (int j = 0; j < number_of_terms; j++)
			{
				log_file << set_of_force_constants[k * number_of_terms + j] << " ";
			}
			log_file << '\n';
		}
	}
	//*/
	/*
	for (int i = 1; i <= trial; i++)
	{
		vector<double> set_of_force_constants = brute_force();
		double error = rmse(set_of_force_constants);
		sum_error += error;
		prev_mean = current_mean;
		current_mean = (prev_mean * (i - 1) + error)/i;
		M += (error - current_mean) * (error - prev_mean);
		
		log_file << "Error: " 
				<< error << ' ' 
				<< simplicity_vs_accuracy(set_of_force_constants) << ' ' 
				<< interaction_correlation(set_of_force_constants) << ' ' 
				<< exp(-error) * variance_of_main_multiplicity(set_of_force_constants) << '\n';
		log_file << "Offset constant: " << set_of_force_constants[set_of_force_constants.size() - 1] << '\n';
		
		for (int k = 0; k < number_of_angles; k++)
		{
			log_file << "Force constants of angle " << k << ": ";
			for (int j = 0; j < number_of_terms; j++)
			{
				log_file << set_of_force_constants[k * number_of_terms + j] << " ";
			}
			log_file << '\n';
		}
	}
	//*/
	log_file << "Error average: " << sum_error/trial << '\n';
	log_file << "Standard deviation: " << M/(trial-1) << '\n';
    log_file << "95% confidence interval:" << 1.960 * M/(trial-1)/sqrt(trial) << '\n';
	auto t2 = chrono::high_resolution_clock::now();
	log_file << "Total runtime: " << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0 << '\n';
	//log_file << "Average runtime: " << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0/trial << '\n'; 
	
	log_file.close();
}