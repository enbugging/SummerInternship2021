#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "src/ExtremaFinding/GlobalMinimumFinder.h"

const int number_of_angles = 9;
//string quantum_mechanic_data = "soft_dihe_C1_C2.dat";
string quantum_mechanic_data = "soft_dihe_C2_C3.dat";
//string quantum_mechanic_data = "propane_soft_dihe_C1_C2.dat";
int
	number_of_terms = 3,
	number_of_data_points = 36, 
	number_of_distinct_angles, 
	nbrs_of_central_atom_1,
	nbrs_of_central_atom_2,
	angles_id[number_of_angles], 
	multiplicities[5] = {1, 2, 3, 4, 6};
vector<vector<double> > 
	angles, 
	interaction;
vector<double> 
	energy,  
	pearson_rank, 
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
int principle_multiplicity(
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

double compare(
	double a, 
	double b)
{
	// discontinuous version
	return (a > b);
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
		string s;
		getline(quantum_mechanics_file, s);
		while(s[0] != '|') s.erase(0, 1);
		s.erase(0, 1);
		while(s[0] == ' ') s.erase(0, 1);

		
		// get id of the angle
		if (angles_id_dict.find(s) != angles_id_dict.end())
		{
			angles_id[i] = angles_id_dict[s];
		}
		else
		{
			string t = s;
			reverse(t.begin(), t.end());
			angles_id_dict[s] = number_of_distinct_angles;
			angles_id_dict[t] = number_of_distinct_angles++;
			angles_id[i] = angles_id_dict[s];
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
	pearson_rank.resize(number_of_angles * (number_of_terms - 1));
	interact_mag.resize(number_of_angles * (number_of_terms - 1));
	vector<double> temp(number_of_angles * (number_of_terms - 1));
	int cnt = 0;
	for (int i = 0; i < number_of_angles; i++)
	{
		double sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;
		// preparing 2 sequences of interaction and cos(multiplicity * angle)
		vector<double> a(number_of_data_points);
		for (int j = 0; j < number_of_data_points; j++)
		{
			a[j] = interaction[j][i];
			sumA += a[j];
			sumA2 += a[j] * a[j];
		}
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] != principle_multiplicity(nbrs_of_central_atom_1, nbrs_of_central_atom_2))
			{
				sumB = 0, sumB2 = 0, sumAB = 0;
				for (int k = 0; k < number_of_data_points; k++)
				{
					double b = cos(multiplicities[j] * angles[k][i] * M_PI / 180.0);
					sumB += b;
					sumB2 += b * b;
					sumAB += a[k] * b;
				}

				temp[cnt++] = 
				abs(number_of_data_points * sumAB - sumA * sumB) / 
				sqrt(
					(number_of_data_points * sumA2 - sumA * sumA) * 
					(number_of_data_points * sumB2 - sumB * sumB)
				);
			}
		}
	}
	for (int i = 0; i < cnt; i++)
	{
		pearson_rank[i] = 1;
		for (int j = 0; j < cnt; j++)
		{
			pearson_rank[i] += compare(temp[i], temp[j]);
		}
	}

	for (int i = 0; i < number_of_angles; i++)
	{
		double max_a = numeric_limits<double>::min(), min_a = numeric_limits<double>::max(), delta;
		for (int j = 0; j < number_of_data_points; j++)
		{
			min_a = min(min_a, interaction[j][i]);
			max_a = max(max_a, interaction[j][i]);
		}
		delta = max_a - min_a;
		cnt = 0;
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] != principle_multiplicity(nbrs_of_central_atom_1, nbrs_of_central_atom_2))
			{
				interact_mag[cnt++] = delta;
			}
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
	// rank force constants
	int cnt = 0;
	vector<double> 
		temp(number_of_angles * (number_of_terms - 1)), 
		rank(number_of_angles * (number_of_terms - 1));
	for (int i = 0; i < number_of_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] != principle_multiplicity(nbrs_of_central_atom_1, nbrs_of_central_atom_2))
			{
				temp[cnt++] = set_of_force_constants[number_of_terms * i + j];
			}
		}
	}
	for (int i = 0; i < cnt; i++)
	{
		rank[i] = 1;
		for (int j = 0; j < cnt; j++)
		{
			rank[i] += compare(temp[i], temp[j]);
		}
	}
	// calculate Spearman's rank correlation coefficient
	double sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0, spearman_coeff, pearson_coeff;
	for (int i = 0; i < cnt; i++)
	{
		sumA += pearson_rank[i];
		sumB += rank[i];
		sumAB += pearson_rank[i] * rank[i];
		sumA2 += pearson_rank[i] * pearson_rank[i];
		sumB2 += rank[i] * rank[i];
	}
	if ((cnt * sumA2 - sumA * sumA) <= 0 && (cnt * sumB2 - sumB * sumB) <= 0)
	{
		spearman_coeff = 0;
	}
	else if ((cnt * sumA2 - sumA * sumA) <= 0 || (cnt * sumB2 - sumB * sumB) <= 0)
	{
		spearman_coeff = 1;
	}
	else{
		spearman_coeff = 
			(1 - 
			abs(cnt * sumAB - sumA * sumB) / 
			sqrt(
				(cnt * sumA2 - sumA * sumA) * 
				(cnt * sumB2 - sumB * sumB)
			));
	}
	// calculate Pearson's correlation coefficient
	sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;
	for (int i = 0; i < cnt; i++)
	{
		sumA += interact_mag[i];
		sumB += temp[i];
		sumAB += interact_mag[i] * temp[i];
		sumA2 += interact_mag[i] * interact_mag[i];
		sumB2 += temp[i] * temp[i];
	}
	if ((cnt * sumA2 - sumA * sumA) <= 0 && (cnt * sumB2 - sumB * sumB) <= 0)
	{
		pearson_coeff = 0;
	}
	else if ((cnt * sumA2 - sumA * sumA) <= 0 || (cnt * sumB2 - sumB * sumB) <= 0)
	{
		pearson_coeff = 1;
	}
	else{
		pearson_coeff = 
			(1 - 
			abs(cnt * sumAB - sumA * sumB) / 
			sqrt(
				(cnt * sumA2 - sumA * sumA) * 
				(cnt * sumB2 - sumB * sumB)
			));
	}
	return sqrt(spearman_coeff * pearson_coeff);
}

double objective_function(
	vector<double>& set_of_force_constants_reduced)
{
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
	return rmse(set_of_force_constants);// + interaction_correlation(set_of_force_constants);
}

double test_function(
	vector<double>& set_of_force_constants_reduced)
{
	// expand the set of force constants
	double result = 0;
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
	double sum_of_square_error = 0;
	for (int i = 0; i < number_of_data_points; i++)
	{
		double error = energy[i] - set_of_force_constants[set_of_force_constants.size()-1];
		log_file << energy[i] << ' ' << error << '\n';
		for (int j = 0; j < number_of_angles; j++)
		{
			double angle = angles[i][j];
			for (int k = 0; k < number_of_terms; k++)
			{
				error -= 
				set_of_force_constants[j * number_of_terms + k] * 
				cos(multiplicities[k] * angle / 180.0 * M_PI);
				
				//*
				log_file << "Minus contribution of mult " << multiplicities[k] << ", angle " << j << ": " << 
				set_of_force_constants[j * number_of_terms + k] << " x " << 
				cos(multiplicities[k] * angle / 180.0 * M_PI) << " = " << 
				set_of_force_constants[j * number_of_terms + k] * 
				cos(multiplicities[k] * angle / 180.0 * M_PI) << '\n';
				//*/
			}
		}
		sum_of_square_error += error * error;
		log_file << "ADD: " << error << '\n';
	}
	result += sqrt(sum_of_square_error/number_of_data_points);
	/*
	int cnt = 0;
	vector<double> 
		temp(number_of_angles * (number_of_terms - 1)), 
		rank(number_of_angles * (number_of_terms - 1));
	for (int i = 0; i < number_of_angles; i++)
	{
		for (int j = 0; j < number_of_terms; j++)
		{
			if (multiplicities[j] != principle_multiplicity(nbrs_of_central_atom_1, nbrs_of_central_atom_2))
			{
				temp[cnt++] = set_of_force_constants[number_of_terms * i + j];
			}
		}
	}
	for (int i = 0; i < cnt; i++)
	{
		rank[i] = 1;
		for (int j = 0; j < cnt; j++)
		{
			rank[i] += compare(temp[i], temp[j]);
		}
	}
	// calculate Spearman's rank correlation coefficient
	double sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0, spearman_coeff, pearson_coeff;
	for (int i = 0; i < cnt; i++)
	{
		sumA += pearson_rank[i];
		sumB += rank[i];
		sumAB += pearson_rank[i] * rank[i];
		sumA2 += pearson_rank[i] * pearson_rank[i];
		sumB2 += rank[i] * rank[i];
	}
	if ((cnt * sumA2 - sumA * sumA) <= 0 && (cnt * sumB2 - sumB * sumB) <= 0)
	{
		spearman_coeff = 0;
	}
	else if ((cnt * sumA2 - sumA * sumA) <= 0 || (cnt * sumB2 - sumB * sumB) <= 0)
	{
		spearman_coeff = 1;
	}
	else{
		spearman_coeff = 
			(1 - 
			abs(cnt * sumAB - sumA * sumB) / 
			sqrt(
				(cnt * sumA2 - sumA * sumA) * 
				(cnt * sumB2 - sumB * sumB)
			));
	}
	// calculate Pearson's correlation coefficient
	sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;
	for (int i = 0; i < cnt; i++)
	{
		sumA += interact_mag[i];
		sumB += temp[i];
		sumAB += interact_mag[i] * temp[i];
		sumA2 += interact_mag[i] * interact_mag[i];
		sumB2 += temp[i] * temp[i];
	}
	if ((cnt * sumA2 - sumA * sumA) <= 0 && (cnt * sumB2 - sumB * sumB) <= 0)
	{
		pearson_coeff = 0;
	}
	else if ((cnt * sumA2 - sumA * sumA) <= 0 || (cnt * sumB2 - sumB * sumB) <= 0)
	{
		pearson_coeff = 1;
	}
	else{
		pearson_coeff = 
			(1 - 
			abs(cnt * sumAB - sumA * sumB) / 
			sqrt(
				(cnt * sumA2 - sumA * sumA) * 
				(cnt * sumB2 - sumB * sumB)
			));
	}
	log_file << "Correlation: " << sqrt(spearman_coeff * pearson_coeff) << '\n';
	result += sqrt(spearman_coeff * pearson_coeff);
	//*/
	return result;
}


int main()
{
	log_file.open("log.txt");

	input();
	//correlation_preparation();
	int trial = 2;
	double sum_error = 0;
	auto t1 = chrono::high_resolution_clock::now();
	double prev_mean = 0, current_mean = 0, M = 0;
	for (int i = 1; i <= trial; i++)
	{
		vector<double> set_of_force_constants_reduced = simulated_annealing(objective_function, number_of_distinct_angles * number_of_terms + 1, 1.0);
		double error = test_function(set_of_force_constants_reduced);
		sum_error += error;
		prev_mean = current_mean;
		current_mean = (prev_mean * (i - 1) + error)/i;
		M += (error - current_mean) * (error - prev_mean);

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
		
		log_file << "Error: " << error << '\n';
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
	log_file << "Error average: " << sum_error/trial << '\n';
	log_file << "Standard deviation: " << M/(trial-1) << '\n';
    log_file << "95% confidence interval:" << 1.960 * M/(trial-1)/sqrt(trial) << '\n';
	auto t2 = chrono::high_resolution_clock::now();
	log_file << "Total runtime: " << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0 << '\n';
	//log_file << "Average runtime: " << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()/1000.0/trial << '\n'; 
	
	log_file.close();
}