#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "src/ExtremaFinding/GlobalMinimumFinder.h"

int multiplicities[6] = {0, 1, 2, 3, 4, 6};
string quantum_mechanic_data = "soft_dihe_C1_C2.dat";
int
	number_of_terms = 6,
	number_of_data_points = 36, 
	number_of_angles = 9, 
	number_of_neighbors_of_central_atom_1 = 4,
	number_of_neighbors_of_central_atom_2 = 4;
vector<vector<double> > angles, interaction;
vector<double> energy, set_of_force_constants;
vector<int> rank_by_pearson;

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

void preprocess()
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
	rank_by_pearson.resize(number_of_angles * number_of_terms);

	// read input
	ifstream quantum_mechanics_file;
	string dummy;
	quantum_mechanics_file.open(quantum_mechanic_data);
	int cnt = 0;
	while(cnt < 3)
	{
		getline(quantum_mechanics_file, dummy);
		if (dummy[0] == '$')
		{
			cnt++;
		}
	}
	for (int i = 0; i < number_of_data_points; i++)
	{
		for (int j = 0; j < number_of_angles; j++)
		{
			quantum_mechanics_file >> angles[i][j] >> interaction[i][j];
		}
		quantum_mechanics_file >> energy[i];
	}
	quantum_mechanics_file.close();

	vector<tuple<double, int> > temp(number_of_angles * number_of_terms);
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
			sumB = 0, sumB2 = 0, sumAB = 0;
			for (int k = 0; k < number_of_data_points; k++)
			{
				double b = cos(multiplicities[j] * angles[k][i] * M_PI / 180.0);
				sumB += b;
				sumB2 += b * b;
				sumAB += a[k] * b;
			}

			temp[number_of_terms * i + j] = 
			{
				(number_of_data_points * sumAB - sumA * sumB) / 
				sqrt(
					(number_of_data_points * sumA2 - sumA * sumA) * 
					(number_of_data_points * sumB2 - sumB * sumB)
				),
				number_of_terms * i + j
			};
		}
	}
	sort(temp.begin(), temp.end());
	for (int i = 0; i < (int) temp.size(); i++)
	{
		rank_by_pearson[get<1>(temp[i])] = i;
	}
}

double rmse(
	vector<double>& set_of_force_constants)
{
	double sum_of_square_error = 0;
	for (int i = 0; i < number_of_data_points; i++)
	{
		double error = energy[i];
		for (int j = 0; j < number_of_angles; j++)
		{
			double angle = angles[i][j];
			for (int k = 0; k < number_of_terms; k++)
			{
				error -= 
				set_of_force_constants[j * number_of_terms + k] * 
				(cos(multiplicities[k] * angle / 180.0 * M_PI));
			}
		}
		sum_of_square_error += error * error;
	}
	return sqrt(sum_of_square_error/number_of_data_points);
}

double objective_function(
	vector<double>& set_of_force_constants)
{
	double result = rmse(set_of_force_constants);
	
	// rank force constants
	vector<tuple<double, int> > temp(set_of_force_constants.size());
	for (int i = 0; i < (int) temp.size(); i++)
	{
		temp[i] = {set_of_force_constants[i], i};
	}
	sort(temp.begin(), temp.end());
	// calculate Spearman's rank correlation coefficient
	double ans = 0;
	for (int i = 0; i < (int) temp.size(); i++)
	{
		ans += 6 * (i - rank_by_pearson[i]) * (i - rank_by_pearson[i]);
	}
	result += -(1 - ans/(temp.size() * (temp.size() * temp.size() - 1)));

	return result;
}

int main()
{
	preprocess();
	set_of_force_constants = simulated_annealing(objective_function, number_of_angles * number_of_terms, 3.0);
	cerr << "Error: " << rmse(set_of_force_constants);
	for (int i = 0; i < number_of_angles; i++)
	{
		cerr << "Force constants of angle " << i << ": ";
		for (int j = 0; j < number_of_terms; j++)
		{
			cerr << set_of_force_constants[i * number_of_terms + j] << " ";
		}
		cerr << '\n';
	}
}