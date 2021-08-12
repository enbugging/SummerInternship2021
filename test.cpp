#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "src/ExtremaFinding/GlobalMinimumFinder.h"
#include "src/ExtremaFinding/ForceConstantFinder.h"

#define random_step(rng,min,max) min + (double)rng() / UINT_MAX * (max - min)
normal_distribution<double> gaussian (0.0,1.0);

double 
    cutoff,
    upper_limit = 3.0;
int
	number_of_terms = 4,
	angle_step = 10,
    number_of_tests = 1,
    result = 0;
vector<double>
	angles,
	test_points,
	test_weights,
	test_force_constants,
    force_constants;

/* 
Auxiliary functions
 */
int 
    multiplicities[6] = {0, 1, 2, 3, 4, 6};
double 
    dummy_cutoff, 
    simplicity_accuracy_trading[6] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.35};

double force_field_calculate(
    vector<double>& force_constants,
    double angle)
{
    double E = force_constants[0];
    for (int i = 1; i < (int)force_constants.size(); i++)
    {
        E += force_constants[i] * (cos(multiplicities[i] * angle / 180.0 * M_PI));
    }
    return E;
}

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
    return non_analytic_smooth_function(x) / (non_analytic_smooth_function(x) + non_analytic_smooth_function(1 - x));
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
    //else return cubical((abs(t) - cutoff)/border_width);
    //else return (abs(t) - cutoff)/border_width;
    else return smooth_transition_function((abs(t) - cutoff) / border_width);
    //return 0.9999;
}

double rmse(
    vector<double>& force_constants)
{
    if (angles.empty()) return 0;
    // if the function is unweighted
    if (test_weights.empty())
    {
        test_weights.assign(angles.size(), 1.0);
    }

    double sum_of_squares_of_error = 0, sum_of_weights = 0;
    for (int i = 0; i < (int) angles.size(); i++)
    {
        double error = force_field_calculate(force_constants, angles[i]) - test_points[i];
        sum_of_squares_of_error += test_weights[i] * error * error;
        sum_of_weights += test_weights[i];
    }
    return sqrt(sum_of_squares_of_error/sum_of_weights);
}

double rmse_with_cutoff_and_simplicity_accuracy_trading(
    vector<double>& force_constants)
{
    if (angles.empty()) return 0;
    // if the function is unweighted
    if (test_weights.empty())
    {
        test_weights.assign(angles.size(), 1.0);
    }

    double sum_of_squares_of_error = 0, sum_of_weights = 0;
    for (int i = 0; i < (int) angles.size(); i++)
    {
        double error = force_field_calculate(force_constants, angles[i]) - test_points[i];
        sum_of_squares_of_error += test_weights[i] * error * error;
        sum_of_weights += test_weights[i];
    }
    sum_of_squares_of_error = sqrt(sum_of_squares_of_error/sum_of_weights);
    for (int i = 0; i < (int) force_constants.size(); i++)
    {
        sum_of_squares_of_error += simplicity_accuracy_trading[i] * coefficient(force_constants[i], cutoff);
    }
    return sum_of_squares_of_error;
}

void preprocess1()
{
    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count(); // 17 for testing
    mt19937 rng(seed);

    cutoff = 0.0;
    test_force_constants.resize(number_of_terms);
    for (int i = 0; i < number_of_terms; i++)
    {
        if (rng() & 1)
        {
            if (rng() & 1)
            {
                test_force_constants[i] = random_step(rng, cutoff, upper_limit);
            }
            else
            {
                test_force_constants[i] = random_step(rng, -upper_limit, -cutoff);
            }
        }
        else
        {
            test_force_constants[i] = random_step(rng, -cutoff, cutoff);
        }
    }
    angles.clear();
    for (int i = -180 + angle_step; i <= 180; i += angle_step)
    {
        angles.push_back(i);
    }
    test_points.resize(angles.size());
    for (int i = 0; i < (int) angles.size(); i++)
    {
        test_points[i] = force_field_calculate(test_force_constants, angles[i]);// + 2.0e-6 * gaussian(rng);
    }
}

int summary1()
{
    /*
    printf("-----------------------------------------------------\n");
    printf("Origi constant: ");
    for (int i = 0; i < number_of_terms; i++)
    {
        printf("%.6lf ", test_force_constants[i]);
    }
    printf("\n");
    printf("Force constant: ");
    for (int i = 0; i < number_of_terms; i++)
    {
        printf("%.6lf ", force_constants[i]);
    }
    printf("\n");
    */
    bool match = true;
    for (int i = 0; i < number_of_terms; i++)
	{
        //cerr << "FINAL " << i << " : " << test_force_constants[i] << ' ' << force_constants[i] << '\n';
        match &= (abs(force_constants[i] - (abs(test_force_constants[i]) > cutoff ? test_force_constants[i] : 0.0)) <= 1.0e-6);
	}
    if (not match)
    {
        for (int i = 0; i < number_of_terms; i++)
        {
            cerr << "WRONG " << i << " : " << test_force_constants[i] << ' ' << force_constants[i] << '\n';
        }
        cerr << '\n';
    }
    return match;
}

void preprocess2()
{
    // random seed, can be fixed for rerun of experiments
    unsigned int seed = 1777732000; //1777732000; 2315039500; 1137751000; //3981461068; chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(seed);
    cerr << seed << '\n';

    cutoff = 0.0;
    test_force_constants.resize(number_of_terms);
    for (int i = 0; i < number_of_terms; i++)
    {
        if (rng() & 1)
        {
            if (rng() & 1)
            {
                test_force_constants[i] = random_step(rng, cutoff, upper_limit);
            }
            else
            {
                test_force_constants[i] = random_step(rng, -upper_limit, -cutoff);
            }
        }
        else
        {
            test_force_constants[i] = random_step(rng, -cutoff, cutoff);
        }
    }
    angles.clear();
    for (int i = -180 + angle_step; i <= 180; i += angle_step)
    {
        angles.push_back(i);
    }
    test_points.resize(angles.size());
    for (int i = 0; i < (int) angles.size(); i++)
    {
        test_points[i] = force_field_calculate(test_force_constants, angles[i]);
    }

    test_force_constants = find_with_simplicity_accuracy_trading(angles, test_points, number_of_terms);
}

int summary2()
{
    /*
    printf("-----------------------------------------------------\n");
    printf("Origi constant: ");
    for (int i = 0; i < number_of_terms; i++)
    {
        printf("%.6lf ", test_force_constants[i]);
    }
    printf("\n");
    printf("Force constant: ");
    for (int i = 0; i < number_of_terms; i++)
    {
        printf("%.6lf ", force_constants[i]);
    }
    printf("\n");
    */
    bool match = true;
    for (int i = 0; i < number_of_terms; i++)
	{
        match &= (abs(force_constants[i] - (abs(test_force_constants[i]) >= cutoff ? test_force_constants[i] : 0.0)) <= 1.0e-5);
	}
    //if (not match)
    if (1)
    {
        for (int i = 0; i < number_of_terms; i++)
        {
            cerr << "WRONG " << i << " : " << test_force_constants[i] << ' ' << force_constants[i] << ' '  << (abs(force_constants[i] - (abs(test_force_constants[i]) >= cutoff ? test_force_constants[i] : 0.0)) <= 1.0e-6) << '\n';
        }
        cerr << '\n';
    }
    return match;
}

int main()
{
	for (int test = 0; test < 0/*number_of_tests*/; test++)
    {
        // preprocessing
        preprocess1();

        // GLOBAL MINIMUM FINDING
        // simulated annealing
        force_constants = simulated_annealing(rmse, number_of_terms, upper_limit);
        //force_constants = threshold_accepting(angles, test_points, number_of_terms, 2000, 1.0, upper_limit, cutoff);

        // grading
        result += summary1();
    }
    //printf("Test 1 - Cutoff. Correct: %d/%d\n", result, number_of_tests);

    result = 0;
    for (int test = 0; test < number_of_tests; test++)
    {
        // preprocessing
        preprocess2();
        ///*
        ofstream plot;
        plot.open("plot.txt");
        force_constants.resize(number_of_terms);
        for (double i = -2; i <= 2; i += 1e-1)
        {
            for (double j = -2; j <= 2; j += 1e-1)
            {
                force_constants[1] = i;
                force_constants[2] = j;
                plot << i << " " << j << " " << rmse_with_cutoff_and_simplicity_accuracy_trading(force_constants) << "\n";
            }
        }
        plot.close();
        //*/
        // GLOBAL MINIMUM FINDING
        // simulated annealing
        force_constants = simulated_annealing(rmse_with_cutoff_and_simplicity_accuracy_trading, number_of_terms, upper_limit);

        // exact algorithm
        //force_constants = find_with_simplicity_accuracy_trading(angles, test_points, number_of_terms);
        // grading
        result += summary2();
    }
    printf("Test 2 - Simplicity-accuracy tradding. Correct: %d/%d", result, number_of_tests);
}


