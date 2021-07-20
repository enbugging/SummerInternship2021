#define _USE_MATH_DEFINES
#include <bits/stdc++.h>

using namespace std;

#include "src/ExtremaFinding/GlobalMinimumFinder.h"
#include "src/ExtremaFinding/LocalMinimaFinder.h"

#define random_step(rng,min,max) min + (double)rng() / UINT_MAX * (max - min)

double 
    cutoff,
    upper_limit = 3.0;
int
	number_of_terms = 6,
	angle_step = 10,
    number_of_tests = 1000,
    result = 0;
vector<double>
	angles,
	test_points,
	test_weights,
	test_force_constants,
    force_constants;

void preprocess1()
{
    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    cutoff = 0.1;
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
    for (int i = 0; i < 360; i += angle_step)
    {
        angles.push_back(i);
    }
    test_points.resize(angles.size());
    for (int i = 0; i < (int) angles.size(); i++)
    {
        test_points[i] = force_field_calculate(test_force_constants, angles[i]);
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
        match &= (abs(force_constants[i] - (abs(test_force_constants[i]) >= cutoff ? test_force_constants[i] : 0.0)) <= 1.0e-6);
	}
    return match;
}

void preprocess2()
{
    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
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
    for (int i = 0; i < 360; i += angle_step)
    {
        angles.push_back(i);
    }
    test_points.resize(angles.size());
    for (int i = 0; i < (int) angles.size(); i++)
    {
        test_points[i] = force_field_calculate(test_force_constants, angles[i]);
    }

    test_force_constants = simulated_annealing_brute_force(angles, test_points, number_of_terms, 5000, 1.0, upper_limit);
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
        match &= (abs(force_constants[i] - (abs(test_force_constants[i]) >= cutoff ? test_force_constants[i] : 0.0)) <= 1.0e-6);
	}
    return match;
}

int main()
{
	for (int test = 0; test < number_of_tests; test++)
    {
        // preprocessing
        preprocess1();

        // GLOBAL MINIMUM FINDING
        // simulated annealing
        force_constants = simulated_annealing(angles, test_points, number_of_terms, 5000, 1.0, upper_limit, cutoff);
        //force_constants = threshold_accepting(angles, test_points, number_of_terms, 5000, 50.0, 3.0);

        // grading
        result += summary1();
    }
    printf("Test 1 - Cutoff. Correct: %d/%d\n", result, number_of_tests);

    result = 0;
    for (int test = 0; test < 0/*number_of_tests*/; test++)
    {
        // preprocessing
        preprocess2();

        // GLOBAL MINIMUM FINDING
        // simulated annealing
        //force_constants = simulated_annealing(angles, test_points, number_of_terms, 5000, 1.0, upper_limit, cutoff);
        force_constants = simulated_annealing_with_simplicity_accuracy_trading(angles, test_points, number_of_terms, 5000, 1.0, upper_limit, cutoff);
        //force_constants = threshold_accepting(angles, test_points, number_of_terms, 5000, 50.0, 3.0);
            
        // grading
        result += summary2();
    }
    printf("Test 2 - Simplicity-accuracy tradding. Correct: %d/%d", result, number_of_tests);
}
