#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "../nomad/interfaces/CInterface/NomadStdCInterface.h"

#include "LocalMinimaFinder.h"

long long power(
    int a, 
    int p)
{
    if (p == 0)
    {
        return 1;
    }
    if (p == 1)
    {
        return a;
    }
    long long b = power(a, p>>1);
    b *= b;
    if (p&1)
    {
        b *= a;
    }
    return b;
}

struct LocalMinimaFinder::Point
{
    vector<double> x;
    double value;
    bool clustered;
    bool operator<( const Point& other) const {
        return value < other.value;
    }
};

double LocalMinimaFinder::random_step(
    double min,
    double max)
{
    double f = (double)rng() / UINT_MAX;
    return min + f * (max - min);
}

long double norm_squared(
    vector<double> a, 
    vector<double> b
)
{
    long double sum = 0;
    for (int i = 0; i < min((int) a.size(), (int) b.size()); i++)
    {
        double d = a[i] - b[i];
        sum += d * d;
    }
    return sum;
}

double LocalMinimaFinder::norm(
    Point a, 
    Point b)
{
    long double sum = 0;
    for (int i = 0; i < min((int) a.x.size(), (int) b.x.size()); i++)
    {
        double d = a.x[i] - b.x[i];
        sum += d * d;
    }
    return sqrt(sum);
}

void LocalMinimaFinder::clustering(
    vector<Point>& C_twiddle, 
    vector<Point>& clustered,
    double r)
{
    bool clusterable = true;
    while (clusterable)
    {
        clusterable = false;
        for (int ind1 = 0; ind1 < (int) C_twiddle.size(); ind1++)
        {
            if (not C_twiddle[ind1].clustered)
            {
                for (int ind2 = 0; ind2 < (int) clustered.size(); ind2++)
                {
                    if (clustered[ind2].value < C_twiddle[ind1].value && norm(C_twiddle[ind1], clustered[ind2]) <= r)
                    {
                        clustered.push_back(C_twiddle[ind1]);
                        C_twiddle[ind1].clustered = true;
                        clusterable = true;
                    }
                }
            }
        }
    }
}

bool objective_function_for_Nomad(
    int nb_inputs, 
    double *x, 
    int nb_outputs, 
    double *bb_outputs, 
    bool *count_eval, 
    NomadUserDataPtr data)
{
     bool eval_ok = true;

    vector<double> x_copy(x, x + nb_inputs);
    bb_outputs[0] = LocalMinimaFinder::objective_function(x_copy);

    *count_eval = true;
    return eval_ok;
}
/*
bool LocalMinimaFinder::objective_function_for_Nomad(
    int nb_inputs, 
    double *x, 
    int nb_outputs, 
    double *bb_outputs, 
    bool *count_eval, 
    NomadUserDataPtr data)
{
    bool eval_ok = true;

    vector<double> x_copy(x, x + nb_inputs);
    bb_outputs[0] = objective_function(x_copy);

    *count_eval = true;
    return eval_ok;
}
*/

vector<double> LocalMinimaFinder::MADS(
    vector<double>& x
)
{
    // indispensable parameters to create the problem
    int 
        nb_inputs = x.size(), 
        nb_outputs = 1;
    
    // create Nomad problem
    NomadProblem nomad_pb = createNomadProblem(objective_function_for_Nomad,
                                               nb_inputs,
                                               nb_outputs);
    
    double granularity_values[10] = {0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001,
                                     0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001};
    
    // ix parameters using NOMAD convention

    // fix important parameters
    addNomadParam(nomad_pb, (char*)"BB_OUTPUT_TYPE PB PB OBJ EB");

    // fix some external parameters
    addNomadArrayOfDoubleParam(nomad_pb, (char*)"GRANULARITY",  granularity_values);

    addNomadParam(nomad_pb, (char*)"DISPLAY_DEGREE 2");
    addNomadParam(nomad_pb, (char*)"DISPLAY_STATS EVAL ( SOL ) OBJ CONS_H H_MAX");
    addNomadParam(nomad_pb, (char*)"DISPLAY_ALL_EVAL true");
    addNomadParam(nomad_pb, (char*)"DISPLAY_UNSUCCESSFUL false");

    // for reproducibility
    addNomadValParam(nomad_pb, (char*)"NB_THREADS_OPENMP", 1);

    // and the number of blackbox allowed
    addNomadParam(nomad_pb, (char*)"MAX_BB_EVAL 1000");

    // run problem
    double x0[nb_inputs]; // starting point
    copy(x.begin(), x.end(), x0);
    double x_feas_sol[nb_inputs]; // feasible solution
    memset(x_feas_sol, 0, sizeof(x_feas_sol));
    double x_inf_sol[nb_inputs]; // infeasible solution
    memset(x_inf_sol, 0, sizeof(x_inf_sol));
    double outputs_feas_sol[nb_outputs]; // feasible solution outputs
    memset(outputs_feas_sol, 0, sizeof(outputs_feas_sol));
    double outputs_inf_sol[nb_outputs]; // infeasible solution outputs
    memset(outputs_inf_sol, 0, sizeof(outputs_inf_sol));

    bool exists_feas, exists_infeas = false; // flag which indicates if the solution exists or not

    solveNomadProblem(nomad_pb, 1, x0,
                      &exists_feas, x_feas_sol, outputs_feas_sol,
                      &exists_infeas, x_inf_sol, outputs_inf_sol,
                      NULL);
    
    vector <double> ans(x_feas_sol, x_feas_sol + nb_inputs);
    return ans;
}

/*
Implementation follows Armstrong, J.C., Favorite, J.A. Using a derivative-free 
optimization method for multiple solutions of inverse transport problems. 
Optim Eng 17, 105–125 (2016). https://doi.org/10.1007/s11081-015-9306-x

Relevant variables' name all follow original paper. Here, we assume the function
is well defined over the considered domain.
*/
vector<LocalMinimaFinder::Point> LocalMinimaFinder::MLSL_MADS()
{
    if (calculated)
    {
        return X_star;
    }

    // C0: INITIALIZE
    // 1. Choose termination conditions for MLSL and MADS
    // 2. Choose N_twiddle > 0, N_hat > 0, 0 < gamma < 1, alpha > 0
    bool i;
    int N_hat = 20;
    double 
        gamma = 0.5,
        alpha = 1.0;
    // 3. C <- empty, X_star <- empty, X_hat <- empty, N_0 <- 0, k <- 0
    vector<vector<double> > C;
    vector<Point> X_hat;
    vector<pair<double, int> > values_in_C_bar;
    int N = 0, k = 0;
    C.resize(n);
    values_in_C_bar.resize(n);
    for (int ind = 0; ind < N; ind++)
    {
        C[ind].resize(n);
    }

    while (true)
    {
    // C1: GLOBAL PHASE

    // C1-0: Setting
        // (a) C_twiddle <- empty, j <- 0, k <- k + 1
        vector<Point> C_twiddle;
        int j = 0;
        k++;
        // (b) Go to C1-1.

    // C1-1: Sampling
        // (a) C_bar <- empty, j <- j + 1
        set<Point> C_bar;
        j++;
        // (b) Draw N_hat points with uniform distribution on Theta, 
        // evluate objective funciton vaues at samples, and add sample points 
        // to C.
        for (int ind1 = 0; ind1 < N_hat; ind1++)
        {
            for(int ind2 = 0; ind2 < n; ind2++)
            {
                C[ind1][ind2] = random_step(-l, l);
            }
            values_in_C_bar[ind1] = {objective_function(C[ind1]), ind1};
        }
        // (c) Choose [gamma * k * N_hat] points in C with lowest objective 
        // function values and place them in C_bar.
        int gamma_k_N_hat = gamma * k * N_hat;
        sort(values_in_C_bar.begin(), values_in_C_bar.end());
        for (int ind = 0; ind < gamma_k_N_hat; ind++)
        {
            C_bar.insert({
                C[values_in_C_bar[ind].second], 
                values_in_C_bar[ind].first, 
                false
            });
        }
        // (d) If all points in C_bar have finite objective function values, then go C1-2
        // (e) If j < N_twiddle, then go to C1-1, else go to C1-2.

    // C1-2: Processing
        // (a) Remove all points whose objective function values are infinity from C_bar.
        // (b) If C_bar = empty, go to C3.
        // (c) Choose all distinct points in C_bar and place them in C_twiddle.
        C_twiddle.resize(C_bar.size());
        int ind2 = 0;
        for (set<Point>::iterator ind1 = C_bar.begin(); ind1 != C_bar.end(); ind1++)
        {
            C_twiddle[ind2++] = (*ind1);
        }
        // (d) N_k <- N_(k-1) + j * N_hat; compute r_k
        N = N + j * N_hat;
        double r = 1/sqrt(M_PI) * pow(tgamma(1.0 + n/2.0) * power(l, n) * (1 - pow(alpha, 1.0/(N - 1))), 1.0/n);

        // (e) If X_star != empty, apply the single linkage clustering procedure to points in
        // C_twiddle that have not been assigned to clusters with seed points from X_star U X_hat.
        if (not X_star.empty())
        {
            // the set of seed points
            vector<Point> clustered;
            clustered.resize((int) (X_star.size() + X_hat.size()));
            for (int ind = 0; ind < (int) X_star.size(); ind++)
            {
                clustered[ind] = X_star[ind];
            }
            for (int ind = 0; ind < (int) X_hat.size(); ind++)
            {
                clustered[ind + (int) X_star.size()] = X_hat[ind];
            }
            clustering(C_twiddle, clustered, r);
        }

        // (f) If every point in C_twiddle has been assigned to a cluster, then go to C3, else
        if (C_twiddle.empty())
        {
            break;
        }
        // (C_hat <- 0, i <- 0, go to C2).
        else 
        {
            i = 0;
        }

    // C2: LOCAL PHASE
        for (int ind1 = 0; ind1 < (int) C_twiddle.size(); ind1++)
        {
    // 1. Choose x_twiddle in C_twiddle \ C_hat, and add x_twiddle to C_hat
    // 2. If x_twiddle has not been assigned to a cluster, then 
            if (not C_twiddle[ind1].clustered)
            {
        // (a) Apply MADS algorithm with empty SEARCH step using x_twiddle as a start point to find 
        // a local minimizer x_star.
        
                // TODO: x_star finding usign MADS
                vector<Point> x_hat;
                vector<double>x_star;
        // (b) If x_star not in X_star, then (x_hat <- x_star, add x_star to X_star, i <- 1)
        //      else (x_hat <- x_twiddle, add x_twiddle to X_hat.
                bool in_X_star = false;
                for (int ind2 = 0; ind2 < (int) X_star.size(); ind2++)
                {
                    if (X_star[ind2].x == x_star)
                    {
                        in_X_star = true;
                        break;
                    }
                }
                if (not in_X_star)
                {
                    X_hat.push_back({x_star, objective_function(x_star), true});
                    x_hat.push_back({x_star, objective_function(x_star), true});
                    i = 1;
                }
                else
                {
                    X_hat.push_back(C_twiddle[ind1]);
                    x_hat.push_back(C_twiddle[ind1]);
                }

        // (c) Apply the single linkage clustering procedure to C_twiddle using x_hat as a seed point.
                clustering(C_twiddle, x_hat, r);
    // 3. If size(C_twiddle) > size(C_hat), go to C2.
            }
        }
    // 4. If i != 1, go to C3.
        if (i != 1)
        {
            break;
        }
    // 5. If any MLSLS termination criterion is satisfied, then go to C3. Else, go to C1.
        if (k > 100 || X_star.size() > 100 || X_hat.size() > 100)
        {
            break;
        }
    }

    // C3: FINALIZE
    // 1. If X_star != empty, order points in X_star according to their objective function values.
    // 2. Stop MLSL-MADS
    calculated = true;
    return X_star;
}

/*
#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "LocalMinimaFinder.h"
#include "../nomad/interfaces/CInterface/NomadStdCInterface.h"

long long power(
    int a, 
    int p)
{
    if (p == 0)
    {
        return 1;
    }
    if (p == 1)
    {
        return a;
    }
    long long b = power(a, p>>1);
    b *= b;
    if (p&1)
    {
        b *= a;
    }
    return b;
}

struct LocalMinimaFinder::Point
{
    array<double, n> x;
    double value;
    bool clustered;

    Point(array<int, n> y, double val)
    {
        x = y;
        value = val;
        clustered = false;
    }

    bool operator<( const Point& other) const {
        return value < other.value;
    }
};

double LocalMinimaFinder::random_step(
    double min,
    double max)
{
    double f = (double)rng() / UINT_MAX;
    return min + f * (max - min);
}

long double norm_squared(
    vector<double> a, 
    vector<double> b
)
{
    long double sum = 0;
    for (int i = 0; i < (int) min(a.x.size(), b.x.size()); i++)
    {
        double d = a[i] - b[i];
        sum += d * d;
    }
    return sum;
}

double norm(
    LocalMinimaFinder::Point a, 
    LocalMinimaFinder::Point b)
{
    long double sum = 0;
    for (int i = 0; i < (int) min(a.x.size(), b.x.size()); i++)
    {
        double d = a.x[i] - b.x[i];
        sum += d * d;
    }
    return sqrt(sum);
}

void clustering(
    vector<LocalMinimaFinder::Point>& C_twiddle, 
    vector<LocalMinimaFinder::Point>& clustered,
    double r)
{
    bool clusterable = true;
    while (clusterable)
    {
        clusterable = false;
        for (int ind1 = 0; ind1 < (int) C_twiddle.size(); ind1++)
        {
            if (not C_twiddle[ind1].clustered)
            {
                for (int ind2 = 0; ind2 < (int) clustered.size(); ind2++)
                {
                    if (clustered[ind2].value < C_twiddle[ind1].value && norm(C_twiddle[ind1], clustered[ind2]) <= r)
                    {
                        clustered.push_back(C_twiddle[ind1]);
                        C_twiddle[ind1].clustered = true;
                        clusterable = true;
                    }
                }
            }
        }
    }
}

bool LocalMinimaFinder::objective_function_for_Nomad(
    int nb_inputs, 
    double *x, 
    int nb_outputs, 
    double *bb_outputs, 
    bool *count_eval, 
    NomadUserDataPtr data)
{
    bool eval_ok = true;

    bb_outputs[0] = objective_function(x);

    *count_eval = true;
    return eval_ok;
}

vector<double> LocalMinimaFinder::MADS(
    array<int, n> x0
)
{
    // indispensable parameters to create the problem
    int 
        nb_inputs = n, 
        nb_outputs = 1;
    
    // create Nomad problem
    NomadProblem nomad_pb = createNomadProblem(objective_function,
                                               nb_inputs,
                                               nb_outputs);
    
    array<int, n> granularity_values;
    granularity_values.fill(0.0000001);
    
    // ix parameters using NOMAD convention

    // fix important parameters
    addNomadParam(nomad_pb,"BB_OUTPUT_TYPE PB PB OBJ EB");

    // fix some external parameters
    addNomadArrayOfDoubleParam(nomad_pb, "GRANULARITY",  granularity_values);

    addNomadParam(nomad_pb, "DISPLAY_DEGREE 2");
    addNomadParam(nomad_pb, "DISPLAY_STATS EVAL ( SOL ) OBJ CONS_H H_MAX");
    addNomadParam(nomad_pb, "DISPLAY_ALL_EVAL true");
    addNomadParam(nomad_pb, "DISPLAY_UNSUCCESSFUL false");

    // for reproducibility
    addNomadValParam(nomad_pb, "NB_THREADS_OPENMP", 1);

    // and the number of blackbox allowed
    addNomadParam(nomad_pb, "MAX_BB_EVAL 1000");

    // run problem
    array<int, n> x_feas_sol;
    memset(x_feas_sol, 0, sizeof(x_feas_sol)); // feasible solution

    array<int, n> x_inf_sol;
    memset(x_inf_sol, 0, sizeof(x_inf_sol)); // infeasible solution

    array<int, 1> outputs_feas_sol;
    memset(outputs_feas_sol, 0, sizeof(outputs_feas_sol)); // feasible solution outputs

    array<int, 1> outputs_inf_sol;
    memset(outputs_inf_sol, 0, sizeof(outputs_inf_sol)); // infeasible solution outputs

    bool exists_feas, exists_infeas = false; // flag which indicates if the solution exists or not

    solveNomadProblem(nomad_pb, 1, x0,
                      &exists_feas, x_feas_sol, outputs_feas_sol,
                      &exists_infeas, x_inf_sol, outputs_inf_sol,
                      NULL);
}

vector<LocalMinimaFinder::Point> LocalMinimaFinder::MLSL_MADS()
{
    if (calculated)
    {
        return X_star;
    }

    // C0: INITIALIZE
    // 1. Choose termination conditions for MLSL and MADS
    // 2. Choose N_twiddle > 0, N_hat > 0, 0 < gamma < 1, alpha > 0
    bool i;
    int 
        N_twiddle = 10,
        N_hat = 20;
    double 
        gamma = 0.5,
        alpha = 1.0;
    // 3. C <- empty, X_star <- empty, X_hat <- empty, N_0 <- 0, k <- 0
    Point C[N_hat];
    vector<Point> X_hat;
    pair<double, int> values_in_C_bar[N_hat];
    int N = 0, k = 0;

    while (true)
    {
    // C1: GLOBAL PHASE

    // C1-0: Setting
        // (a) C_twiddle <- empty, j <- 0, k <- k + 1
        vector<Point> C_twiddle;
        int j = 0;
        k++;
        // (b) Go to C1-1.

    // C1-1: Sampling
        // (a) C_bar <- empty, j <- j + 1
        set<Point> C_bar;
        j++;
        // (b) Draw N_hat points with uniform distribution on Theta, 
        // evluate objective funciton vaues at samples, and add sample points 
        // to C.
        for (int ind1 = 0; ind1 < N_hat; ind1++)
        {
            for(int ind2 = 0; ind2 < n; ind2++)
            {
                C[ind1].x[ind2] = random_step(-l, l);
            }
            values_in_C_bar[ind1] = {objective_function(C[ind1].x), ind1};
        }
        // (c) Choose [gamma * k * N_hat] points in C with lowest objective 
        // function values and place them in C_bar.
        int gamma_k_N_hat = gamma * k * N_hat;
        sort(values_in_C_bar.begin(), values_in_C_bar.end());
        for (int ind = 0; ind < gamma_k_N_hat; ind++)
        {
            C_bar.insert(
                Point(
                    C[values_in_C_bar[ind].second], 
                    values_in_C_bar[ind].first
                    ) 
                );
        }
        // (d) If all points in C_bar have finite objective function values, then go C1-2
        // (e) If j < N_twiddle, then go to C1-1, else go to C1-2.

    // C1-2: Processing
        // (a) Remove all points whose objective function values are infinity from C_bar.
        // (b) If C_bar = empty, go to C3.
        // (c) Choose all distinct points in C_bar and place them in C_twiddle.
        C_twiddle.resize(C_bar.size());
        int ind2 = 0;
        for (set<Point>::iterator ind1 = C_bar.begin(); ind1 != C_bar.end(); ind1++)
        {
            C_twiddle[ind2++] = (*ind1);
        }
        // (d) N_k <- N_(k-1) + j * N_hat; compute r_k
        N = N + j * N_hat;
        double r = 1/sqrt(M_PI) * pow(tgamma(1.0 + n/2.0) * power(l, n) * (1 - pow(alpha, 1.0/(N - 1))), 1.0/n);

        // (e) If X_star != empty, apply the single linkage clustering procedure to points in
        // C_twiddle that have not been assigned to clusters with seed points from X_star U X_hat.
        if (not X_star.empty())
        {
            // the set of seed points
            vector<Point> clustered;
            clustered.resize((int) (X_star.size() + X_hat.size()));
            for (int ind = 0; ind < (int) X_star.size(); ind++)
            {
                clustered[ind] = X_star[ind];
            }
            for (int ind = 0; ind < (int) X_hat.size(); ind++)
            {
                clustered[ind + (int) X_star.size()] = X_hat[ind];
            }
            clustering(C_twiddle, clustered, r);
        }

        // (f) If every point in C_twiddle has been assigned to a cluster, then go to C3, else
        if (C_twiddle.empty())
        {
            break;
        }
        // (C_hat <- 0, i <- 0, go to C2).
        else 
        {
            i = 0;
        }

    // C2: LOCAL PHASE
        for (int ind1 = 0; ind1 < (int) C_twiddle.size(); ind1++)
        {
    // 1. Choose x_twiddle in C_twiddle \ C_hat, and add x_twiddle to C_hat
    // 2. If x_twiddle has not been assigned to a cluster, then 
            if (not C_twiddle[ind1].clustered)
            {
        // (a) Apply MADS algorithm with empty SEARCH step using x_twiddle as a start point to find 
        // a local minimizer x_star.
        
                // TODO: x_star finding usign MADS
                array<int, n> x_hat;
        // (b) If x_star not in X_star, then (x_hat <- x_star, add x_star to X_star, i <- 1)
        //      else (x_hat <- x_twiddle, add x_twiddle to X_hat.
                bool in_X_star = false;
                for (int ind2 = 0; ind2 < (int) X_star.size(); ind2++)
                {
                    if (X_star[ind2].x == x_star)
                    {
                        in_X_star = true;
                        break;
                    }
                }
                if (not in_X_star)
                {
                    X_hat.push_back(
                        Point(
                            x_star, 
                            objective_function(x_star)
                            )
                        );
                    x_hat.push_back(
                        Point(
                            x_star, 
                            objective_function(x_star)
                            )
                        );
                    i = 1;
                }
                else
                {
                    X_hat.push_back(
                        Point(
                            C_twiddle[ind1], 
                            objective_function(C_twiddle[ind1])
                            )
                        );
                    x_hat.push_back(
                        Point(
                            C_twiddle[ind1], 
                            objective_function(C_twiddle[ind1])
                            )
                        );
                }

        // (c) Apply the single linkage clustering procedure to C_twiddle using x_hat as a seed point.
                clustering(C_twiddle, x_hat, r);
    // 3. If size(C_twiddle) > size(C_hat), go to C2.
            }
        }
    // 4. If i != 1, go to C3.
        if (i != 1)
        {
            break;
        }
    // 5. If any MLSLS termination criterion is satisfied, then go to C3. Else, go to C1.
        if (k > 100 || X_star.size() > 100 || X_hat.size() > 100)
        {
            break;
        }
    }

    // C3: FINALIZE
    // 1. If X_star != empty, order points in X_star according to their objective function values.
    // 2. Stop MLSL-MADS
    calculated = true;
    return X_star;
}
*/