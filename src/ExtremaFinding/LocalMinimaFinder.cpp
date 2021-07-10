#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "LocalMinimaFinder.h"
#ifndef _NR_
#define _NR_
#include "../numerical_recipes/nr.h" // numerical recipes
#endif

#define random_step(rng,min,max) min + (double)rng() / UINT_MAX * (max - min) 

function<double(vector<double>&)> objective_function;

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

double norm(
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

void clustering(
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
                    if (clustered[ind2].value - 1e-5 < C_twiddle[ind1].value && norm(C_twiddle[ind1], clustered[ind2]) <= r)
                    {
                        clustered.push_back(C_twiddle[ind1]);
                        C_twiddle[ind1].clustered = true;
                        clusterable = true;
                        break;
                    }
                }
            }
        }
    }
}

DP objective_function_wrapper(Vec_I_DP &x) {
    vector<double> y;
    y.resize(x.size());
    for (int i = 0; i < (int) x.size(); i++)
    {
        y[i] = x[i];
    }
    return objective_function(y);
}

/*
Implementation follows Armstrong, J.C., Favorite, J.A. Using a derivative-free 
optimization method for multiple solutions of inverse transport problems. 
Optim Eng 17, 105–125 (2016). https://doi.org/10.1007/s11081-015-9306-x

Relevant variables' name all follow original paper. Here, we assume the function
is well defined over the considered domain.
*/
vector<Point> MLSL(
    int n, 
    double l, 
    function<double(vector<double>&)> target_function)
{
    vector<Point> X_star;
    objective_function = target_function;
    
    // random seed, can be fixed for rerun of experiments
    unsigned int seed = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 rng(seed);

    // Powell's algorithm, using Numerical Recipes, provided by Prof. Alexandrov
    // setup minimizer
    DP FTOL = 1.0e-5;
    const int IMAXSTEP = 10000;
    int iter;
    DP fret;
    Vec_DP p(0.0, n); // variables
    Mat_DP xi(n, n);

    // C0: INITIALIZE
    // 1. Choose termination conditions for MLSL and MADS
    // 2. Choose N_twiddle > 0, N_hat > 0, 0 < gamma < 1, alpha > 0
    bool i;
    int 
        N_twiddle = 10, 
        N_hat = 20;
    double 
        gamma = 0.5,
        alpha = 0.5;
    // 3. C <- empty, X_star <- empty, X_hat <- empty, N_0 <- 0, k <- 0
    vector<vector<double> > C;
    vector<Point> X_hat;
    vector<double> x;
    vector<pair<double, int> > values_in_C;
    int N = 0, k = 0;
    x.resize(n);

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
        set<Point> C_bar;
        while (j < N_twiddle)
        {
        // (a) C_bar <- empty, j <- j + 1
            C_bar.clear();
            j++;
        // (b) Draw N_hat points with uniform distribution on Theta, 
        // evluate objective funciton vaues at samples, and add sample points 
        // to C.
            for (int ind1 = 0; ind1 < N_hat; ind1++)
            {
                for(int ind2 = 0; ind2 < n; ind2++)
                {
                    x[ind2] = random_step(rng, -l, l);
                }
                C.push_back(x);
                values_in_C.push_back({objective_function(x), C.size() - 1});
            }
    
        // (c) Choose [gamma * k * N_hat] points in C with lowest objective 
        // function values and place them in C_bar.
            sort(values_in_C.begin(), values_in_C.end());
            int gamma_k_N_hat = gamma * k * N_hat;
            for (int ind = 0; ind < gamma_k_N_hat; ind++)
            {
                C_bar.insert({
                    C[values_in_C[ind].second], 
                    values_in_C[ind].first, 
                    false
                });
            }

        // (d) If all points in C_bar have finite objective function values, then go C1-2
        // (e) If j < N_twiddle, then go to C1-1, else go to C1-2.
        }
    
    // C1-2: Processing
        // (a) Remove all points whose objective function values are infinity from C_bar.
        // (b) If C_bar = empty, go to C3.
        // (c) Choose all distinct points in C_bar and place them in C_twiddle.
        int ind2 = (int) C_twiddle.size();
        C_twiddle.resize(C_twiddle.size() + C_bar.size());
        for (set<Point>::iterator ind1 = C_bar.begin(); ind1 != C_bar.end(); ind1++)
        {
            C_twiddle[ind2++] = (*ind1);
        }
        // (d) N_k <- N_(k-1) + j * N_hat; compute r_k
        N = N + j * N_hat;
        double r = 1.0/sqrt(M_PI) * pow(tgamma(1.0 + n/2.0) * power(2.0 * l, n) * (1 - pow(alpha, 1.0/(N - 1))), 1.0/n);
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
                vector<Point> x_hat;
                Point x_star = {C_twiddle[ind1].x, 0.0, false};
                
                // Local minimization using Powell's algorithm
                for (int ind2 = 0; ind2 < n; ind2++)
                {
                    p[ind2] = x_star.x[ind2];
                }
                for (int ind2 = 0; ind2 < n; ind2++)
                {
                    for (int ind3 = 0; ind3 < n; ind3++)
                    {
                        xi[ind2][ind3] = (ind2 == ind3 ? p[ind2] : 0.0);
                    }
                }
                
                NR::powell(p, xi, FTOL, IMAXSTEP, iter, fret, objective_function_wrapper);
                for (int ind2 = 0; ind2 < n; ind2++)
                {
                    x_star.x[ind2] = p[ind2];
                }
                x_star.value = objective_function(x_star.x);
                x_star.clustered = true;
        // (b) If x_star not in X_star, then (x_hat <- x_star, add x_star to X_star, i <- 1)
        //      else (x_hat <- x_twiddle, add x_twiddle to X_hat.
                bool in_X_star = false;
                for (int ind2 = 0; ind2 < (int) X_star.size(); ind2++)
                {
                    if (norm(X_star[ind2], x_star) <= 1.0e-4)
                    {
                        in_X_star = true;
                        break;
                    }
                }
                if (not in_X_star)
                {
                    X_star.push_back(x_star);
                    x_hat.push_back(x_star);
                    i = 1;
                }
                else
                {
                    X_hat.push_back(C_twiddle[ind1]);
                    C_twiddle[ind1].clustered = true;
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
    return X_star;
}