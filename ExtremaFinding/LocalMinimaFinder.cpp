#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "LocalMinimaFinder.h"

long long power(
    int a, 
    int n)
{
    if (n == 0)
    {
        return 1;
    }
    if (n == 1)
    {
        return a;
    }
    long long b = power(a, n>>1);
    b *= b;
    if (n&1)
    {
        b *= a;
    }
    return b;
}

struct MachineNumber
{
    int a, b, b0;

    MachineNumber()
    {
        a = 1, b = 0, b0 = 0;
    }

    void increase()
    {
        if(a == 1)
        {
            a = 2;
        }
        else if(a == 2)
        {
            a = 5;
        }
        else if(a == 5)
        {
            a = 1, b++;
        }
    }

    void decrease()
    {
        if(a == 2)
        {
            a = 1;
        }
        else if(a == 5)
        {
            a = 2;
        }
        else if(a == 1)
        {
            a = 5, b--;
        }
    }

    void inverse()
    {
        if (a == 2)
        {
            a = 5;
        }
        else if (a == 5)
        {
            a = 2;
        }

        b = -b - (a == 1 ? 0 : 1);
    }
};

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

/*
Implemetation follows Audet, Charles & Le Digabel, Sébastien & Tribes, 
Christophe. (2019). The Mesh Adaptive Direct Search Algorithm for Granular 
and Discrete Variables. SIAM Journal on Optimization. 29. 1164-1189. 
doi:10.1137/18M1175872.

Relevant variables' name all follow original paper.
*/
struct LocalMinimaFinder::MADS
{
    bool calculated;
    int n;
    double a_t;
    vector<double> x0, x;
    vector<MachineNumber> Delta, delta;
    vector<int> rho;
    function <double(vector<double>)> objective_function;

    MADS(
        vector<double> x_0, 
        function <double(vector<double>)> target_function)
    {
        x0 = x_0;
        objective_function = target_function;
        x = x0;
        calculated = false;

        n = x0.size();
        Delta.resize(n);
        delta.resize(n);
    }

    vector<double> MADS_search()
    {
        // if the result has been calculated before, return immediately
        if (calculated) return x;

        // initial scaling
        for (int i = 0; i < n; i++)
        {
            double alpha = x0[i]/10;
            while (alpha < 1.0)
            {
                Delta[i].b--,
                alpha *= 10;
            }
            if (alpha <= 1.5)
            {
                Delta[i].a = 1;
            }
            else if (alpha <= 2.5)
            {
                Delta[i].a = 2;
            }
            else if (alpha <= 7.5)
            {
                Delta[i].a = 5;
            }
            else 
            {
                Delta[i].a = 1, Delta[i].b++;
            }
            Delta[i].b0 = Delta[i].b;
        }

        while (true)
        {
            // set delta^k = 10^(b^k - |b^k - b^0|)
            for (int i = 0; i < n; i++)
            {
                delta[i].b = Delta[i].b - abs(Delta[i].b - Delta[i].b0);
            }

            // set rho^k = diag(delta^k)^-1 * Delta^k
            for (int i = 0; i < n; i++)
            {
                MachineNumber p = delta[i];
                p.inverse(), 
                p.a = Delta[i].a;
                p.b += Delta[i].b;

                //rho[i] = p;
            }

            bool successful = false;
            // Search step: test finitely main trial points on the mesh M^k (Omitted)

            // TODO: Poll step: Test P^k = {x^k + diag(delta^k) * d : d in D^k} 
            // where D^k in Z^(n * n) is a positive spanning set satisfying
            // -psi^k <= d <= psi^k for every d in D^k

            // if the search or pool step was successful at a trial point t 
            // in M^k, the set x^(k+1) = t and 
            // - Delta^(k+1) = increase(Delta^k)
            // if |d_i|/p_i^k > a_t, or if delta_i^k < delta_i^0 and psi_i^k >
            // (psi_l^k) for some l in N, and for every index i, and where d in
            // Z^n is the direction that led too success: 
            // x^(k+1) = x^k + diag(delta^k) * d.
            // - Delta^(k+1) = decrease(Delta^k) otherwise
            if (successful)
            {
               // x = t;

                // TODO: check for sufficient condition
                bool satisfying = false;

                if (satisfying)
                {
                    for (int i = 0; i < n; i++)
                    {
                        Delta[i].increase();
                    }
                }
            }
            // else set x^(k+1) = x^k and Delta^(k+1) = decrease(Delta^k)
            else
            {
                for (int i = 0; i < n; i++)
                {
                    Delta[i].decrease();
                }
            }
        }

        calculated = true;
        return x;
    }
};

double norm(
    LocalMinimaFinder::Point a, 
    LocalMinimaFinder::Point b)
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
    int 
        N_twiddle = 10,
        N_hat = 20;
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
                MADS new_run = MADS(C_twiddle[ind1].x, objective_function);
                vector<double> x_star = new_run.MADS_search();
                vector<Point> x_hat;
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