#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 sim_sequencing.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

////////////////////////////////////////
// parameters
////////////////////////////////////////

// secondary infections distribution - negative binomial
#define kappa 0.57                  // dispersion parameter
#define R0 1.1                       // mean number of secondary infections
#define p_neg kappa/(kappa+R0)      // success probability for the negative-binomial      

// Infectiousness distribution - gamma distribution (mean = 6 days, SD = 2.5 days)
#define shape_inf 6.6
#define scale_inf 0.833

// hospitalization time - gamma distribution
#define shape_detect 12.
#define scale_detect 7./12.

// hospitalization probability
#define p_detect 0.042

// Repetitions
#define repeats 10000                              // number of repetitions

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

////////////////////////////////////////
// Main code
////////////////////////////////////////

int RUN(int);
// int argc,char *argv[]
int main()
{
    int r_ind = 0;                                    // repetition index
    int success = 0;                                  // epidemic counter
        
    while (success < repeats)
    {
        gsl_rng_set(r,r_ind);                  // setting the seed
        success += RUN(r_ind);                 // Run stochastic simulation
        r_ind++;
    }
        
    return(0);
}

////////////////////////////////////////
// Stochastic simulation
////////////////////////////////////////
int RUN(int r_ind)
{
    double samp_min = 1000.;               // minimal sampling time -- 1000 is a dummy number!
    double t=0;                             // time
    double dt = 0.1;                        // time step !!! WARNING: Only values of format 10^(-k) allowed
    int length = (int)samp_min/dt;         // array length
    int return_value = 0;
    int index = 0;
    
    int I_work[length+1] = {0};             // number of infecteds working vector
    int I[(int)samp_min+1] = {0};          // number of infecteds at a certain time (days!)
    
    // Initialization = 1 infected individual at time 0
    I_work[0] = 1;
    I[0] = 1;
    
    // Random distributions necessary throughout
    //unsigned int gsl_ran_negative_binomial(const gsl_rng * rrest, double p, double n);  // p = success probability, n = sample size
    unsigned int gsl_ran_poisson(const gsl_rng * r, double mu);  // mu = mean
    double gsl_ran_gamma(const gsl_rng * rrest, double a, double b);        // a = shape, b = scale
    double gsl_ran_lognormal(const gsl_rng * r, double zeta, double sigma); // zeta = mean, sigma = standard deviation
    unsigned int gsl_ran_binomial(const gsl_rng * rbin, double p, unsigned int n); // p = prob, n = size
    
    int work_ind = 0;
    while (t < samp_min)
    {   
        // check for sampling of individuals added at time t
        int n_samp = gsl_ran_binomial(r, p_detect, I_work[work_ind]);
        
        // update minimal sampling time
        for (int i = 1; i <= n_samp; i++)
        {
            double samp_time = t + gsl_ran_gamma(r,shape_detect,scale_detect);
            if (samp_time < samp_min)
            {
                samp_min = samp_time;
                index = work_ind;
            }
        }
        
        // no. of offspring of individuals added at time t
        int offspring = 0;
        for (int i = 1; i <= I_work[work_ind]; i++)
        {
            offspring += gsl_ran_poisson(r,R0);
        }
        
        // loop over offsprings from individuals that were added at time t
        for (int i = 1; i <= offspring; i++)
        {   
            double inf_time = t+gsl_ran_gamma(r,shape_inf,scale_inf);
            // IMPORTANT THAT dt = 10^(-k)! otherwise adapt this rounding formula!
            int inf_time_ind = roundf(inf_time/dt);
            int I_ind = ceilf(inf_time);
            
            if (inf_time_ind < length)
            {
                I_work[inf_time_ind]++;
            }
            
            if (I_ind < (int)samp_min)
            {
                I[I_ind]++;
            }            
        }
        
        t += dt;
        work_ind++;
    }
       
    // If epidemic happened increase success counter
    if (samp_min < 1000)
        return_value = 1;
    
    // If epidemic happened, save final pop size and time
    if (return_value ==1)
    {        
        ofstream file ("size_time.txt", ios::app);
        file << samp_min;
        file << ",";
        file << accumulate(I_work,I_work+work_ind,0);
        file << "\n";
        file.close();
    }

    return(return_value);
}





