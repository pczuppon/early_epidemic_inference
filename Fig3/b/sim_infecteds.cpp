#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 sim_infecteds.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

////////////////////////////////////////
// parameters
////////////////////////////////////////

// R0 distribution
#define R0 2.9                      // mean number of secondary infections

// Infectiousness distribution - gamma distribution (mean = 6 days, SD = 2.5 days)
#define shape_inf 6.6
#define scale_inf 0.833

// Repetitions
#define repeats 10000                              // number of repetitions

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

////////////////////////////////////////
// Main code
////////////////////////////////////////

int RUN(int,int);

int main(int argc,char *argv[])
{
    int tfin = atoi(argv[1]);
    int r_ind = 0;                                    // repetition index
    int success = 0;                                  // epidemic counter
    
    while (success < repeats)
    {
        gsl_rng_set(r,r_ind);                         // setting the seed
        success += RUN(tfin,success);                 // Run stochastic simulation
        r_ind++;
    }
        
    return(0);
}

////////////////////////////////////////
// Stochastic simulation
////////////////////////////////////////
int RUN(int tfin,int counter)
{
    double t=0;                         // time
    double dt = 0.1;                    // time step !!! WARNING: Only values of format 10^(-k) allowed
    int length = (tfin+1)/dt;           // array length
    int return_value = 0;
    
    int I_work[length] = {0};           // number of infecteds working vector
    int I[tfin+1] = {0};                // number of infecteds at a certain time (days!)
    int I_detected[tfin+1] = {0};       // number of detectable (days!)
    
    // Initialization = 1 infected individual at time 0
    I_work[0] = 1;
    I[0] = 1;
    
    // Random distributions necessary throughout
    unsigned int gsl_ran_poisson(const gsl_rng * r, double mu);  // mu = mean
    double gsl_ran_gamma(const gsl_rng * r, double a, double b);        // a = shape, b = scale
    
    int work_ind = 0;
    while (t < (double)tfin)
    {   
        // no. of offspring of individuals added at time t
        int offspring = 0;
        for (int i = 1; i <= I_work[work_ind]; i++)
        {
            offspring = gsl_ran_poisson(r,R0);
            
            if (offspring > 0)
            {                
                // loop over offsprings from individuals that were added at time t
                for (int i = 1; i <= offspring; i++)
                {   
                    double inf_time = t + gsl_ran_gamma(r,shape_inf,scale_inf);
                    
                    // IMPORTANT THAT dt = 10^(-k)! otherwise adapt this rounding formula!
                    int inf_time_ind = roundf(inf_time/dt);
                    int I_ind = ceilf(inf_time);
                    
                    if (inf_time_ind < length)
                    {
                        I_work[inf_time_ind]++;
                    }
                    
                    if (I_ind <= tfin)
                    {
                        I[I_ind]++;
                    }            
                }
            }
        }            
        
        t += dt;
        work_ind++;
    }
       
    // If epidemic happened increase success counter
    if (accumulate(I+tfin-5,I+tfin,0) > 0)
        return_value = 1;
    
    // If epidemic happened, save trajectory
    if (return_value ==1)
    {
        ofstream file ("trajectory_tfin_" + to_string(tfin) + "_id_" + to_string(counter) + ".txt", ios::app); 
        for (int ind = 0; ind <= tfin; ind++)
        {
            file << ind;
            file << ",";
            file << I[ind];
            file << "\n";
        }
        file.close();
    }

    return(return_value);
}

