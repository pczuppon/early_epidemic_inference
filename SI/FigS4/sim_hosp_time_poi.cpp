#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 sim_hosp_time_poi.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

////////////////////////////////////////
// parameters
////////////////////////////////////////

// secondary infections distribution - negative binomial
#define R0 1.3                       // mean number of secondary infections     

// Infectiousness distribution - gamma distribution (mean = 6 days, SD = 2.5 days)
#define shape_inf 6.6
#define scale_inf 0.833

// hospitalization time - gamma distribution
#define shape_hosp 31.
#define scale_hosp 0.463

// hospitalization probability
#define p_hosp 0.029

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
        gsl_rng_set(r,r_ind);                    // setting the seed
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
    double hosp_min = 1000.;               // minimal hosp time -- 1000 is a dummy number!
    double t=0;                             // time
    double dt = 0.1;                        // time step !!! WARNING: Only values of format 10^(-k) allowed
    int length = (int)hosp_min/dt;         // array length
    int return_value = 0;
    int index = 0;
    
    int I_work[length+1] = {0};             // number of infecteds working vector
    int I[(int)hosp_min+1] = {0};          // number of infecteds at a certain time (days!)
    
    // Initialization = 1 infected individual at time 0
    I_work[0] = 1;
    I[0] = 1;
    
    // Random distributions necessary throughout
    unsigned int gsl_ran_poisson(const gsl_rng * r, double mu);  // mu = mean
    double gsl_ran_gamma(const gsl_rng * r, double a, double b);        // a = shape, b = scale
    unsigned int gsl_ran_binomial(const gsl_rng * r, double p, unsigned int n); // p = prob, n = size
    
    int work_ind = 0;
    while (t < hosp_min)
    {   
        // check for hosp of individuals added at time t
        int n_hosp = gsl_ran_binomial(r, p_hosp, I_work[work_ind]);
        
        // update minimal hospitalization time
        for (int i = 1; i <= n_hosp; i++)
        {
            double hosp_time = t + gsl_ran_gamma(r,shape_hosp,scale_hosp);
            if (hosp_time < hosp_min)
            {
                hosp_min = hosp_time;
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
            
            if (I_ind < (int)hosp_min)
            {
                I[I_ind]++;
            }            
        }
        
        t += dt;
        work_ind++;
    }
       
    // If epidemic happened increase success counter
    if (hosp_min < 1000)
        return_value = 1;
    
    // If epidemic happened, save trajectory and final pop size separately
    if (return_value ==1)
    {
        //ofstream file ("trajectory_hosp_id_" + to_string(r_ind) + ".txt", ios::app); 
        //for (int ind = 0; ind < work_ind; ind++)
        //{
        //    file << ind;
        //    file << ",";
        //    file << accumulate(I,I+ind+1,0);
        //    file << "\n";
        //}
        //file.close();
        
        ofstream file2 ("size_time.txt", ios::app);
        file2 << hosp_min;
        file2 << ",";
        file2 << accumulate(I_work,I_work+work_ind,0);
        file2 << ",";
        file2 << accumulate(I_work,I_work+index,0);
        file2 << "\n";
        file2.close();
    }

    return(return_value);
}





