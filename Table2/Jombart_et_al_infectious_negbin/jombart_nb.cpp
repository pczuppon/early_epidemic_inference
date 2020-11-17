#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 jombart_nb.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

////////////////////////////////////////
// parameters
////////////////////////////////////////

#define kappa 0.57                      // dispersion parameter
#define R0 2.9                      // mean number of secondary infections
#define p_neg kappa/(kappa+R0)       // success probability for the negative-binomial      

// Infectiousness distribution - gamma distribution (mean = 6 days, SD = 2.5 days)
#define shape_inf 6.6
#define scale_inf 0.833

// hospitalization time - gamma distribution
#define shape_hosp 31.0
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
    // Random distributions necessary throughout
    unsigned int gsl_ran_negative_binomial(const gsl_rng * r, double p, double n);  // p = success probability, n = sample size
    double gsl_ran_gamma(const gsl_rng * r, double a, double b);        // a = shape, b = scale
    unsigned int gsl_ran_geometric(const gsl_rng * r, double p);        // p = success probability
    
    // Draw hospitalization time
    double hosp_time = gsl_ran_gamma(r, shape_hosp,scale_hosp);
    
    // initialize epidemic               
    double t=0;                             // time
    double dt = 0.1;                        // time step !!! WARNING: Only values of format 10^(-k) allowed
    int length = roundf(hosp_time/dt);         // array length
    int length_int = roundf(hosp_time);
       
    int I_work[length+1] = {0};             // number of infecteds working vector
    int I[length_int+1] = {0};          // number of infecteds at a certain time (days!)
    int I_removed[length_int+1] = {0};  // removed individuals
        
    // Initialization = geometric number of cases at t=0
    I_work[0] = gsl_ran_geometric(r,p_hosp);
    I[0] = I_work[0];
    
    int work_ind = 0;
    while (t < hosp_time)
    {   
              
        // no. of offspring of individuals added at time t
        int offspring = 0;
        for (int i = 1; i <= I_work[work_ind]; i++)
        {
            offspring = gsl_ran_negative_binomial(r,p_neg,kappa);
        
            if (offspring == 0)
            {
                int rem_ind = ceilf(work_ind*dt);
                I_removed[rem_ind]++;
            }
            
            else
            {
                double max_inf_time = 0.;
                
                // loop over offsprings from individuals that were added at time t
                for (int i = 1; i <= offspring; i++)
                {   
                    double inf_time = t+gsl_ran_gamma(r,shape_inf,scale_inf);
                    max_inf_time = max(inf_time,max_inf_time);
                    
                    // IMPORTANT THAT dt = 10^(-k)! otherwise adapt this rounding formula!
                    int inf_time_ind = roundf(inf_time/dt);
                    int I_ind = ceilf(inf_time);
                    
                    if (inf_time_ind < length)
                    {
                        I_work[inf_time_ind]++;
                    }
                    
                    if (I_ind < length_int)
                    {
                        I[I_ind]++;
                    }   
                }
                
                int I_rem_ind = ceilf(max_inf_time);

                if (I_rem_ind <= length_int)
                {
                    I_removed[I_rem_ind]++;
                }
            }
        }        
        t += dt;
        work_ind++;
    }
       
    // Increase success counter
    int return_value = 1;
    
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
        int ind = ceilf(work_ind*dt);
        
        ofstream file2 ("size.txt", ios::app);
        file2 << accumulate(I,I+ind,0);
        file2 << ",";
        file2 << accumulate(I_removed,I_removed+ind,0);
        file2 << "\n";
        file2.close();
    }

    return(return_value);
}





