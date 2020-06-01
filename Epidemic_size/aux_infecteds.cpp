#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// g++ -std=c++11 aux_infecteds.cpp `pkg-config --libs gsl` command for compiling

using namespace std;

////////////////////////////////////////
// parameters
////////////////////////////////////////

// secondary infections distribution - negative binomial
#define kappa 0.91                      // dispersion parameter
#define R0 3.41                      // mean number of secondary infections
#define p_neg kappa/(kappa+R0)       // success probability for the negative-binomial      

// Infectiousness distribution - gamma distribution (mean = 6 days, SD = 2.5 days)
#define shape_inf 5.76
#define scale_inf 1.04167         

// Time
#define tfin 100.
#define dt 0.01

// Extinction probability
#define ext 0.314

// Random number generation with Mersenne Twister
gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);

////////////////////////////////////////
// Main code
///////////////////////////////////////

int RUN(int,int);

int main(int argc,char *argv[])
{
    int length = int(tfin/dt)+1;
    double I[length] = {0};
    
    double gsl_ran_gamma_pdf(double x, double a, double b); // a = shape, b = scale 
    
    I[0] = 1.;
    
    while (I[0] <= 1)
    {
       
        double t_ind = 0;
        int ind = 1;
        int j = 0;
    
        ofstream file ("infectivity_" + to_string((int)I[0]) + ".txt", ios::app);
    
        while (ind < length)
        {
            t_ind += dt;
            double inf = 0;
            int j = 0;
            double t = 0;
            
            while (j < ind)
            {
                inf += I[j] * dt * R0 * gsl_ran_gamma_pdf(t_ind-t, shape_inf, scale_inf);
                t += dt;
                j++;
            }
            
            I[ind] = inf * (1.-pow(ext,accumulate(I,I+ind,0.)))/(1.-ext);
            
            file << I[ind];
            file << "\n";
            
            ind++;
        }
        
        file.close();
        
        I[0]++;
    
    }       
    return(0);
}


