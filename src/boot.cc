// boot.cc
// c++ source code for creating boot.oct file using mkoctfile in Octave 
//
// boot.oct is a function file for generating balanced bootstrap sample indices
//
// USAGE
// bootsam = boot (n, nboot)
// bootsam = boot (n, nboot, u)
//
// INPUT VARIABLES
// n (short integer, int16) is the number of rows (of the data vector)
// nboot (integer, int32) is the number of bootstrap resamples
// u (boolean) for unbiased: false (for bootstrap) or true (for bootknife)
//
// u is an optional input argument. The default is false.
//
// If u is true then the sample index for omission in each bootknife resample 
// is selected systematically. If the remaining number of bootknife resamples 
// is not divisible by the sample size (n), then the sample index omitted is  
// selected randomly. 
//
// OUTPUT VARIABLE
// bootsam (short integer, int16) is an n x nboot matrix of bootstrap resamples
//
// Uniform random numbers are generated using the Mersenne Twister 19937 generator
//
// Author: Andrew Charles Penn (2022)

#include <octave/oct.h>
#include <random>

DEFUN_DLD (boot, args, , 
           " Function file (boot.oct) for generating balanced bootstrap sample indices \n"\
           " \n"\
           " USAGE \n"\
           " bootsam = boot (n, nboot) \n"\
           " bootsam = boot (n, nboot, u) \n"\
           " \n"\
           " INPUT VARIABLES \n"\
           " n (short integer, int16) is the number of rows (of the data vector) \n"\
           " nboot (integer, int32) is the number of bootstrap resamples \n"\
           " u (boolean) for unbiased: false (for bootstrap) or true (for bootknife) \n"\
           " \n"\
           " u is an optional input argument. The default is false. \n"\
           " \n"\
           " If u is true then the sample index for omission in each bootknife resample \n"\
           " is selected systematically. If the remaining number of bootknife resamples \n"\
           " is not divisible by the sample size (n), then the sample index omitted is \n"\
           " selected randomly. \n"\
           " \n"\
           " OUTPUT VARIABLE \n"\
           " bootsam (short integer, int16) is an n x nboot matrix of bootstrap resamples \n"\
           " \n"\
           " Uniform random numbers are generated using the Mersenne Twister 19937 generator \n"\
           " \n"\
           " Author: Andrew Charles Penn (2022)")
 {

    // Input variables
    if (args.length () < 2) {
        print_usage ();
    }
    const short int n = args(0).int_value ();
    const int nboot = args(1).int_value ();
    bool u;
    if (args.length () < 3) {
        u = false;
    } else {
        u = args(2).bool_value ();
    }
    
    // Declare variables
    dim_vector dv (n, nboot); 
    int16NDArray bootsam (dv);       // Array of bootstrap sample indices
    int d;                           // Counter for cumulative sum calculations
    int c[n];                        // Counter for each of the sample indices
    for (int i = 0; i < n ; i++) {   
        c[i] = nboot;                // Set each element in c to nboot
    }
    int N = n * nboot;               // Total counts of all sample indices
    int k;                           // Variable to store random number
    bool LOO = false;                // Leave-one-out (LOO) flag for the current bootstrap iteration (remains false if u is false)
    int r = -1;                      // Sample index for LOO (remains -1 and is ignored if u is false)
    int m = 0;                       // Counter for LOO sample index r (remains 0 if u is false) 

    // Create pointer so that we can more rapidly access elements of bootsam
    octave_int<short> *ptr = bootsam.fortran_vec ();
    
    // Initialize random number generator
    std::random_device rd;
    std::seed_seq seed {rd(), rd(), rd(), rd()};
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> distr (0, n - 1);
    
    // Perform balanced sampling
    for (int b = 0; b < nboot ; b++) { 
        if (u) {    
            if ((b / n) == (nboot / n)) {
                r = distr (rng);      // random
            } else {
                r = b - (b / n) * n;  // systematic
            }
        }
        for (int i = 0; i < n ; i++) {
            if (u) {
                if (c[r] < N) {   // Only LOO if sample index r doesn't account for all remaining sampling counts
                    m = c[r];
                    c[r] = 0;
                    LOO = true;
                }
            }
            std::uniform_int_distribution<int> distk (0, N - m - 1);
            k = distk (rng); 
            d = c[0];
            for (int j = 0; j < n ; j++) { 
                if (k < d) {
                    *(ptr + b * n + i) = j + 1;
                    c[j] -= 1;
                    N -= 1;
                    break;
                } else {
                    d += c[j + 1];
                }
            }
            if (LOO == true) {
                c[r] = m;
                m = 0;
                LOO = false;
            }
        }   
    }    

    return octave_value (bootsam);
} 
