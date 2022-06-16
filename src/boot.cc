// boot.cc
// c++ source code for creating boot.oct file using mkoctfile in Octave 
//
// boot.oct is a function file for generating balanced bootstrap sample indices
//
// bootsam = boot (n, B, u)
//
// INPUT VARIABLES
// n (double) is the number of rows (of the data vector)
// B (double) is the number of bootstrap resamples
// u (boolean) false (for bootstrap) or true (for bootknife)
//
// OUTPUT VARIABLE
// bootsam (short integer, int16) is an n x B matrix of bootstrap resamples
//
// Uniform random numbers are generated using the Mersenne Twister 19937 generator
//
// Author: Andrew Charles Penn (2022)

#include <cstdlib>
#include <stdio.h>
#include <random>
#include <octave/oct.h>

DEFUN_DLD (boot, args, , 
           "Function file (boot.oct) for generating balanced bootstrap sample indices \n\n"\
           "bootsam = boot (n, B, u) \n\n"\
           "INPUT VARIABLES \n"\
           "n (double) is the number of rows (of the data vector) \n"\
           "B (double) is the number of bootstrap resamples \n"\
           "u (boolean) false (for bootstrap) or true (for bootknife) \n\n"\
           "OUTPUT VARIABLE \n"\
           "bootsam (short integer, int16) is an n x B matrix of bootstrap resamples \n\n"\
           "Uniform random numbers are generated using the Mersenne Twister 19937 generator \n\n"\
           "Author: Andrew Charles Penn (2022)")
 {

    if (args.length () != 3)
        print_usage ();
    
    // Declare variables
    const short int n = args(0).int_value ();
    const int B = args(1).int_value ();
    bool u = args(2).bool_value ();
    dim_vector dv (n, B); 
    int16NDArray bootsam (dv);       // Array of bootstrap sample indices
    int d;                           // Counter for cumulative sum calculations
    int c[n];                        // Counter for each of the sample indices
    for (int i = 0; i < n ; i++) {
        c[i] = B;
    }
    int N = n * B;                   // Total counts of all sample indices
    int k;                           // Variable to store random number
    bool LOO = false;                // Leave-one-out (LOO) flag (ignored if u is false)
    int r = -1;                      // Sample index for LOO (ignored if u is false)
    int m = 0;                       // Counter for LOO sample index r (ignored if u is false) 

    // Initialize random number generator
    std::random_device seed;
    std::mt19937 rng(seed());
    std::uniform_real_distribution<float> dist(0,1);
    
    // Perform balanced sampling
    for (int b = 0; b < B ; b++) {     
        if (u) {
            r = b - (b / n) * n;
        }
        for (int i = 0; i < n ; i++) {
            if (u) {
                if (c[r] != N) { 
                    m = c[r];
                    c[r] = 0;
                    LOO = true;
                }
            }
            k = dist(rng) * (N - m - 1); 
            d = c[0];
            for (int j = 0; j < n ; j++) { 
                if (k < d) {
                    bootsam(i, b) = j + 1;
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
