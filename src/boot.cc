// Function file for generating balanced bootstrap sample indices
//
// bootsam = boot (n, B, w)
//
// INPUT VARIABLES
// n (double) is the number of rows (of data vector)
// B (double) is the number of bootstrap resamples
// w (double) should be set to 0 (for bootstrap) or 1 (for bootknife)
//
// OUTPUT VARIABLE
// bootsam (short integer, int16) is an n x B matrix of bootstrap resamples
//
// Author: Andrew Charles Penn (2022)

#include <cstdlib>
#include <stdio.h>
#include <octave/oct.h>

DEFUN_DLD (boot, args, , "bootsam = boot (n, B, w)")

 {

    if (args.length () != 3)
        print_usage ();
    
    // Initialise
    const short int n = args(0).int_value ();
    const int B = args(1).int_value ();
    const short int w = args(2).int_value ();
    dim_vector dv (n, B); 
    int16NDArray bootsam (dv);
    int r;
    int c[n];
    for (int i = 0; i < n ; i++) {
        c[i] = B;
    }
    int d[n];
    int k;
    float m;          // Must be float to avoid integer division
    float N = n * B;  // Must be float to avoid integer division

    // Perform balanced sampling
    for (int b = 0; b < B ; b++) 
    {   
        r = b - (b / n) * n;
        for (int i = 0; i < n ; i++) {
            if (c[r] != N && w > 0) {
                // Choose which row of the data to exclude for this bootknife sample
                m = c[r];
                c[r] = 0;
            } else {
                m = 0;
            } 
            k = rand() * (N - m) / RAND_MAX; 
            d[0] = c[0];
            for (int j = 0; j <= n ; j++) {
                if (j < n) {
                    if (k >= d[j]) {
                        d[j+1] = d[j] + c[j+1];
                    } else {
                        bootsam(i,b) = j+1;
                        c[j] -= 1;
                        break;
                    }   
                } else {
                    bootsam(i,b) = j-1;
                    c[j-1] -= 1;    
                }
            }
            if (c[r] != N && w > 0) {
                c[r] = m;
            }
            N -= 1;
        }      
    }  

    return octave_value (bootsam);
} 
