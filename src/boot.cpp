// boot.cpp
// c++ source code for creating boot.mex file using mex as follows:
// 
// mex -compatibleArrayDims boot.cpp
//
// boot.mex is a function file for generating balanced bootstrap sample indices
//
// USAGE
// BOOTSAM = boot (N, NBOOT)
// BOOTSAM = boot (X, NBOOT)
// BOOTSAM = boot (..., NBOOT, UNBIASED)
// BOOTSAM = boot (..., NBOOT, UNBIASED, WEIGHTS)
// BOOTSAM = boot (..., NBOOT, UNBIASED, WEIGHTS, SEED)
//
// INPUT VARIABLES
// N (double) is the number of rows (of the data vector)
// X (double) is a data vector intended for resampling
// NBOOT (double) is the number of bootstrap resamples
// UNBIASED (boolean) for unbiased resampling: false (for bootstrap) or true
//   (for bootknife)
// WEIGHTS (double) is a weight vector of length N
// SEED (double) is a seed used to initialise the pseudo-random number generator
//
// OUTPUT VARIABLE
// BOOTSAM (double) is an N x NBOOT matrix of sample indices (N) or resampled
//   data (X)
//
// NOTES
// UNBIASED is an optional input argument. The default is false. If UNBIASED is
// true, then the sample index for omission in each bootknife resample is
// selected systematically. When the remaining number of bootknife resamples is
// not divisible by the sample size (N), then the sample index omitted is
// selected randomly. 
// WEIGHTS is an optional input argument. If WEIGHTS is empty or not provided,
// the default is a vector of each element equal to NBOOT (i.e. uniform
// weighting). Each element of WEIGHTS is the number of times that the
// corresponding index is represented in bootsam. Therefore, the sum of WEIGHTS
// should equal N * NBOOT. Note that the mex function compiled from this source
// code is not thread safe. Below is an example of a line of code one can run in
// Octave/Matlab before attempting parallel operation of boot.mex in order to
// ensure that the initial random seeds of each thread are unique:
//
// In Octave:
// >> pararrayfun(nproc, @boot, 1, 1, false, [], 1:nproc)
// In Matlab:
// >> ncpus = feature('numcores'); parfor i = 1:ncpus; boot (1, 1, false, [], i); end;
//
// Author: Andrew Charles Penn (2022)


#include "mex.h"
#include <vector>
using namespace std;


void mexFunction (int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[]) 
{

    // Input variables
    if ( nrhs < 2 ) {
        mexErrMsgTxt ("function requires at least 2 input arguments");
    }
    // First input argument (n or x)
    double *x = (double *) mxGetData (prhs[0]);
    int n = mxGetNumberOfElements (prhs[0]);
    bool isvec;
    if ( n > 1 ) {
        const mwSize *sz = mxGetDimensions (prhs[0]);
        if ( sz[0] > 1 && sz[1] > 1 ) {
            mexErrMsgTxt ("the first input argument must be a scalar or a vector");
        }
        isvec = true;
    } else {
        isvec = false;
        n = *(mxGetPr (prhs[0]));
        if ( !mxIsFinite (n) ) {
            mexErrMsgTxt ("the first input argument cannot be NaN or Inf");    
        }
        if ( mxIsComplex (prhs[0]) ) {
            mexErrMsgTxt ("the first input argument cannot contain an imaginary part");
        }
    }
    if ( !mxIsClass (prhs[0], "double") ) {
        mexErrMsgTxt ("the first input argument must be of type double");
    }
    // Second input argument (nboot)
    const int nboot = *(mxGetPr (prhs[1]));
    if ( mxGetNumberOfElements (prhs[1]) > 1 ) {
        mexErrMsgTxt ("the second input argument (nboot) must be scalar");
    }
    if ( !mxIsClass (prhs[1], "double") ) {
        mexErrMsgTxt ("the second input argument (nboot) must be of type double");
    }
    if ( mxIsComplex (prhs[1]) ) {
        mexErrMsgTxt ("the second input argument (nboot) cannot contain an imaginary part");
    }
    if ( nboot <= 0 ) {
        mexErrMsgTxt ("the second input argument (nboot) must be a positive integer");
    }
    if ( !mxIsFinite (nboot) ) {
        mexErrMsgTxt ("the second input argument (nboot) cannot be NaN or Inf");    
    }
    // Third input argument (u)
    bool u;
    if ( nrhs < 3 ) {
        u = false;
    } else {
        if (mxGetNumberOfElements (prhs[2]) > 1 || !mxIsClass (prhs[2], "logical")) {
            mexErrMsgTxt ("the third input argument (u) must be a logical scalar value");
        }
        u = *(mxGetLogicals (prhs[2]));
    }
    // Fourth input argument (w) - error checking is handled later (see below)
    // Fifth input argument (s)
    unsigned int seed;
    if ( nrhs > 4 ) {
        if ( mxGetNumberOfElements (prhs[4]) > 1 ) {
            mexErrMsgTxt ("the fifth input argument (s) must be a scalar value");
        }
        seed = *(mxGetPr(prhs[4]));
        if ( !mxIsFinite (seed) ) {
            mexErrMsgTxt ("the fifth input argument (s) cannot be NaN or Inf");    
        }
        srand (seed);
    }  

    // Output variables
    if (nlhs > 1) {
        mexErrMsgTxt ("the function (boot) can only return a single output argument");
    }

    // Declare variables
    mwSize dims[2] = {n, nboot};
    plhs[0] = mxCreateNumericArray(2, dims, 
                mxDOUBLE_CLASS, 
                mxREAL);           // Prepare array for bootstrap sample indices
    long long int N = n * nboot;   // Total counts of all sample indices
    long long int k;               // Variable to store random number
    long long int d;               // Counter for cumulative sum calculations
    vector<long long int> c;       // Counter for each of the sample indices
    c.reserve (n);
    if ( nrhs > 3 && !mxIsEmpty (prhs[3]) ) {
        // Assign user defined weights (counts)
        if ( !mxIsClass (prhs[3], "double") ) {
            mexErrMsgTxt ("the fourth input argument (weights) must be of type double");
        }
        if ( mxIsComplex (prhs[3]) ) {
            mexErrMsgTxt ("the fourth input argument (weights) cannot contain an imaginary part");
        }
        double *w = (double *) mxGetData (prhs[3]);
        if ( mxGetNumberOfElements (prhs[3]) != n ) {
            mexErrMsgTxt ("the fourth input argument (weights) must be a vector of length n");
        }
        long long int s = 0; 
        for ( int i = 0; i < n ; i++ )  {
            if ( !mxIsFinite (w[i]) ) {
                mexErrMsgTxt ("the fourth input argument (weights) cannot contain NaN or Inf");    
            }
            if ( w[i] < 0 ) {
                mexErrMsgTxt ("the fourth input argument (weights) must contain only positive integers");
            }
            c.push_back (w[i]); // Set each element in c to the specified weight    
            s += c[i];
        }
        if ( s != N ) {
            mexErrMsgTxt ("the elements in the forth input argument (weights) must sum to n * nboot");
        }
    } else {
        // Assign weights (counts) for uniform sampling
        for ( int i = 0; i < n ; i++ ) {   
            c.push_back (nboot);      // Set each element in c to nboot
        }
    }
    bool LOO = false;     // Leave-one-out (LOO) flag
    long long int m = 0;  // Counter for LOO sample index r
    int r = -1;           // Sample index for LOO

    // Create pointer so that we can access elements of bootsam (i.e. plhs[0])
    double *ptr = (double *) mxGetData(plhs[0]);

    // Perform balanced sampling
    for ( int b = 0; b < nboot ; b++ ) { 
        if (u == true) {    
            if ( (b / n) == (nboot / n) ) {
                r = rand () / (RAND_MAX / n + 1);  // random
            } else {
                r = b - (b / n) * n;               // systematic
            }
        }
        for ( int i = 0; i < n ; i++ ) {
            if (u == true) {
                // Only LOO if sample index r doesn't account for all remaining 
                // sampling counts
                if (c[r] < N) {
                    m = c[r];
                    c[r] = 0;
                    LOO = true;
                }
            }
            k = rand () / (RAND_MAX / (N - m) + 1); 
            d = c[0];
            for ( int j = 0; j < n ; j++ ) { 
                if ( k < d ) {
                    if (isvec) {
                        ptr[b * n + i] = x[j];
                    } else {
                        ptr[b * n + i] = j + 1;
                    }
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

    return;

}