// smoothmedian.cpp
// c++ source code for creating smoothmedian.mex file as follows:
//
// mex -compatibleArrayDims smoothmedian.cpp
//
// smoothmedian.mex is a function file for calculating a smooth 
// version of the median [1]
//
// M = smoothmedian (x)
// M = smoothmedian (x, dim)
// M = smoothmedian (x, dim, Tol)      
//
// INPUT VARIABLES
// x (double) is the data vector or matrix
// dim (double) is the dimension (1 for columnwise, 2 for rowwise)
// Tol (double) sets the step size that will to stop optimization 
//
// OUTPUT VARIABLE
// M (double) is vector of the smoothed median
//
// If x is a vector, find the univariate smoothed median (M) of x.
// If x is a matrix, compute the univariate smoothed median value
// for each column and return them in a row vector (default). The 
// argument dim defines which dimension to operate along. Arrays  
// of more than two dimensions are not currently supported. Tol 
// configures the stopping criteria, in terms of the absolute 
// change in the step size. By default, Tol = RANGE * 1e-4.
//
// The smoothed median is a slightly smoothed version of the ordinary 
// median and is an M-estimator that is both robust and efficient:
//
// | Asymptotic            | Mean |    Median  |    Median  |
// | properties            |      | (smoothed) | (ordinary) |
// |-----------------------|------|------------|------------|
// | Breakdown point       | 0.00 |      0.341 |      0.500 |
// | Pitman efficacy       | 1.00 |      0.865 |      0.637 |
//
// Smoothing the median is achieved by minimizing the following
// objective function:
//
//      S (M) = sum (((x(i) - M).^2 + (x(j) - M).^2).^ 0.5)
//             i < j
// 
// where i and j refers to the indices of the Cartesian product 
// of each column of x with itself. 
//
// This function minimizes the above objective function by finding 
// the root of the first derivative using a fast, but reliable, 
// Newton-Bisection hybrid algorithm. The tolerance (Tol) is the 
// maximum step size that is acceptable to break from optimization.
//
// Bibliography:
// [1] Brown, Hall and Young (2001) The smoothed median and the
//      bootstrap. Biometrika 88(2):519-534
//
// Author: Andrew Charles Penn (2022)



#include "mex.h"
#include <stdlib.h>
#include <vector>
using namespace std;

void mexFunction (int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[]) 
{
    
    // Input variables
    if (nrhs < 1) {
        mexErrMsgTxt("function requires at least 1 input argument");
    }
    
    // Input variable declaration
    double *x = (double *) mxGetData (prhs[0]);
    short int dim;
    if (nrhs < 2) {
        dim = 1;
    } else {
        dim = *(mxGetPr (prhs[1]));
    }
    double Tol;
    if (nrhs > 2) {
        Tol = *(mxGetPr (prhs[2]));
    }

    // Get data dimensions and prepare output vector
    int ndims = (int) mxGetNumberOfDimensions (prhs[0]);
    const mwSize *sz = mxGetDimensions (prhs[0]);
    if (sz[0] == 1) {
        dim = 2;
    }
    int m, n;
    if (dim == 1) {
        m = sz[0];
        n = sz[1];
        plhs[0] = mxCreateDoubleMatrix (1, n, mxREAL);
    } else if (dim == 2) {
        n = sz[0];
        m = sz[1]; 
        plhs[0] = mxCreateDoubleMatrix (n, 1, mxREAL);
    }
    int N = mxGetNumberOfElements (prhs[0]);
    float mid = 0.5 * m;
    double *M = (double *) mxGetData(plhs[0]);
    
    // Declare temporary variables needed for the optimization step
    vector<double> xvec;
    xvec.reserve (m);
    double a, b, range, T, v, U, D, R, step, nwt, p;
    
    // Loop through the data and apply smoothing to the median
    int MaxIter = 500;
    for (int k = 0; k < n ; k++) {

        // Copy the next row/column of the data to temporary vector and sort it
        if (dim == 1) {
            for (int j = 0; j < m ; j++) xvec.push_back ( x[k * m + j] );
        } else if (dim == 2) { 
            for (int j = 0; j < m ; j++) {int i = j * n; xvec.push_back ( x[i + k] );};
        }
        sort(xvec.begin(), xvec.end());
        
        // Set the (ordinary) median as the starting value
        M[k] = xvec[int(mid)];           // Median when m is odd
        if ( mid == int(mid) ) {      
            M[k] += xvec [int(mid) - 1];
            M[k] *= 0.5;                 // Median when m is even
        }
        
        // Set initial bracket bounds
        a = xvec[0];                     // Minimum
        b = xvec[m - 1];                 // Maximimum
        
        // Calculate range (and set stopping criteria if not specified)
        range = b - a;                   // Range
        if (nrhs < 3) {
            Tol = range * 1e-4; 
        }
        
        // Set initial value of free parameter to the midrange
        p = M[k];   
        
        // Start iterations
        for (int Iter = 0; Iter < MaxIter ; Iter++) {
            
            // Break from iterations if the range of the x values is zero
            if (range == 0) {
                break;
            }   
            
            // Perform computations for the current data vector
            T = 0;
            v = 0;
            U = 0;
            for (int j = 0; j < m ; j++) {
                                
                if ( !mxIsFinite(xvec[j]) ) {
                    mexErrMsgTxt ("x cannot contain NaN or Inf");
                }
                
                for (int i = 0; i < j ; i++) {
                    
                    // Calculate first derivative (T)
                    D = pow (xvec[i] - p, 2) + pow (xvec [j] - p, 2);
                    R = sqrt(D);
                    T += (2 * p - xvec[i] - xvec [j]) / R;
                    
                    // Calculate second derivative (U)
                    U += pow (xvec[i] - xvec [j], 2) * R / pow (D, 2);
                    
                }
            }
            
            // Compute Newton step (fast quadratic convergence but unreliable)
            step = T / U;
            
            // Evaluate convergence
            if (abs (step) < Tol) {
                break; // Break from optimization when converged to tolerance 
            } else {
                
                // Update bracket bounds for Bisection
                if (T < -Tol) {
                    a = p;
                } else if (T > +Tol) {
                    b = p;
                }
                
                // Preview new value of the smoothed median
                nwt = p - step;
                
                // Choose which method to use to update the smoothed median
                if (nwt > a && nwt < b) {
                    // Use Newton step if it is within bracket bounds
                    p = nwt;
                } else {
                    // Compute Bisection step (slow linear convergence but very safe)
                    p = 0.5 * (a + b);
                }
            }
            
            if (Iter == MaxIter) {
                mexWarnMsgTxt ("Warning: Root finding failed to reach the specified tolerance");
            }
            
        }
        
        // Assign parameter value that optimizes the objective function for the smoothed median
        M[k] = p;
        
        // Clear the temporary vector for the next cycle of the loop
        xvec.clear(); 
    }
    
    return;
    
}