// smoothmedian.cpp
// c++ source code for creating smoothmedian.mex file in Matlab as follows:
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
    mxArray *xarr = mxDuplicateArray (prhs[0]);
    double *x = (double *) mxGetData (xarr);
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
    
    // Get data dimensions
    int ndims = (int) mxGetNumberOfDimensions (xarr);
    const mwSize *sz = mxGetDimensions (xarr);
    int m = sz[0];
    int n = sz[1];
    int N = mxGetNumberOfElements (xarr);

    // Prepare output vector
    mwSize dims[2] = {1,n};
    plhs[0] = mxCreateNumericArray (2, dims, 
                mxDOUBLE_CLASS, 
                mxREAL); 
    double *M = (double *) mxGetData(plhs[0]);
    
    // Calculate basic statistics for each column of the data
    float mid = 0.5 * m;
    std::vector<double> xmin (n);
    std::vector<double> xmax (n);
    std::vector<double> range (n);
    for (int k = 0; k < n ; k++) {
        sort(x + k * m, x + k * m + m);
        M[k] = x[k * m + int (mid)];    // Median when m is odd
        if ( mid == int (mid) ) {      
            M[k] += x [k * m + int (mid) - 1];
            M[k] *= 0.5;                // Median when m is even
        }
        xmin[k] = x[k * m];
        xmax[k] = x[k * m + m - 1];
        range[k] = xmax[k] - xmin[k];
    }
        
    // Declare variables that we update in the loop with math assignment operators
    double T;
    double v;
    double U;
    
    // Loop through each column of the data and apply smoothing
    int MaxIter = 500;
    for (int k = 0; k < n ; k++) {
        if (nrhs < 3) {
            Tol = range[k] * 1e-4; 
        }
        // Using the (ordinary) median as the starting value, find the smoothed median
        // Set initial bracket bounds
        double a = xmin[k]; 
        double b = xmax[k];
        // Set initial value of free parameter to the midrange
        double p = M[k];    
        // Start iterations
        for (int Iter = 0; Iter < MaxIter ; Iter++) {
            // Break from iterations if the range of the x values is zero
            if (range[k] == 0) {
                break;
            }   
            T = 0;
            v = 0;
            U = 0;
            for (int j = 0; j < m ; j++) {
                double xj = x [k * m + j];
                for (int i = 0; i < j ; i++) {
                    double xi = x [k * m + i];
                    // Calculate first derivative (T)
                    double D = pow (xi - p, 2) + pow (xj - p, 2);
                    double R = sqrt(D);
                    T += (2 * p - xi - xj) / R;
                    // Calculate second derivative (U)
                    U += pow (xi - xj, 2) * R / pow (D, 2);
                }
            }
            // Compute Newton step (fast quadratic convergence but unreliable)
            double step = T / U;
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
                double nwt = p - step;
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
    }
    
    // Clean up
    mxDestroyArray (xarr);
    
    return;
    
}