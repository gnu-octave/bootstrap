// smoothmedian.cpp
// c++ source code for creating smoothmedian.mex file as follows:
//
// mex smoothmedian.cpp
//
// smoothmedian.mex is a function file for calculating a smooth version of the
// median [1]
//
// M = smoothmedian (X)
// M = smoothmedian (X, DIM)
// M = smoothmedian (X, DIM, TOL)
//
// INPUT VARIABLES
// X (double) is the data vector or matrix.
// DIM (double) is the dimension (1 for columnwise, 2 for rowwise).
// TOL (double) sets the step size that will to stop optimization.
//
// OUTPUT VARIABLE
// M (double) is a scalar or vector of the smoothed median(s)
//
// If X is a vector, find the univariate smoothed median (M) of X. If X is a
// matrix, compute the univariate smoothed median value for each column and
// return them in a row vector (default). The argument DIM defines which
// dimension to operate along. Arrays of more than two dimensions are not
// currently supported. TOL configures the stopping criteria, in terms of the
// absolute change in the step size. By default, TOL = RANGE * 1e-4.
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
// Smoothing the median is achieved by minimizing the following objective
// function:
//
//      S (M) = sum (((X(i) - M).^2 + (X(j) - M).^2).^ 0.5)
//             i < j
// 
// where i and j refers to the indices of the Cartesian product of each column
// of X with itself.
//
// This function minimizes the above objective function by finding the root of
// the first derivative using a fast, but reliable, Newton-Bisection hybrid
// algorithm. The tolerance (TOL) is the maximum step size that is acceptable to
// break from optimization. Data (X) values that are equal to NaN are ignored.
//
// The smoothing works by slightly reducing the breakdown point of the median.
// Bootstrap confidence intervals using the smoothed median have good coverage
// for the ordinary median of the population distribution and can be used to
// obtain second order accurate intervals with Studentized bootstrap and
// calibrated percentile bootstrap methods [1]. When the population distribution
// is thought to be strongly skewed, coverage errors can be reduced by improving
// symmetry through appropriate data transformation. Unlike kernel-based
// smoothing approaches, bootstrapping smoothmedian does not require explicit
// choice of a smoothing parameter or a probability density function.
//
// Bibliography:
// [1] Brown, Hall and Young (2001) The smoothed median and the
//      bootstrap. Biometrika 88(2):519-534
//
// Author: Andrew Charles Penn (2022)


#include "mex.h"         // for mex functions
#include <vector>        // for vectpr function
#include <cmath>         // for pow function
#include <algorithm>     // for nth_element function
using namespace std;


void mexFunction (int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[]) 
{

    // Input variables
    if ( nrhs < 1 ) {
        mexErrMsgTxt ("At least one input argument is required.");
    }
    // First input argument (x)
    if ( !mxIsClass (prhs[0], "double") ) {
        mexErrMsgTxt ("The first input argument (X) must be of type double");
    }
    if ( mxIsComplex (prhs[0]) ) {
        mexErrMsgTxt ("The first input argument (X) cannot contain an imaginary part");
    }
    double *x = (double *) mxGetData (prhs[0]);
    // Second input argument (dim)
    short int dim;
    if ( nrhs > 1 && !mxIsEmpty (prhs[1]) ) {
        if ( mxGetNumberOfElements (prhs[1]) > 1 ) {
            mexErrMsgTxt ("The second input argument (DIM) must be scalar");
        }
        if ( !mxIsClass (prhs[1], "double") ) {
            mexErrMsgTxt ("The second input argument (DIM) must be of type double");
        }
        if ( mxIsComplex (prhs[1]) ) {
            mexErrMsgTxt ("The second input argument (DIM) cannot contain an imaginary part");
        }
        dim = *(mxGetPr (prhs[1]));
    } else {
        dim = 1;
    }
    if ( dim != 1 && dim != 2) {
        mexErrMsgTxt ("The second input argument (DIM) must be 1 (column-wise) or 2 (row-wise)");
    }
    // Third input argument (Tol)
    double Tol;
    if ( nrhs > 2 && !mxIsEmpty (prhs[2]) ) {
        if ( mxGetNumberOfElements (prhs[2]) > 1 ) {
            mexErrMsgTxt ("The third input argument (TOL) must be scalar");
        }
        if ( !mxIsClass (prhs[2], "double") ) {
            mexErrMsgTxt ("The third input argument (TOL) must be of type double");
        }
        if ( mxIsComplex (prhs[2]) ) {
            mexErrMsgTxt ("The third input argument (TOL) cannot contain an imaginary part");
        }
        Tol = *(mxGetPr (prhs[2]));
        if ( !mxIsFinite (Tol) ) {
            mexErrMsgTxt ("The third input argument (TOL) cannot be NaN or Inf");    
        }
        if ( Tol < 0 ) {
            mexErrMsgTxt ("The third input argument (TOL) must be a positive value");
        }
    }

    // Get data dimensions and prepare output vector
    int ndims = (int) mxGetNumberOfDimensions (prhs[0]);
    const mwSize *sz = mxGetDimensions (prhs[0]);
    if ( sz[0] == 1 ) {
        dim = 2;
    }
    int m, n, l;
    if ( dim == 1 ) {
        m = sz[0];
        n = sz[1];
        plhs[0] = mxCreateDoubleMatrix (1, n, mxREAL);
    } else if ( dim == 2 ) {
        n = sz[0];
        m = sz[1]; 
        plhs[0] = mxCreateDoubleMatrix (n, 1, mxREAL);
    }
    int N = mxGetNumberOfElements (prhs[0]);
    double *M = (double *) mxGetData(plhs[0]);

    // Declare temporary variables needed for the optimization step
    vector<double> xvec;
    xvec.reserve (m);
    double a, b, mid, range, S, T, U, D, R, step, nwt;

    // Loop through the data and apply smoothing to the median (maximum 25 iterations)
    int MaxIter = 24;
    for ( int k = 0; k < n ; k++ ) {

        // Copy the next row/column of the data to temporary vector
        if ( dim == 1 ) {
            for ( int j = 0; j < m ; j++ ) xvec.push_back ( x[k * m + j] );
        } else if ( dim == 2 ) { 
            for ( int j = 0; j < m ; j++ ) {int i = j * n; xvec.push_back ( x[i + k] );};
        }

        // Omit NaN values and calculate the length of the resulting vector
        xvec.erase (remove_if (xvec.begin(), xvec.end(), mxIsNaN), xvec.end());
        l = xvec.size ();

        // Set the (ordinary) median as the starting value
        mid = 0.5 * l;
        nth_element (xvec.begin(), xvec.begin() + int (mid), xvec.end());
        // After running nth_element, none of the elements in xvec preceding the
        // nth are greater than it, and none of the elements after it are less.
        if ( mid == int (mid) ) {
            // Median when l is even
            M[k] = xvec[mid];
            M[k] += *max_element (xvec.begin(), xvec.begin() + mid);
            M[k] *= 0.5;
        } else {
            // Median when l is odd
            mid = int (mid);
            M[k] = xvec[mid];
        }

        // Set initial bracket bounds to the minimum and maximum data values
        a = *min_element (xvec.begin(), xvec.begin() + mid);
        b = *max_element (xvec.begin() + mid, xvec.end());

        // Calculate range
        range = b - a;
        
        // Set stopping criteria (if Tol is not already specified)
        if ( nrhs < 3 || mxIsEmpty (prhs[2]) ) {
            Tol = range * 1e-4; 
        }  

        // Start iterations
        for ( int Iter = 0; Iter <= MaxIter ; Iter++ ) {
            
            // Break from iterations if the distance between the bracket bounds 
            // < Tol since the smoothed median will be equal to the median 
            if ( range <= Tol ) {
                break;
            }   

            // Calculate derivatives of the objective function for Newton-Raphson method
            //S = 0;
            T = 0;
            U = 0;
            for ( int j = 0; j < l ; j++ ) {        
                for ( int i = 0; i < j ; i++ ) {
                    D = pow (xvec[i] - M[k], 2) + pow (xvec [j] - M[k], 2);
                    R = sqrt (D);
                    // Objective function (S)
                    //S += R;
                    if ( D != 0 ) {
                        // First derivative (T)
                        T += (2 * M[k] - xvec[i] - xvec [j]) / R;
                        // Second derivative (U)
                        U += pow (xvec[i] - xvec [j], 2) * R / pow (D, 2);
                    }
                }
            }

            // Compute Newton step (fast quadratic convergence but unreliable)
            step = T / U;

            // Evaluate convergence
            if ( abs (step) <= Tol ) { 
                break; // Break from optimization when converged to tolerance 
            } else {
                // Update bracket bounds for Bisection method
                if ( step < 0 ) {
                    a = M[k] + Tol;
                } else if ( step > 0 ) {
                    b = M[k] - Tol;
                }
                // Update the range with the distance between the bracket bounds
                range = b - a;
                // Preview new value of the smoothed median
                nwt = M[k] - step;
                // Choose which method to use to update the smoothed median
                if ( nwt > a && nwt < b ) {
                    // Use Newton step if it is within bracket bounds
                    M[k] = nwt;
                } else {
                    // Compute Bisection step (slow linear convergence but very safe)
                    M[k] = 0.5 * (a + b);
                }
            }

            if ( Iter == MaxIter ) {
                if (dim == 1) {
                    mexPrintf ("warning: Root finding failed to reach tolerance for column %d \n", k+1);
                } else {
                    mexPrintf ("warning: Root finding failed to reach tolerance for row %d \n", k+1);
                }
            }
            
        }

        // Clear the temporary vector for the next cycle of the loop
        xvec.clear(); 
    }

    return;

}
