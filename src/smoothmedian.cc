// smoothmedian.cc
// c++ source code for creating smoothmedian.oct file using mkoctfile 
// in Octave 
//
// smoothmedian.oct is a function file for calculating a smooth 
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

#include <octave/oct.h>

DEFUN_DLD (smoothmedian, args, , 
           " smoothmedian.oct is a function file for calculating a smoothed \n"\
           " version of the median [1] \n"\
           " \n"\
           " M   = smoothmedian (x) \n"\
           " ... = smoothmedian (x, dim) \n"\
           " ... = smoothmedian (x, dim, Tol) \n"\
           " \n"\
           " INPUT VARIABLES \n"\
           " x (double) is the data vector or matrix \n"\
           " dim (double) is the dimension (1 for columnwise, 2 for rowwise) \n"\
           " Tol (double) sets the step size that will to stop optimization \n"\
           " \n"\
           " OUTPUT VARIABLE \n"\
           " M (double) is vector of the smoothed median \n"\
           " \n"\
           " If x is a vector, find the univariate smoothed median (M) of x. \n"\
           " If x is a matrix, compute the univariate smoothed median value \n"\
           " for each column and return them in a row vector. The argument dim \n"\
           " defines which dimension to operate along. Arrays of more than two \n"\
           " dimensions are not currently supported. Tol configures the stopping \n"\
           " criteria, in terms of the change in the step size. \n"\
           " \n"\
           " The smoothed median is a slightly smoothed version of the ordinary \n"\
           " median and is an M-estimator that is both robust and efficient: \n"\
           " \n"\
           " | Asymptotic            | Mean |    Median  |    Median  | \n"\
           " | properties            |      | (smoothed) | (ordinary) | \n"\
           " |-----------------------|------|------------|------------| \n"\
           " | Breakdown point       | 0.00 |      0.341 |      0.500 | \n"\
           " | Pitman efficacy       | 1.00 |      0.865 |      0.637 | \n"\
           " \n"\
           " Smoothing the median is achieved by minimizing the following \n"\
           " objective function: \n\n"\
           "      S (M) = sum (((x(i) - M).^2 + (x(j) - M).^2).^ 0.5) \n"\
           "             i < j \n\n"\
           " \n"\
           " where i and j refers to the indices of the Cartesian product \n"\
           " of each column of x with itself. \n"\
           " \n"\
           " This function minimizes the above objective function by finding \n"\
           " the root of the first derivative using a fast, but reliable, \n"\
           " Newton-Bisection hybrid algorithm. The tolerance (Tol) is the \n"\
           " maximum step size that is acceptable to break from optimization. \n"\
           " \n"\
           " Bibliography: \n"\
           " [1] Brown, Hall and Young (2001) The smoothed median and the \n"\
           "      bootstrap. Biometrika 88(2):519-534 \n"\
           " \n"\
           " Author: Andrew Charles Penn (2022)")
 {

    // Input variable declaration
    Matrix x;
    short int dim;
    double Tol;
        
    // Input variables
    if (args.length () < 1) {
        print_usage ();
    } else {
        x = args(0).matrix_value();
    }
    if (args.length () < 2) {
        dim = 1;
    } else {
        dim = args(1).int_value();
    }
    if (args.length () > 2) {
        Tol = args(2).double_value();
    }
    
    // Check that there are no inf or nan values in the data
    if (x.any_element_is_inf_or_nan ()) {
        error("x cannot contain Inf or NaN values \n");
    }
    
    // If applicable, switch dimension
    if (dim > 1) {
        x = x.transpose ();
    }
      
    // Check data dimensions
    const dim_vector sz = x.dims ();
    int m;
    int n;
    if (sz(0) == 1) {
        // If row vector, switch dimensions
        m = sz(1);
        n = sz(0);
        x = x.transpose ();
    } else {
        m = sz(0);
        n = sz(1);
    }
    
    // Calculate basic statistics for each column of the data
    Matrix xsort = x.sort ();
    float mid = 0.5 * m;
    Matrix M;
    M = xsort.extract (int (mid), 0, int (mid), n - 1);   // Median when m is odd
    if ( mid == int (mid) ) {
        // Median when m is even
        M += xsort.extract (int (mid) - 1, 0, int (mid) - 1, n - 1);
        M *= 0.5;
    }
    Matrix xmin = xsort.extract (0, 0, 0, n - 1);   // Minimum
    Matrix xmax = xsort.extract (m - 1, 0, m - 1, n - 1);   // Maximimum
    Matrix range = (xmax - xmin);   // Range

    // Create pointers so that we can more rapidly access elements of the matrices
    double *ptrX = x.fortran_vec ();
    double *ptrM = M.fortran_vec ();
    double *ptrXMIN = xmin.fortran_vec ();
    double *ptrXMAX = xmax.fortran_vec ();
    double *ptrRANGE = range.fortran_vec ();

    // Declare variables that we update in the loop with math assignment operators
    double T;
    double v;
    double U;
    
    // Loop through each column of the data and apply smoothing
    int MaxIter = 500;
    for (int k = 0; k < n ; k++) {
        if (args.length () < 3) {
            Tol = *(ptrRANGE + k) * 1e-4; 
        }
        // Using the (ordinary) median as the starting value, find the smoothed median
        // Set initial bracket bounds
        double a = *(ptrXMIN + k); 
        double b = *(ptrXMAX + k);
        // Set initial value of free parameter to the midrange
        double p = *(ptrM + k);    
        // Start iterations
        for (int Iter = 0; Iter < MaxIter ; Iter++) {
            // Break from iterations if the range of the x values is zero
            if (*(ptrRANGE + k) == 0) {
                break;
            }   
            T = 0;
            v = 0;
            U = 0;
            for (int j = 0; j < m ; j++) {
                double xj = *(ptrX + k * m + j);
                for (int i = 0; i < j ; i++) {
                    double xi = *(ptrX + k * m + i);
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
                octave_stdout << "Warning: Root finding failed to reach the specified tolerance \n";
            }
        }
        // Assign parameter value that optimizes the objective function for the smoothed median
        *(ptrM + k) = p;
    }
    // If applicable, switch dimension
    if (dim > 1) {
        M = M.transpose ();
    }
    return octave_value(M);
} 
