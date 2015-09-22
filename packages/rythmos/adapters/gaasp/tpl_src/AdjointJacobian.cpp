////////////////////////////////////////////////////////////////////////////////
// function: double** adjointJacobian(int const n, double *xin, double *xp, 
//			double *xold, double tnew, double told, int flag)
//
// Written by: Jeff Sandelin
//    Colorado State University
//    Dept. of Mathematics
//    summer, 2004
//
// Computes a finite difference approximation of the Jacobian.
//
// Input: n  - number of equations in the system
//   xin  - point of approximation
//   xp  - predicted values
//   tnew - new time
//   told - old time
//   flag - 0 = forward, 1 = dual
// Output: **x  - two dimension array of approximate Jacobian
////////////////////////////////////////////////////////////////////////////////

#include <limits>
#include "MemoryMgmt.h"
#include "GAdjointSolve.h"
#include <cmath>

namespace GAASP {

using std::sqrt;
using std::numeric_limits;

double** GAdjointSolve::adjointJacobian( int const n, double *xin, double *xp, double *xold, 
		double const tnew, double const told, int const flag )
{

    // Allocate memory for the evaluation point.
    double *xeval = CreateArray<double>( n );	// Evaluation point
    double *xout1 = CreateArray<double>( n );	// Output of function evaluation
    double *xout2 = CreateArray<double>( n );	// Output from non-perturbed function evaluation
    double **jac = CreateArray<double>( n, n );	// Jacobian

    // Compute the approximate Jacobian by column.

    int icol = 0;			// Column counter
    for ( ; icol < n; ++icol )
    {
        for ( int j = 0;j < n; ++j )
            xeval[ j ] = xp[ j ];

        // double dh = 1.01 * pow( machEpsilon, 0.5 );	// Perturbation
        
        double dh = 1.01 * sqrt(machEpsilon);
        dh = double( long( dh * 10e+12 ) / 10e+12 );

        xeval[ icol ] += dh;

        modelPtr->calcDerivs( xp, xout1, tnew );
        modelPtr->calcDerivs( xeval, xout2, tnew );

        // Transfer output to the jacobian array.
        for ( int j = 0;j < n; ++j )
            jac[ j ][ icol ] = ( xout2[ j ] - xout1[ j ] ) / dh;

    }

    DeleteArray ( xout1 );
    DeleteArray ( xout2 );
    DeleteArray ( xeval );
    return jac;
}

} // namespace GAASP
