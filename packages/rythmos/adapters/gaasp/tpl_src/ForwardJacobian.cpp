#include "GForwardSolve.h"
#include "MemoryMgmt.h"
#include <cmath>

namespace GAASP {

using std::sqrt;

double** GForwardSolve::forwardJacobian( int const n, double *xin, double *xp, double *xold, 
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

        fMethod( xold, xp, xout1, tnew, told, n );
        fMethod( xold, xeval, xout2, tnew, told, n );

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

