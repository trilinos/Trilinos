//*****************************************************************************
//	function: newtonStep(double**, double *)
//	Written by:	Jeff Sandelin
//	Solves F'(x)H = -F(x) for H, the Newton step.
//
//	Input:	N		-	number of equations in the system
//			fPrime	-	finite difference approximation of the Jacobian
//			y		-	vector of f evaluated at x
//	Output:	*H		-	vector of the Newton step values
//*****************************************************************************

//	Include files for the solve functions.

#include <stdlib.h>
#include "dense.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"
#include "MemoryMgmt.h"

namespace GAASP {

using namespace CVODE;

//	Defines for the linear system solver.

#define Ith(v,i)	N_VIth(v,i-1)
#define IJth(A,i,j)	DENSE_ELEM(A,i,j)

double *newtonStep( int const N, double **fPrime, double *y )
{
    DenseMat M = DenseAllocMat( N );		// Matrix to solve
    DenseMat J = DenseAllocMat( N );		// Saved Jacobian
    int *pivots = DenseAllocPiv( N );		// Pivots
    N_Vector b = N_VNew( N, NULL );		// RHS

    int i;			// Loop counter
    int row, col;		// Loop counters
    int ier;			// Error code
    double *bv;			// b vector for return
    double junk;

    // gamma = 1.0;
    DenseZero( M );		// Initialize the matrix with zeros
    DenseCopy( M, J );	// Copy the initialized matrix into the Jacobian matrix

    bv = CreateArray<double>(N);

    // Populate the Jacobian with the values from the finite difference approximation.
    // Changed to run 0 to < N so don't need to subtract indexing - RM!!!

    for ( int row = 0; row < N; ++row )
    {
		for ( col = 0; col < N; ++col )
		{
			junk = fPrime[ row ][ col ];
			IJth( M, row, col ) = junk;
		}
    }

    // Populate the load vector with the function evaluation.

    for ( row = 1; row <= N; ++row )
	Ith( b, row ) = y[ row - 1 ];

    // Scale and add I to get M = I - gamma*J
    // It is unclear at this point whether these function calls are
    // needed or not.
    //	DenseScale(-gamma, M);
    //	DenseAddI(M);

    // Do the LU factorization of M.
    ier = DenseFactor( M, pivots );

    // Matrix solve.
    DenseBacksolve( M, pivots, b );

    // Free allocated memory.
    DenseFreeMat( M );
    DenseFreeMat( J );
    DenseFreePiv( pivots );

    // Transfer results back to a c array.
    for ( i = 0; i < N; ++i )
	bv[ i ] = Ith( b, i + 1 );
    N_VFree ( b );

    return bv;
}

} // namespace GAASP

