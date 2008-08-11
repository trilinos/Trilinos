#include <stdlib.h>
//*****************************************************************************
//	function: matrixSolve(double**, double *)	
//	Written by:	Jeff Sandelin
//	Solves Ax = b for x.
//
//	Input:	sysSize	-	number of equations in the system
//		Ma	-	2-D array for matrix A
//		lv	-	load vector b
//	Output:	*bv	-	values for the solution x
//*****************************************************************************

//	Include files for the solve functions.

#include "dense.h"
#include "llnltyps.h"
#include "nvector.h"
#include "llnlmath.h"

namespace GAASP {

using namespace CVODE;

//	Defines for the linear system solver.

#define Ith(v,i)	N_VIth(v,i-1)	
#define IJth(A,i,j)	DENSE_ELEM(A,i-1,j-1)


double *matrixSolve(int const sysSize, double **Ma, double *lv)
{

	DenseMat M;		// Matrix to solve
	N_Vector b;		// RHS
	int *pivots;		// Pivots
	int i;			// Loop counter
	int row, col;		// Loop counters
	int ier;		// Error code
	double *bv;		// b vector for return
	double junk;

	// Initialize the matrix to be solved.

	M = DenseAllocMat(sysSize);
	pivots = DenseAllocPiv(sysSize);
	b = N_VNew(sysSize, NULL);

	DenseZero(M);		// Initialize the matrix with zeros

	bv = new double[sysSize];	// Allocate memory for the solution

	// Populate the matrix M with the input values.

	for (row = 1; row <= sysSize; ++row )
	{
		for (col = 1; col <= sysSize; ++col )
		{
			junk = Ma[row-1][col-1];
			IJth(M,row,col) = junk;
		}
	}

	// Populate the load vector with the function evaluation.

	for (row = 1; row <= sysSize; ++row )
	{
		Ith(b,row) = lv[row-1]; 
	}

	// Do the LU factorization of M.

	ier = DenseFactor(M, pivots);

	// Matrix solve.

	DenseBacksolve(M, pivots, b);

	// Free allocated memory.

	DenseFreeMat(M);
	DenseFreePiv(pivots);

	// Transfer results back to a c array and return.

	for (i = 0; i < sysSize; ++i )
	{
		bv[i] = Ith(b, i+1);
	}
	N_VFree (b);
	return bv;
}

} // namespace GAASP

