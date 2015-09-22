///////////////////////////////////////////////////////////////////////////////
//	function: void spline(double *x, double *y, int n, double yp1, double ypn,
//							double *y2)
//	Compute the second derivative values for the interpolating function.
//	Input:	*x		-	tabulated x values
//			*y		-	tabulated y=f(x) values
//			n		-	size of the arrays
//			yp1		-	first derivative at point 1
//			ypn		-	first derivative at point n
//	Output:	*y2		-	array of second derivatives
//
//	Slightly modified code from Numerical Recipes in C, Section 3.3
///////////////////////////////////////////////////////////////////////////////

#ifndef INC_Spline_h
#define INC_Spline_h

#include "MemoryMgmt.h"
#include <cmath>

namespace GAASP {

//void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
void spline ( double *x, double *y, int n, double * const y2 );

//*****************************************************************************
//	function: int splint(double *xa, double *ya, double *y2, int n,
//							double x, double *y)
//	Compute the second derivative values for the interpolating function.
//	Input:	*xa		-	tabulated x values
//		*ya		-	tabulated y=f(x) values
//		*y2		-	output from spline()
//		n		-	size of the arrays
//		x		-	input x value
//	Output:	*y		-	interpolated y value
//	Returns:		0 = success
//				1 = error
//
//	Slightly modified code from Numerical Recipes in C, Section 3.3
//*****************************************************************************
int splint (
	double  const * const xa,
	double  const * const ya,
	double  const * const y2,
	int  const n,
	double  const x,
	double * const y );

} // namespace GAASP

#endif // INC_Spline_h

