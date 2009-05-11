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

#include "MemoryMgmt.h"
#include "Spline.h"
#include <cmath>

namespace GAASP {

using std::pow;

//void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
void spline ( double *x, double *y, int n, double * const y2 )
{
    int i, k;		// Loop counters
    double p, qn, sig, un, *u = 0, yp1, ypn;

    //	Allocate memory for the vector.
    u = CreateArray<double>(n);

    //	Compute first derivatives at point 1 and point n.
    yp1 = ( y[1] - y[0] ) / ( x[1] - x[0] );
    ypn = ( y[n-1] - y[n-2] ) / ( x[n-1] - x[n-2] );

    //	Spline computation.
    if ( yp1 > 0.99e30 )
	y2[0] = u[0] = 0.0;
    else
    {
	y2[0] = -0.5;
	u[0] = ( 3.0 / ( x[1] - x[0] ) ) * ( ( y[1] - y[0] ) / ( x[1] - x[0] ) - yp1 );
    }

    for ( i = 1;i <= n - 2; ++i )
    {
	sig = ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] );
	p = sig * y2[i-1] + 2.0;
	y2[i] = ( sig - 1.0 ) / p;
	u[i] = ( y[i+1] - y[i] ) / ( x[i+1] - x[i] ) - ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
	u[i] = ( 6.0 * u[i] / ( x[i+1] - x[i-1] ) - sig * u[i-1] ) / p;
    }

    if ( ypn > 0.99e30 )
	qn = un = 0.0;
    else
    {
	qn = 0.5;
	un = ( 3.0 / ( x[n-1] - x[n-2] ) ) * ( ypn - ( y[n-1] - y[n-2] ) / ( x[n-1] - x[n-2] ) );
    }

    y2[n-1] = ( un - qn * u[n-2] ) / ( qn * y2[n-2] + 1.0 );

    for ( k = n - 2;k >= 0;k-- )
    {
	y2[k] = y2[k] * y2[k+1] + u[k];
    }

    //	Free the memory.
    DeleteArray(u);

}

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
	double * const y )
{
    int klo, khi, k;
    double h, b, a;

    klo = 1;
    khi = n;

    while ( khi - klo > 1 )
    {
	k = ( khi + klo ) >> 1;
	if ( xa[k] > x )
	    khi = k;
	else
	    klo = k;
    }

    h = xa[khi] - xa[klo];
    if ( h == 0.0 ) return 1;
    a = ( xa[khi] - x ) / h;
    b = ( x - xa[klo] ) / h;
    *y = a * ya[klo] + b * ya[khi] + ( ( pow ( a, 3 ) - a ) * y2[klo] + ( pow ( b, 3 ) - b ) * y2[khi] ) * ( h * h ) / 6.0;

    return 0;
}

} // namespace GAASP

