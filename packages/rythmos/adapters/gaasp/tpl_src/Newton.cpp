#include "GForwardSolve.h"
#include "MemoryMgmt.h"

namespace GAASP {
double *newtonStep( int const, double **, double * );
double l2Norm( double * const, int const );

bool GForwardSolve::newton ( double *yc, double *yp, double *yout, double *yold, 
	      double const tnew, double const told, int const n, int const currentStep )
{
    bool converged = false;	// Return value
    double const nTOL = 1e-10;	// Acceptance tolerance
    int const MAXITER = 30;	// Maximum number of iterations
    double r = 0;

    // Allocate memory.  Note that for the higher order method (METHOD=1) the
    // system doubles in dimension to account for the intermediate nodes.

    int tmpSizeMultiplier = 1;

    // Change if higher order
    if (MethodFlag == Method_DG1) tmpSizeMultiplier = 2;

    int const sysSize = n * tmpSizeMultiplier;  // Size of the system

    double *b = CreateArray<double>(sysSize);
    double *y = CreateArray<double>(sysSize);

    // Evaluate the Jacobian at current step.

    // Jacobian
    double **J = forwardJacobian( sysSize, yc, yp, yold, tnew, told, 0 );

    // Determine the Newton step.

    for ( int j = 0;j < MAXITER; ++j )
    {
        // Compute F at current step and change the sign.

        fMethod( yold, yp, b, tnew, told, sysSize );
		r = l2Norm( b, sysSize );
		double rchk = r;

        for ( int i = 0;i < sysSize; ++i )
            b[ i ] *= -1.0;

        // Solve for the Newton direction.

		double *H = newtonStep( sysSize, J, b );     // Newton step

		double lambda = 1.0;			// Newton step length
		double oldlambda = 2.0;
		while ( r >= rchk )
        {
 			if ( oldlambda - lambda < 10e-10 )
            {   // Clean up
                DeleteArray (H);
                DeleteArray (J, sysSize);	// delete Jacobian matrix
                DeleteArray (b);
                DeleteArray (y);
				return 1;
            }

		    // Compute the new x values.
		    for ( int i = 0;i < sysSize; ++i )
			y[ i ] = yp[ i ] + lambda * H[ i ];

		    // Compute the function evaluation at the new point
		    // and determine the residual.
			fMethod( yold, y, b, tnew, told, sysSize );
            r = l2Norm( b, sysSize );
            if ( r < nTOL )
                break;
            oldlambda = lambda;
            lambda *= 0.5;
        }

        DeleteArray (H);
        if ( r < nTOL )
        {
            converged = true;
            break;
        }

        for ( int i = 0; i < sysSize; ++i )
            yp[ i ] = y[ i ];

    }

    DeleteArray (J, sysSize);		// delete Jacobian matrix

    // Transfer to output.
    switch ( MethodFlag )
    {
	case Method_DG0:        // dG(0)

	    for(int i=0; i<n; ++i)
			yout[i] = y[i];

	    break;

	case Method_DG1:        // dG(1)

	    for(int i=0; i<n; ++i)
			derivatives[i][currentStep] = (y[n+i] - y[i]) / ((tnew - told)*2.0/3.0);

		for(int i=0; i<n; ++i)
			yout[i] = y[n+i];

	    break;
    }

    // Clean up.
    DeleteArray (b);
    DeleteArray (y);

    return converged;
}

} // namespace GAASP

