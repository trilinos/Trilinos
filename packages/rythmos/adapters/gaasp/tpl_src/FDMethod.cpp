///////////////////////////////////////////////////////////////////////////////
//	function: void fdMethod(double *yold, double *xin, double *x, double tnew,
//								double told, double T, int n)
//
//	Written by:	Jeff Sandelin
//				Colorado State University
//				Dept. of Mathematics
//				summer, 2004
//
//	Evaluates the integral for the dual problem using requested order.
//
//	Input:	yold	-	function values at current node
//			xin		-	new point for function evaluation
//			tnew	-	time at new node
//			told	-	time at current node
//			n		-	size of the system
//	Output:	*sln	-	array of new y values
///////////////////////////////////////////////////////////////////////////////

#include "GAdjointSolve.h"
#include "Spline.h"

namespace GAASP {

void transpose( double **, int const, int const, double ** );
double *matrixSolve( int const, double **, double * );

double **GAdjointSolve::fdMethod( double *yold, double *xin, double tnew,
                   double told, double T, int dsteps, int sysDim )
{

    int i, k, j;
    int err;
    double **J;
    double *jy1, *jy2;					// Matrix multiplication results
    double yspln;

    //	Constants used in two stage order four IRK

    double C1 = 0.2113248654051871;	// .5 - sqrt(3)/6
    double C2 = 0.7886751345948130;	// .5 + sqrt(3)/6
    double A11 = 0.25;
    double A12 = -0.0386751345948129;	// .25 - sqrt(3)/6
    double A21 = 0.5386751345948129;	// .25 + sqrt(3)/6
    double A22 = 0.25;

    double dt = told - tnew;		// Change in t

    double tau1;			// Interior time node
    double tau2;			// Interior time node

    //	Set matrix size.

    int const s = sysDim * ( MethodFlag + 1 );	// Dimension of matrix

    //	Allocate memory.

    double **J1;			// Jacobian
    double **J2;
    double **Jt1;			// Transpose of the Jacobian
    double **Jt2;
    Jt1 = new double * [ s ];
    for ( i = 0;i < s; ++i )
    {
        Jt1[ i ] = new double[ s ];
    }
    Jt2 = new double * [ s ];
    for ( i = 0;i < s; ++i )
    {
        Jt2[ i ] = new double[ s ];
    }
    J = new double * [ s ];
    for ( i = 0;i < s; ++i )
    {
        J[ i ] = new double[ s ];
    }

    double *yn;
    double *ynplus1;
    double *lv;							// Load vector
    yn = new double[ sysDim ];
    ynplus1 = new double[ sysDim ];
    lv = new double[ sysDim ];
    double **sln1;				// Solution METHOD = 1
    //	Note: use sln1 to return values to calling function in either case.
    //	This just means we have to transfer sln0 to sln1 before return.
    sln1 = new double * [ sysDim ];
    j = MethodFlag == 0 ? 1 : 3;
    for ( i = 0; i < sysDim; ++i )
    {
        sln1[ i ] = new double[ j ];
    }

    //	Do the requested order.

    switch ( MethodFlag )
    {
        case 0:     		// cG(1)
            {

                //	Get the spline value for y at s(sysDim).

                if ( told == tNodes[ 0 ] )
                {
                    for ( i = 0; i < sysDim; ++i )
                    {
                        yn[ i ] = spln[ i ][ dsteps ];
                    }
                }
                else
                {
                    for ( i = 0; i < sysDim; ++i )
                    {
                        err = splint( tNodes, fsoln[ i ], y2[ i ], dsteps, told, &yspln );
                        yn[ i ] = yspln;
                    }
                }

                //	Get the spline value for y at s(n+1).

                for ( i = 0; i < sysDim; ++i )
                {
                    err = splint( tNodes, fsoln[ i ], y2[ i ], dsteps, tnew, &yspln );
                    ynplus1[ i ] = yspln;
                }

                //*** Note: "told" in the 6th argument of both calls to fdjac
                //***       is simply a placeholder and not used for the dual
                //***       solve.  Also, the 0 in the last place makes no
                //***		difference.

                //	Compute the Jacobian at t(sysDim+1) using Y(sysDim+1) and transpose.

                J1 = adjointJacobian( sysDim, xin, ynplus1, yn, tnew, told, 1 );
                transpose( J1, sysDim, sysDim, Jt1 );

                //	Compute the Jacobian at t(sysDim) using Y(sysDim)

                J2 = adjointJacobian( sysDim, xin, yn, yn, told, told, 1 );
                transpose( J2, sysDim, sysDim, Jt2 );

                //	Build the matrix.
                //	=> [I-f'(sysDim+1)^T*dt/2]

                for ( i = 0; i < sysDim; ++i )
                {
                    for ( k = 0; k < sysDim; ++k )
                    {
                        J[ i ][ k ] = 0.0;
                        J[ i ][ k ] -= Jt1[ i ][ k ] * dt / 2.0;
                    }
                    J[ i ][ i ] += 1.0;
                }

                //	Evaluate psi at told and tnew.
		//	!!! - Becky, dsteps and time are not used in computation of psi

                double *pt1;		// psi at old time or tau1
                double *pt2;		// psi at new time or tau2
				pt1 = psi();
				pt2 = psi();

                //	Compute the load vector.
                //	[I + f'(sysDim)^T*dt/2]*phi(sysDim) + (psi(sysDim+1)+psi(sysDim))*dt/2

                for ( i = 0; i < sysDim; ++i )
                {
                    lv[ i ] = 0.0;
                    for ( k = 0; k < sysDim; ++k )
                    {
                        lv[ i ] += Jt2[ i ][ k ] * yold[ k ] * dt / 2.0;
                    }
                    lv[ i ] += yold[ i ] + ( pt2[ i ] + pt1[ i ] ) * dt / 2.0;
                }

                //	Pass matrix and load vector to the solver.

                double *sln0 = matrixSolve( sysDim, J, lv );	// Solution METHOD = 0

                //	Return results.

                for ( i = 0; i < sysDim; ++i )
                {
                    sln1[ i ][ 0 ] = sln0[ i ];
                }

                delete [] pt1;
                delete [] pt2;
                delete [] sln0;
                break;
            }

        case 1:     		// cG(2)
            {

                //	This is a two stage, order four implicit Runge Kutta method.
                //	First step is to build the matrix and solve for the
                //	intermediate y-tilde values.  These are then plugged into
                //	the equation to solve for y(sysDim+1).

                //Compute the interior nodes.

                tau1 = told - C1 * dt;			// Interior time node
                tau2 = told - C2 * dt;			// Interior time node

                //	Initialize the matrix.

                for ( i = 0;i < s; ++i )
                {
                    for ( k = 0;k < s;k++ )
                    {
                        J[ i ][ k ] = 0.0;
                    }
                }

                //	Get the spline value for y at tau1.

                for ( i = 0; i < sysDim; ++i )
                {
                    err = splint( tNodes, fsoln[ i ], y2[ i ], dsteps, tau1, &yspln );
                    yn[ i ] = yspln;
                }

                //	Get the spline value for y at tau2.

                for ( i = 0; i <  sysDim; ++i )
                {
                    err = splint( tNodes, fsoln[ i ], y2[ i ], dsteps, tau2, &yspln );
                    ynplus1[ i ] = yspln;
                }

                //	Compute the Jacobian at tau1 using Y(tau1) and transpose.

                J1 = adjointJacobian( sysDim, xin, yn, yn, tau1, told, 1 );
                transpose( J1, sysDim, sysDim, Jt1 );

                //	Compute the Jacobian at tau2 using Y(tau2)

                J2 = adjointJacobian( sysDim, xin, ynplus1, yn, tau2, told, 1 );
                transpose( J2, sysDim, sysDim, Jt2 );

                //	Build the matrix.

                for ( i = 0; i < sysDim; ++i )     	// Upper left block
                {
                    for ( k = 0;k < sysDim; ++k )
                    {
                        J[ i ][ k ] -= Jt1[ i ][ k ] * dt * A11;
                    }
                    J[ i ][ i ] += 1.0;
                }

                for ( i = 0; i < sysDim; ++i )     	// Upper right block
                {
                    for ( k = sysDim; k < s; ++k )
                    {
                        J[ i ][ k ] -= Jt2[ i ][ k - sysDim ] * dt * A12;
                    }
                }

                for ( i = sysDim; i < s; ++i )     	// Lower left block
                {
                    for ( k = 0; k < sysDim; ++k )
                    {
                        J[ i ][ k ] -= Jt1[ i - sysDim ][ k ] * dt * A21;
                    }
                }

                for ( i = sysDim; i < s; ++i )     	// Lower right block
                {
                    for ( k = sysDim; k < s; ++k )
                    {
                        J[ i ][ k ] -= Jt2[ i - sysDim ][ k - sysDim ] * dt * A22;
                    }
                    J[ i ][ i ] += 1.0;
                }

                //	Evaluate psi at tau1 and tau2.
				//	!!! - Becky, dsteps and tau# are not currently used in computation of psi

                double *pt1;		// psi at old time or tau1
                double *pt2;		// psi at new time or tau2
				pt1 = psi();
                pt2 = psi();

                //	Compute the load vector.

                delete [] lv;
                lv = new double[ s ];
                for ( i = 0; i < sysDim; ++i )
                {
                    lv[ i ] = yold[ i ] + dt * A11 * pt1[ i ] + dt * A12 * pt2[ i ];
                    lv[ i + sysDim ] = yold[ i ] + dt * A21 * pt1[ i ] + dt * A22 * pt2[ i ];
                }

                //	Solve for ytilde.
                double *ytilde = matrixSolve( s, J, lv );

		// !!! Tom - dimensions are s and sysDim, ytilde, yt1, 2

                //	Break the ytilde solution into the component vectors.
                // Y-tilde vectors for function input
                double *yt1 = new double[ sysDim ];
                double *yt2 = new double[ sysDim ];

                for ( i = 0; i < sysDim; ++i )
                {
                    yt1[ i ] = ytilde[ i ];
                    yt2[ i ] = ytilde[ sysDim + i ];
                }

                //	Compute new Y values.
                //	This computation requires Jt1 * yt1 and Jt2 * yt2.
                //	Matrix multiplications first.

//				double *jy1, *jy2;	// Matrix multiplication results
				jy1 = new double[ sysDim ];
                jy2 = new double[ sysDim ];

                for ( i = 0; i < sysDim; ++i )
                {
                    jy1[ i ] = jy2[ i ] = 0.0;
                    for ( k = 0;k < sysDim; ++k )
                    {
                        jy1[ i ] += ( yt1[ k ] * Jt1[ i ][ k ] );
                        jy2[ i ] += ( yt2[ k ] * Jt2[ i ][ k ] );
                    }
                }

                double *y = new double[ sysDim ];

                //	Y value computation.

                for ( i = 0; i < sysDim; ++i )
                {
                    sln1[ i ][ 1 ] = yt1[ i ];
                    sln1[ i ][ 2 ] = yt2[ i ];
                }
                for ( i = 0; i < sysDim; ++i )
                {
                    //	sln1[i][0] = yold[i] + (dt*jy1[i]/2.0) + (dt*jy2[i]/2.0);
                    sln1[ i ][ 0 ] = yold[ i ] + ( dt / 2.0 ) * ( jy1[ i ] + pt1[ i ] ) + ( dt / 2.0 ) * ( jy2[ i ] + pt2[ i ] );
                }

                //	Clean up.

                delete [] ytilde;  ytilde = 0;
                delete [] yt1;     yt1 = 0;
                delete [] yt2;     yt2 = 0;
                delete [] y;       y = 0;
                delete [] pt1;     pt1 = 0;
				delete [] pt2;     pt2 = 0;
				delete [] jy1;     jy1 = 0;
				delete [] jy2;     jy2 = 0;
                break;
            }

        default:     	// Invalid method specified.

            break;

    }

    //	Clean up.

    for ( i = 0; i < sysDim; ++i )
    {
        delete [] J1[ i ];
        delete [] J2[ i ];
    }
    delete [] J1;
    delete [] J2;

    for ( i = 0;i < s; ++i )
    {
        delete [] J[ i ];
        delete [] Jt1[ i ];
        delete [] Jt2[ i ];
    }
    delete [] J;
    delete [] Jt1;
    delete [] Jt2;

    delete [] yn;
    delete [] ynplus1;
    delete [] lv;

    //	Done.

    return sln1;

}

void transpose ( double **inArray, int const rows, int const cols, double **outArray)
{
	for ( int i = 0; i < rows; ++i )
	{
		for ( int j = 0; j<cols; ++j )
		{
			outArray[j][i] = inArray[i][j];
		}
	}
}

} // namespace GAASP

