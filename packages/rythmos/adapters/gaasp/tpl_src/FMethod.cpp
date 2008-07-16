#include "GForwardSolve.h"

namespace GAASP {

void GForwardSolve::fMethod ( double *xold, double *xin, double *x, 
	       double const tnew, double const told, int const n )
{
	int i, hn;
	double dt;
	double *xtmp1, *xtmp2;
	double tau1, tau2;
	double alpha, beta;

	dt = tnew - told;

	switch ( MethodFlag )
	{
	case Method_DG0:	// dG(0)

		xtmp1 = new double[n];
		modelPtr->calcDerivs(xin, xtmp1, tnew);
		for( i = 0; i < n; ++i )
		{
			x[i] = xin[i] - xtmp1[i]*dt - xold[i];
		}
		delete xtmp1;
		break;

	case Method_DG1:	// dG(1)

		xtmp1 = new double[n];
		xtmp2 = new double[n];
		hn = n/2;
		tau1 = told + dt/3.0;		
		tau2 = tnew;
		alpha = 5.0/12.0;
		beta = 1.0/12.0;

//		(*f)(&xin[0], xtmp1, tau1);		// Evaluate f at z1,tau1
//		(*f)(&xin[hn], xtmp2, tau2);	// Evaluate f at z2,tau2
		modelPtr->calcDerivs( &xin[0], xtmp1, tau1);	// Evaluate f at z1,tau1
		modelPtr->calcDerivs( &xin[hn], xtmp2, tau2);	// Evaluate f at z2,tau2

		for( i = 0; i < hn; ++i )
		{
			x[i] = xin[i] - xold[i] - alpha*dt*xtmp1[i] + beta*dt*xtmp2[i];
		}

		for( i = hn; i < n; ++i )
		{
			x[i] = xin[i] - xold[i-hn] - 0.75*dt*xtmp1[i-hn] - 0.25*dt*xtmp2[i-hn];
		}

		delete [] xtmp1;
		delete [] xtmp2;
		break;

	}

	return;
}

} // namespace GAASP

