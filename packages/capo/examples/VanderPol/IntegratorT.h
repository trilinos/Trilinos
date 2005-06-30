/***************************************************************************
                               IntegratorT.h
                             -------------------
    written by           : Blake Ashby
    last updated         : Nov 15, 2002
    email                : bmashby@stanford.edu

This is the parent class for StiffIntegratorT (implicit integrator based on
RADAU5) and NonStiffIntegratorT (explicit integrator based on DOPRI5). The
code has been modified from the code written originally in Fortran by

         E. Hairer and G. Wanner
         Universite de Geneve, Dept. de Mathematiques
         Ch-1211 Geneve 24, Switzerland
         E-mail:  ernst.hairer@math.unige.ch
                  gerhard.wanner@math.unige.ch

See the header files StiffIntegratorT.h and NonStiffIntegratorT.h for
more details.

 ***************************************************************************/

#ifndef _INTEGRATOR_T_H_
#define _INTEGRATOR_T_H_

#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include "decsol.h"

void Function(double x, double *y, double *f, double eps);

class IntegratorT
{

public:

	// Constructor	
	IntegratorT(const int nin, double yin[], double xin, const double xendin, double dxin,
		int itolerin = 0, double *rtolerin = 0, double *atolerin = 0, const int ioutin = 1,
		double hin = 0.0, double hmaxin = 0.0, int nmaxin = 0, double uroundin = 0.0,
		double safein = 0.0, double faclin = 0.0, double facrin = 0.0);

	// Still need to implement copy constructor and default constructor

	// Destructor
	virtual ~IntegratorT();

	virtual void Integrate(double eps) = 0;

	// Function that controls the output of the results.
	//Modify this routine according to your needs
	int SolutionOutput();

	// get number of function evaluations
	int NumFunction() const { return nfcn; }
	// get number of attempted steps
	int NumStep() const { return nstep; }
	// get number of accepted steps
	int NumAccept() const { return naccpt; }
	// get number of rejected steps
	int NumReject() const { return nrejct; }

protected:

	// CoreIntegrator
	virtual int CoreIntegrator(double eps) = 0;

	virtual double ContinuousOutput(unsigned i) = 0;

// Member variables

	// dimension of system
	const int n;
	// vector for y values
	double *y;
	// independent variable (usually time)
	double x;
	// final value for independent variable
	const double xend;
	// time step for intermediate output
	double dx;
	// switch for rtol and atol; if itol = 0, rtol and atol are scalars
	// if itol = 1, rtol and atol are vectors
	int itoler;
	// relative error tolerance
	double *rtoler;
	// absolute error tolerance
	double *atoler;
	// variables that indicate whether rtoler and atoler are NULL on input
	// or not--needed for proper memory management in destructor
	bool rtolerNULL;
	bool atolerNULL;
	// routine for dense output at every time step is called if iout = 1
	const int iout;
	// integration step length
	double h;

// Derived variables

	// maximal step size
	double hmax;
	// maximal number of steps
	int nmax;
	// smallest number satisfying 1.0 + uround > 1.0
	double uround;
	// safety factor in step size prediction
	double safe;
	// facl, facr--parameters for step size selection
	double facl;
	double facr;

// Counting variables

	// number of function evaluations (not counting those in numerical
	// Jacobian calculations)
	int nfcn;
	// number of attempted steps
	int nstep;
	// number of accepted steps
	int naccpt;
	// number of rejected steps
	int nrejct;

	// stores past value of x
	double xold;
	// stores past value of h
	double hold;
	// x at discrete points specified by dx interval
	double xd;

};

#endif // _INTEGRATOR_T_H_
