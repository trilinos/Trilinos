/***************************************************************************
                          StiffIntegratorT.h
                          -------------------
    written by           : Blake Ashby
    last updated         : Nov 15, 2002
    email                : bmashby@stanford.edu

This is the child class to IntegratorT.h. It is modified for C++ from the
code RADAU5 originally written in FORTRAN  (version of July 9, 1996, latest
small correction: January 18, 2002) by:

         E. Hairer and G. Wanner
         Universite de Geneve, Dept. de Mathematiques
         Ch-1211 Geneve 24, Switzerland
         E-mail:  ernst.hairer@math.unige.ch
                  gerhard.wanner@math.unige.ch

RADAU5 is described in the book:

         E. Hairer and G. Wanner, Solving Ordinary Differential
         Equations II. Stiff and Differential-Algebraic Problems.
         Springer Series in Computational Mathematics 14,
         Springer-Verlag 1991, Second Edition 1996.

This code computes the numerical solution of a stiff (or differential algebraic)
system of first order ordinary differential equations

                         M*y'=f(x,y).

The system can be (linearly) implicit (mass-matrix M != I) or explicit (M = I).
the method used is an implicit Runge-Kutta method (RADAU IIA) of order 5 with
step size control and continuous output. (See Section IV.8 of Hairer and Wanner).

USER PROVIDED FUNCTIONS:
----------------------

Function	Function defining the differential equation. It must have the
			following prototype:

				void Function(double x, double *y, double *f);

			where 'x' is the value of the independent variable x at which you want
			'y' integrated. 'y' is the array of the current values of the state
			vector. The array f will be filled with the function result y' = f

Jacobian	Function which computes the partial derivatives of f(x,y) with respect
			to y. (This function is only called if ijac = 1. Supply a dummy
			subroutine in the case ijac = 0).

			For ijac = 1, this function must have the prototype:
			
				void Jacobian(double x, double *y, double **J);

			If (mljac == n) the Jacobian is supposed to be full and the partial
			derivatives are stored in the array J as:

				J[i][j] = partial f[i] / partial y[j]
				
			Else, the Jacobian is taken as banded and the diagonals of the
			Jacobian are stored as rows in J as:

				J[i-j+mujac][j] = partial f[i] / partial y[j]

Mass		Function computing the mass-matrix M.
			If imas = 0, this matrix is assumed to be the identity matrix and
			does not need to be defined. Supply a dummy function in this case.
			If imas = 1, the function has the prototype:
			
				void Mass(double **M);
			
			If (mlmas == n) the mass-matrix is stored as full matrix like:

				M[i][j] = m[i][j]

			Else, the matrix is taken as banded and the diagonals of the mass
			matrix are stored as rows in M as:

				 M[i-j+mumas][j] = m[i][j]

Input Parameters (before call to Integrate())
----------------
	n			Dimension of the system

	x			Initial value of the dependent variable (usually time)

	*y			Initial y values (double y[n])

	xend		Final x-value (xend - x may be positive or negative)

	h			Initial step size guess;
				For stiff equations with initial transient, h = 1.0/(norm of f'),
				usually 1.0e-3 or 1.0e-5, is good. This choice is not very
				important, the step size is quickly adapted.
				(if h = 0.0 on input, the code sets h = 1.0e-6).

	*rtoler		Relative and absolute error tolerances. They can both be scalars or
	*atoler		vectors of length n. (in the scalar case pass the addresses of
				variables where you have placed the tolerance values). If set as NULL
				on input, they are set to 1.0e-7 and itoler is set to 0.

	itoler		Switch for rtoler and atoler:
				itoler = 0:	Both rtoler and atoler are scalars. The code keeps, roughly,
							the local error of y[i] below rtoler*fabs(y[i])+atoler
				itoler = 1: Both rtoler and atoler are vectors. The code keeps the
							local error of y[i] below rtoler[i]*fabs(y[i])+atoler[i].

	ijac		Switch for the computation of the Jacobian:
				ijac = 0:	Jacobian is computed internally by finite differences,
							function "Jacobian" is never called.
				ijac = 1: 	Jacobian is supplied by function "Jacobian".

	mljac		Switch for the banded structure of the Jacobian:
				mljac = n:			Jacobian is a full matrix. The linear algebra is
									done by full-matrix Gauss-elimination.
				0 <= mljac < n: 	mljac is the lower bandwith of Jacobian
									matrix ( >= number of non-zero diagonals below
									the main diagonal).

	mujac		Upper bandwith of Jacobian matrix ( >= number of non-zero diagonals
				above the main diagonal). Need not be defined if mljac = n.

	imas		Gives information on the mass-matrix:
				imas = 0: 	M is supposed to be the identity matrix, Mass is never called.
				imas = 1: 	mass-matrix is supplied in function Mass.

	mlmas		Switch for the banded structure of the mass-matrix:
				mlmas = n: 		The full matrix case. The linear algebra is done by
								full-matrix Gauss-elimination.
				0 <= mlmas < n:	mlmas is the lower bandwith of the matrix ( >= number
								of non-zero diagonals below the main diagonal).
								mlmas is supposed to be <= mljac.

	mumas		Upper bandwith of mass-matrix ( >= number of non-zero diagonals above
				the main diagonal). Need not be defined if mlmas = n. mumas is supposed
				to be <= mujac.

	iout		Switch for calling the function solout:
				iout = 0: function is never called
				iout = 1: function is available for output.

 ----------------------------------------------------------------------

Sophisticated Setting of Parameters
-----------------------------------

Several parameters have a default value (if set to 0) but, to better adapt the
code to your problem, you can specify particular initial values.

	uround		rounding unit, default 1.0e-16.

	nmax		This is the maximal number of allowed steps. The default value
				is 100000.

	nit			The maximum number of Newton iterations for the solution of the
				implicit system in each step. The default value is 7.

	startn		If startn == 0 the extrapolated collocation solution is taken as
				starting value for Newton's method. If startn != 0 zero starting
				values are used. The latter is recommended if Newton's method has
				difficulties with convergence. (This is the case when nstep is larger
				than naccpt + nrejct; see output parameters). Default is startn = 0.

The following 3 parameters are important for differential-algebraic systems of
index > 1. The function-subroutine should be written such that the index 1, 2, 3
variables appear in this order. In estimating the error the index 2 variables are
multiplied by h, the index 3 variables by h^2.

	nind1		Dimension of the index 1 variables (must be > 0). For ODE's this
				equals the dimension of the system. Default nind1 = n.

	nind2		Dimension of the index 2 variables. Default nind2 = 0.

	nind3		Dimension of the index 3 variables. Default nind3 = 0.

	npred		Switch for step size strategy
					If npred = 1  mod. predictive controller (Gustafsson)
					If npred = 2  classical step size control
				The default value (for npred = 0) is npred = 1. The choice
				npred = 1 seems to produce safer results. For simple problems,
				the choice npred = 2 produces often slightly faster runs.

If the differential system has the special structure that

			y[i]' = y[i+m2]   for  i = 1, ... , m1,

with m1 a multiple of m2, a substantial gain in computertime can be achieved by
setting the parameters m1 and m2. e.g., for second order systems
p' = v, v' = g(p,v), where p and v are vectors of dimension n/2, one has to put
m1 = m2 = n/2.

For m1 > 0 some of the input parameters have different meanings:

	Jacobian:	Only the elements of the non-trivial part of the Jacobian have to
				be stored.
				If (mljac == n-m1) the Jacobian is supposed to be full:
					fjac[i][j] = partial f[i+m1] / partial y[j]
						for i = 0, n-m1-1 and j = 0, n-1
				Else, the Jacobian is banded ( m1 = m2 * mm )
					fjac[i-j+mujac][j+k*m2] = partial f[i+m1] / partial y[j+k*m2]
						for i = 0, mljac + mujac and j = 0, m2-1 and k = 0, mm.

	mljac: 		mljac = n - m1:
					If the non-trivial part of the Jacobian is full
                0 <= mljac < n - m1:
					If the (mm+1) submatrices (for k = 0, mm)
					partial f[i+m1] / partial y[j+k*m2],  i, j = 0, m2-1
					are banded, mljac is the maximal lower bandwidth of these
					mm+1 submatrices

	mujac:		Maximal upper bandwidth of these mm + 1 submatrices. Need not be
				defined if mljac = n - m1

	Mass: 		If imas = 0 this matrix is assumed to be the identity and need
				not be defined. Supply a dummy subroutine in this case. It is
				assumed that only the elements of right lower block of dimension
				n - m1 differ from that of the identity matrix.
				If (mlmas == n - m1) this submatrix is supposed to be full
					fmas[i][j] = m[i+m1][j+m1]
						for i = 0, n-m1-1 and j = 0, n-m1-1
				Else, the mass matrix is banded
					fmas[i-j+mumas][j] = m[i+m1][j+m1]

	mlmas: 		mlmas = n - m1: 		if the non-trivial part of m is full
				0 <= mlmas < n - m1: 	lower bandwidth of the mass matrix

	mumas: 		Upper bandwidth of the mass matrix. Need not be defined if mlmas = n-m1

	m1			Default m1 = 0.

	m2			Default m2 = m1.

	safe		The safety factor in step size prediction, default 0.9.

	thet		Decides whether the Jacobian should be recomputed. Increase thet,
				to 0.1 say, when Jacobian evaluations are costly. for small systems
				thet should be smaller (0.001, say). Negative thet forces the
				code to compute the Jacobian after every accepted step. Default 0.001.

	fnewt		Stopping criterion for Newton's method, usually chosen < 1. Smaller
				values of fnewt make the code slower, but safer.
				Default min(0.03, sqrt(rtoler))

	quot1, quot2
				If quot1 < hnew/hold < quot2, then the step size is not changed.
				This saves, together with a large thet, lu-decompositions and
				computing time for large systems. for small systems one may have
				quot1 = 1.0, quot2 = 1.2, for large full systems quot1 = 0.99,
				quot2 = 2.0 might be good. Defaults quot1 = 1.0, quot2 = 1.2.

	hmax		Maximal step size, default xend - x.

	facl, facr	Parameters for step size selection the new step size is chosen subject
				to the restriction 1/fac1 <= hnew/hold <= 1/facr.
				Default values: fac1 = 5.0, facr = 1.0/8.0

	hess		If hess != 0, the code transforms the Jacobian matrix to
				Hessenberg form. This is particularly advantageous for large
				systems with full Jacobian. It does not work for banded Jacobian
				(mljac < n) and not for implicit systems (imas = 1).

-----------------------------------------------------------------------

Output parameters (after call to Integrate())
-----------------
	x		x-value for which the solution has been computed
			(after successful return x = xend).

	*y		Numerical solution at x  (y[n])

	h		Predicted step size of the last accepted step

 ***************************************************************************/

#ifndef _STIFF_INTEGRATOR_T_H_
#define _STIFF_INTEGRATOR_T_H_

#include "IntegratorT.h"

void Jacobian(double x, double *y, double **J, double eps);
void Mass(double **M);

class StiffIntegratorT: public IntegratorT
{

public:
	//Constructors

	StiffIntegratorT(const int nin, double yin[], double xin, double xendin,
		double dxin, int itolerin, double *rtolerin, double *atolerin,
		const int ioutin, double hin, double hmaxin, int nmaxin, double uroundin,
		double safein, double faclin, double facrin, const int ijacin,
		int mljacin, int mujacin, const int imasin, int mlmasin,
		int mumasin, int nitin = 0, bool startnin = false, int nind1in = 0,
		int nind2in = 0, int nind3in = 0, int npredin = 0, int m1in = 0, int m2in = 0,
		bool hessin = false, double fnewtin = 0, double quot1in = 0, double quot2in = 0,
		double thetin = 0, double eps=1.0e-6);

	StiffIntegratorT(int nin, double yin[], double xin, double xendin, double dxin,
		const int ijacin, int mljacin, int mujacin, const int imasin,
		int mlmasin, int mumasin);

	// Still need to implement copy constructor and default constructor
	
	~StiffIntegratorT();

	void Integrate(double eps);
	
	// get number of Jacobian evaluations
	int NumJacobian() const { return njac; }
	
	// get number of lu-decompositions of both matrices
	int NumDecomp() const { return ndec; }
	
	// get number of forward-backward substitutions, of both systems;
	int NumSol() const { return nsol; }

private:

// member routines

	virtual int CoreIntegrator(double eps);

	virtual double ContinuousOutput(unsigned i);
	
	void ComputeJacobian(double eps);

// Linear Algebra routines:
	
	int DecompReal();

	int DecompComplex();

	int LinearSolve();

	int ErrorEstimate(double eps);

// member variables

	// compute the Jacobian analytically (ijac = 1) or not (ijac = 1)
	const int ijac;
	// number of non-zero diagonals below main diagonal of Jacobian matrix
	//(if mljac = dimx, matrix is full)
	int mljac;
	// number of non-zero diagonals above main diagonal of Jacobian matrix
	int mujac;
	// differential equation is in explicit form (imas = 0) or not (imas = 1)
	const int imas;
	// number of non-zero diagonals below main diagonal of mass matrix
	//(if mlmas = dimx, matrix is full)
	int mlmas;
	// number of non-zero diagonals above main diagonal of mass matrix
	int mumas;
	// maximal number of Newton iterations
	int nit;
	// switch for starting values of Newton iterations
	bool startn;
	// parameters for differential-algebraic components
	int nind1;
	int nind2;
	int nind3;
	// step size control
	int npred;
	bool pred;
	// parameters for second order equations
	int m1;
	int m2;
	int nm1;
	// Convert Jacobian to Hessenberg form?
	int hess;
	// stopping criterion for Newton's method, usually chosen < 1
	double fnewt;
	// quot1 and quot2--if quot1 < hnew/hold < quot2, step size = const
	double quot1;
	double quot2;
	// decides whether the Jacobian should be recomputed
	double thet;
	// implicit, banded or not?
	bool implct;
	bool jband;
	// row-dimensions of the 2-D arrays -- fjac, e1, e2, and fmas
	int ldjac;
	int lde1;
	int ldmas;
	// job number (used in linear algebra routines for type of problem--banded or not, etc)
	int ijob;
	
	// number of Jacobian evaluations
	int njac;
	// number of lu-decompositions of both matrices
	int ndec;
	// number of forward-backward substitutions, of both systems;
	// the nstep forward-backward substitutions, needed for step
	//size selection, are not counted
	int nsol;

	// constants that are used in linear algebra routines
	int mle;
	int mue;
	int mbjac;
	int mbb;
	int mdiag;
	int mdiff;
	int mbdiag;

	double fac1;
	double alphn;
	double betan;

	// variables used in program
	double err;
	bool caljac;
	bool calhes;
	bool first;
	bool reject;

	// arrays used in program
	double *z1;
	double *z2;
	double *z3;
	double *y0;
	double *scal;
	double *f1;
	double *f2;
	double *f3;
	double *cont;
	int *ip1;
	int *ip2;
	int *iphes;
	double **e1;
	double **e2r;
	double **e2i;

	// Jacobian matrix
	double **fjac;
	// mass matrix
	double **fmas;

};

#endif // _STIFF_INTEGRATOR_T_H_
