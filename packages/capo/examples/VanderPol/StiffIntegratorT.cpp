/***************************************************************************
                            StiffIntegratorT.cpp
                             -------------------
    written by:          : Blake Ashby
    last modified        : Nov 15, 2002
    email                : bmashby@stanford.edu
 ***************************************************************************/

#include "StiffIntegratorT.h"

// Constructor
StiffIntegratorT::StiffIntegratorT(const int nin, double yin[], double xin, double xendin,
	double dxin, int itolerin, double *rtolerin, double *atolerin, const int ioutin,
	double hin, double hmaxin, int nmaxin, double uroundin, double safein,
	double faclin, double facrin, const int ijacin, int mljacin, int mujacin,
	const int imasin, int mlmasin, int mumasin, int nitin, bool startnin,
	int nind1in, int nind2in,
	int nind3in, int npredin, int m1in, int m2in, bool hessin, double fnewtin,
	double quot1in, double quot2in, double thetin, double eps) :
	IntegratorT(nin, yin, xin, xendin, dxin, itolerin, rtolerin, atolerin, ioutin,
		hin, hmaxin, nmaxin, uroundin, safein, faclin, facrin),
	ijac(ijacin), mljac(mljacin), mujac(mujacin), imas(imasin), mlmas(mlmasin),
	mumas(mumasin), nit(nitin), startn(startnin), nind1(nind1in), nind2(nind2in),
	nind3(nind3in), npred(npredin), m1(m1in), m2(m2in), hess(hessin),
	fnewt(fnewtin), quot1(quot1in), quot2(quot2in), thet(thetin), njac(0), ndec(0),
	nsol(0), fac1(0), alphn(0), betan(0), err(0), caljac(true), calhes(true),
	first(true), reject(false)
{

	// Check and change the tolerances
	if (itoler == 0) {
		if ((atoler[0] <= 0.0) || (rtoler[0] <= 10.0*uround)) {
			cout << " tolerances are too small" << endl;
			throw -1;
		}
		else {
			double quot = atoler[0]/rtoler[0];
			rtoler[0] = 0.1*pow(rtoler[0], 2.0/3.0);
			atoler[0] = rtoler[0]*quot;
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			if ((atoler[i] <= 0.0) || (rtoler[i] <= 10.0*uround)) {
				cout << " tolerances('" << i << ") are too small" << endl;
				throw -1;
			}
			else {
				double quot = atoler[i]/rtoler[i];
				rtoler[i] = 0.1*pow(rtoler[i], 2.0/3.0);
				atoler[i] = rtoler[i]*quot;
			}
		}
	}
	
	// initial step length
	if (fabs(h) < 10.0*uround) h = 1.0e-6;
	
	// facl, facr--parameters for step size selection
	if (facl == 0.0) facl = 5.0;
	if (facr == 0.0) facr = 1.0/8.0;
	if ((facl < 1.0) || (facr > 1.0)) {
		cout << " curious input facl, facr = " << facl << facr << endl;
		throw -1;
	}

	// nit--maximal number of Newton iterations
	if (nit == 0) nit = 7;
	if (nit <= 0) {
		cout << " curious input nit = " << nit << endl;
		throw -1;
	}

	// startn--switch for starting values of Newton iterations
	if (startn == 0) startn = false;
	else startn = true;

	// parameters for differential-algebraic components
	if (nind1 == 0) nind1 = n;
	if (nind1 + nind2 + nind3 != n) {
		cout << " curious input for nind1, nind2, nind3 = " <<
			nind1 << nind2 << nind3 << endl;
		throw -1;
	}

	// pred--step size control
	if (npred <= 1) pred = true;
	else pred = false;

	// parameters for second order equations
	nm1 = n - m1;
	if (m1 == 0) m2 = n;
	if (m2 == 0) m2 = m1;
	if ((m1 < 0) || (m2 < 0) || (m1 + m2 > n)) {
		cout << " curious input for m1, m2 = " << m1 << m2 << endl;
		throw -1;
	}

	// fnewt--stopping criterion for Newton's method, usually chosen < 1
	if (fnewt == 0.0) fnewt = max(10.0*uround/rtoler[0], min(0.03, sqrt(rtoler[0])));
	if (fnewt <= uround/rtoler[0]) {
		cout << " curious input for fnewt = " << fnewt << endl;
		throw -1;
	}

	// quot1 and quot2--if quot1 < hnew/hold < quot2, step size = const
	if (quot1 == 0.0) quot1 = 1.0;
	if (quot2 == 0.0) quot2 = 1.2;
	if ((quot1 > 1.0) || (quot2 < 1.0)) {
		cout << " curious input for quot1, quot2 = " << quot1 << quot2 << endl;
		throw -1;
	}

	// thet--decides whether the Jacobian should be recomputed
	if (thet == 0.0) thet = 0.001;
	if (thet >= 1.0) {
		cout << " curious input for thet = " << thet << endl;
		throw -1;
	}

	// implicit, banded or not?
	implct = (imas != 0);
	jband = (mljac < nm1);

// Computation of the row-dimensions of the 2-D arrays
	// Jacobian and matrices e1, e2
	if (jband) {
		ldjac = mljac + mujac + 1;
		lde1 = mljac + ldjac;
	}
	else {
		mljac = nm1;
		mujac = nm1;
		ldjac = nm1;
		lde1 = nm1;
	}

	// mass matrix
	if (implct) {
		if (mlmas != nm1) {
			ldmas = mlmas + mumas + 1;
   			if (jband) ijob = 4;
			else ijob = 3;
		}
		else {
			mumas = nm1;
			ldmas = nm1;
			ijob = 5;
		}
		// bandwith of "mas" not smaller than bandwith of "jac"
		if ((mlmas > mljac) || (mumas > mujac)) {
			cout << "bandwith of 'mas' not smaller than bandwith of 'jac'" << endl;
			throw -1;
		}
	}
	else {
		ldmas = 0;
		if (jband) {
			ijob = 2;
		}
		else {
			ijob = 1;
			if ((n > 2) && (hess != 0)) ijob = 7;
		}
	}
	ldmas = max(1, ldmas);

	// Hessenberg option only for explicit equations with full Jacobian
	if ((implct || jband) && (ijob == 7)) {
		cout << " Hessenberg option only for explicit equations with " <<
			"full Jacobian" << endl;
		throw -1;
	}

	// for second-order equations increase ijob by 10
	if (m1 > 0) ijob += 10;

	// Define constants used in linear algebra routines
	mle = mljac;
	mue = mujac;
	mbjac = mljac + mujac + 1;
	mbb = mlmas + mumas + 1;
	mdiag = mle + mue;
	mdiff = mle + mue - mumas;
	mbdiag = mumas + 1;

// Allocate memory for 1-D arrays
	z1 = new double[n];
	z2 = new double[n];
	z3 = new double[n];
	y0 = new double[n];
	scal = new double[n];
	f1 = new double[n];
	f2 = new double[n];
	f3 = new double[n];
	cont = new double[4*n];
	ip1 = new int[nm1];
	ip2 = new int[nm1];
	iphes = new int[n];

// Allocate memory for 2-D arrays
	fjac = new double*[ldjac];
	fmas = new double*[ldmas];
	e1 = new double*[lde1];
	e2r = new double*[lde1];
	e2i = new double*[lde1];

	for (int i = 0; i < ldjac; i++) fjac[i] = new double[n];

	for (int i = 0; i < ldmas; i++) fmas[i] = new double[n];

	for (int i = 0; i < lde1; i++) {
		e1[i] = new double[nm1];
		e2r[i] = new double[nm1];
		e2i[i] = new double[nm1];
	}

} // Constructor

StiffIntegratorT::StiffIntegratorT(const int nin, double yin[], double xin, double xendin,
	double dxin, const int ijacin, int mljacin, int mujacin, const int imasin,
	int mlmasin, int mumasin): IntegratorT(nin, yin, xin, xendin, dxin),
	ijac(ijacin), mljac(mljacin), mujac(mujacin), imas(imasin), mlmas(mlmasin),
	mumas(mumasin), nit(0), startn(0), nind1(0), nind2(0), nind3(0), npred(0), m1(0),
	m2(0), hess(0), fnewt(0), quot1(0), quot2(0), thet(0), njac(0), ndec(0), nsol(0),
	fac1(0), alphn(0), betan(0), err(0), caljac(true), calhes(true), first(true),
	reject(false)
{

	// Check and change the tolerances
	if (itoler == 0) {
		if ((atoler[0] <= 0.0) || (rtoler[0] <= 10.0*uround)) {
			cout << " tolerances are too small" << endl;
			throw -1;
		}
		else {
			double quot = atoler[0]/rtoler[0];
			rtoler[0] = 0.1*pow(rtoler[0], 2.0/3.0);
			atoler[0] = rtoler[0]*quot;
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			if ((atoler[i] <= 0.0) || (rtoler[i] <= 10.0*uround)) {
				cout << " tolerances('" << i << ") are too small" << endl;
				throw -1;
			}
			else {
				double quot = atoler[i]/rtoler[i];
				rtoler[i] = 0.1*pow(rtoler[i], 2.0/3.0);
				atoler[i] = rtoler[i]*quot;
			}
		}
	}
	
	// initial step length
	if (fabs(h) < 10.0*uround) h = 1.0e-6;
	
	// facl, facr--parameters for step size selection
	if (facl == 0.0) facl = 5.0;
	if (facr == 0.0) facr = 1.0/8.0;
	if ((facl < 1.0) || (facr > 1.0)) {
		cout << " curious input facl, facr = " << facl << facr << endl;
		throw -1;
	}

	// nit--maximal number of Newton iterations
	if (nit == 0) nit = 7;
	if (nit <= 0) {
		cout << " curious input nit = " << nit << endl;
		throw -1;
	}

	// startn--switch for starting values of Newton iterations
	if (startn == 0) startn = false;
	else startn = true;

	// parameters for differential-algebraic components
	if (nind1 == 0) nind1 = n;
	if (nind1 + nind2 + nind3 != n) {
		cout << " curious input for nind1, nind2, nind3 = " <<
			nind1 << nind2 << nind3 << endl;
		throw -1;
	}

	// pred--step size control
	if (npred <= 1) pred = true;
	else pred = false;

	// parameters for second order equations
	nm1 = n - m1;
	if (m1 == 0) m2 = n;
	if (m2 == 0) m2 = m1;
	if ((m1 < 0) || (m2 < 0) || (m1 + m2 > n)) {
		cout << " curious input for m1, m2 = " << m1 << m2 << endl;
		throw -1;
	}

	// fnewt--stopping criterion for Newton's method, usually chosen < 1
	if (fnewt == 0.0) fnewt = max(10.0*uround/rtoler[0], min(0.03, sqrt(rtoler[0])));
	if (fnewt <= uround/rtoler[0]) {
		cout << " curious input for fnewt = " << fnewt << endl;
		throw -1;
	}

	// quot1 and quot2--if quot1 < hnew/hold < quot2, step size = const
	if (quot1 == 0.0) quot1 = 1.0;
	if (quot2 == 0.0) quot2 = 1.2;
	if ((quot1 > 1.0) || (quot2 < 1.0)) {
		cout << " curious input for quot1, quot2 = " << quot1 << quot2 << endl;
		throw -1;
	}

	// thet--decides whether the Jacobian should be recomputed
	if (thet == 0.0) thet = 0.001;
	if (thet >= 1.0) {
		cout << " curious input for thet = " << thet << endl;
		throw -1;
	}

	// implicit, banded or not?
	implct = (imas != 0);
	jband = (mljac < nm1);

// Computation of the row-dimensions of the 2-D arrays
	// Jacobian and matrices e1, e2
	if (jband) {
		ldjac = mljac + mujac + 1;
		lde1 = mljac + ldjac;
	}
	else {
		mljac = nm1;
		mujac = nm1;
		ldjac = nm1;
		lde1 = nm1;
	}

	// mass matrix
	if (implct) {
		if (mlmas != nm1) {
			ldmas = mlmas + mumas + 1;
   			if (jband) ijob = 4;
			else ijob = 3;
		}
		else {
			mumas = nm1;
			ldmas = nm1;
			ijob = 5;
		}
		// bandwith of "mas" not smaller than bandwith of "jac"
		if ((mlmas > mljac) || (mumas > mujac)) {
			cout << "bandwith of 'mas' not smaller than bandwith of 'jac'" << endl;
			throw -1;
		}
	}
	else {
		ldmas = 0;
		if (jband) {
			ijob = 2;
		}
		else {
			ijob = 1;
			if ((n > 2) && (hess != 0)) ijob = 7;
		}
	}
	ldmas = max(1, ldmas);

	// Hessenberg option only for explicit equations with full Jacobian
	if ((implct || jband) && (ijob == 7)) {
		cout << " Hessenberg option only for explicit equations with " <<
			"full Jacobian" << endl;
		throw -1;
	}

	// for second-order equations increase ijob by 10
	if (m1 > 0) ijob += 10;

	// Define constants used in linear algebra routines
	mle = mljac;
	mue = mujac;
	mbjac = mljac + mujac + 1;
	mbb = mlmas + mumas + 1;
	mdiag = mle + mue;
	mdiff = mle + mue - mumas;
	mbdiag = mumas + 1;

// Allocate memory for 1-D arrays
	z1 = new double[n];
	z2 = new double[n];
	z3 = new double[n];
	y0 = new double[n];
	scal = new double[n];
	f1 = new double[n];
	f2 = new double[n];
	f3 = new double[n];
	cont = new double[4*n];
	ip1 = new int[nm1];
	ip2 = new int[nm1];
	iphes = new int[n];

// Allocate memory for 2-D arrays
	fjac = new double*[ldjac];
	fmas = new double*[ldmas];
	e1 = new double*[lde1];
	e2r = new double*[lde1];
	e2i = new double*[lde1];

	for (int i = 0; i < ldjac; i++) fjac[i] = new double[n];

	for (int i = 0; i < ldmas; i++) fmas[i] = new double[n];

	for (int i = 0; i < lde1; i++) {
		e1[i] = new double[nm1];
		e2r[i] = new double[nm1];
		e2i[i] = new double[nm1];
	}

} // Constructor

// Destructor
StiffIntegratorT::~StiffIntegratorT()
{
	delete [] z1;
	delete [] z2;
	delete [] z3;
	delete [] y0;
	delete [] scal;
	delete [] f1;
	delete [] f2;
	delete [] f3;
	delete [] cont;
	delete [] ip1;
	delete [] ip2;
	delete [] iphes;
	for (int i = 0; i < lde1; i++) {
		delete [] e1[i];
		delete [] e2r[i];
		delete [] e2i[i];
	}
	for (int i = 0; i < ldjac; i++) delete [] fjac[i];
	for (int i = 0; i < ldmas; i++) delete [] fmas[i];
	delete [] e1;
	delete [] e2r;
	delete [] e2i;
	delete [] fjac;
	delete [] fmas;

}

void StiffIntegratorT::Integrate(double eps)
{

	int idid = CoreIntegrator(eps);

	if (idid < 0) {
		cout << " Computation failed " << endl;
		return;
	}

	// restore tolerances
	if (itoler == 0) {
		double quot = atoler[0]/rtoler[0];
		rtoler[0] = pow(10.0*rtoler[0], 1.5);
		atoler[0] = rtoler[0]*quot;
	}
	else {
		for (int i = 0; i < n; i++) {
			double quot = atoler[i]/rtoler[i];
			rtoler[i] = pow(10.0*rtoler[i], 1.5);
			atoler[i] = rtoler[i]*quot;
		}
	}

	// print final solution
	if (iout == 1) {
		cout << "Step " << naccpt << ": t = " << setw(5) <<
				setprecision(2) << xend << "  y = ";
		for (int i = 0; i < n; i++)
			cout << setw(10) << setprecision(8) << y[i] << "  ";
		cout << endl;
	}

	return;

} // Integrate

// core integrator for radau5

// return value of CoreIntegrator
//  1  computation successful,
//  2  comput. successful (interrupted by SolutionOutput)
// -1  error in linear algebra routines,
// -2  larger nmax is needed,
// -3  step size becomes too small,
// -4  matrix is repeatedly singular.
// -5 not enough memory
int StiffIntegratorT::CoreIntegrator(double eps)
{
// constants
 	const double t11 = 9.1232394870892942792e-02;
	const double t12 = -0.14125529502095420843;
	const double t13 = -3.0029194105147424492e-02;
	const double t21 = 0.24171793270710701896;
	const double t22 = 0.20412935229379993199;
	const double t23 = 0.38294211275726193779;
	const double t31 = 0.96604818261509293619;
	const double ti11 = 4.3255798900631553510;
	const double ti12 = 0.33919925181580986954;
	const double ti13 = 0.54177053993587487119;
	const double ti21 = -4.1787185915519047273;
	const double ti22 = -0.32768282076106238708;
	const double ti23 = 0.47662355450055045196;
	const double ti31 = -0.50287263494578687595;
	const double ti32 = 2.5719269498556054292;
	const double ti33 = -0.59603920482822492497;

	const double sq6 = sqrt(6.0);
	const double c1 = (4.0 - sq6)/10.0;
	const double c2 = (4.0 + sq6)/10.0;
	const double c1m1 = c1 - 1.0;
	const double c2m1 = c2 - 1.0;
	const double c1mc2 = c1 - c2;
	const double u1 = 1.0/((6.0 + pow(81.0, 1.0/3.0) - pow(9.0, 1.0/3.0))/30.0);
	double alph = (12.0 - pow(81.0, 1.0/3.0) + pow(9.0, 1.0/3.0))/60.0;
	double beta = (pow(81.0, 1.0/3.0) + pow(9.0, 1.0/3.0))*sqrt(3.0)/60.0;
	const double cno = alph*alph + beta*beta;

	alph = alph/cno;
	beta = beta/cno;

	const double posneg = copysign(1.0, xend-x);
	const double hmaxn = min(fabs(hmax), fabs(xend - x));
	const double cfac = safe*(1 + 2*nit);

	// compute mass matrix for implicit case
	if (implct) Mass(fmas);

	h = min(fabs(h), hmaxn);
	h = copysign(h, posneg);
	hold = h;

	bool last = false;

	if ((x + h*1.0001 - xend)*posneg >= 0.0) {
		h = xend - x;
		last = true;
	}

	double hopt = h;
	double faccon = 1.0;

	if (iout != 0) {
		int irtrn = 1;
		for (int i = 0; i < n; i++) cont[i] = y[i];
		irtrn = SolutionOutput();
		if (irtrn < 0) {
			cout << " exit of RADAU5 at x = " << x << endl;
			return 2;
		}
	}

	if (itoler == 0) {
		for (int i = 0; i < n; i++)
			scal[i] = atoler[0] + rtoler[0]*fabs(y[i]);
	}
	else {
		for (int i = 0; i < n; i++)
			scal[i] = atoler[i] + rtoler[i]*fabs(y[i]);
	}

	Function(x, y, y0,eps);
	nfcn++;

	double hhfac = h;
	double hacc, erracc, thqold;
	int nsing = 0, ier = 0;

// basic integration step
	ComputeJacobian(eps);
	bool loop = true;
	while (loop) {
		loop = false;
		// compute the matrices e1 and e2 and their decompositions
		fac1 = u1/h;
		alphn = alph/h;
		betan = beta/h;

		ier = DecompReal();

		if (ier != 0) {
			if (ier == -1) return -1;
			nsing++;
			if (nsing >= 5) {
				cout << " exit of RADAU5 at x = " << x << endl;
				cout << " matrix is repeatedly singular, ier= " << ier << endl;
				return -4;
			}
			h *= 0.5;
			hhfac = 0.5;
			reject = true;
			last = false;
			if (!caljac) ComputeJacobian(eps);
			loop = true;
			continue;
		}

		ier = DecompComplex();

		if (ier != 0) {
			if (ier == -1) return(-1);
			nsing++;
			if (nsing >= 5) {
				cout << " exit of RADAU5 at x = " << x << endl;
				cout << " matrix is repeatedly singular, ier= " << ier << endl;
				return -4;
			}
			h *= 0.5;
			hhfac = 0.5;
			reject = true;
			last = false;
			if (!caljac) ComputeJacobian(eps);
			loop = true;
			continue;
		}
		ndec++;

		while (true) {
			nstep++;
			if (nstep >= nmax) {
				cout << " exit of RADAU5 at x = " << x << endl;
				cout << " more than nmax = " << nmax << " steps are needed" << endl;
				return -2;
			}

			if (0.1*fabs(h) <= fabs(x)*uround) {
				cout << " exit of RADAU5 at x = " << x << endl;
				cout << " step size too small, h = " << h << endl;
				return -3;
			}

			// check the index of the problem
			if (nind2 != 0) { // is index 2
				for (int i = nind1; i < nind1 + nind2 ; i++)
					scal[i] = scal[i]/hhfac;
			}

			if (nind3 != 0) { // is index 3
				for (int i = nind1 + nind2; i < nind1 + nind2 + nind3; i++)
					scal[i]=scal[i]/(hhfac*hhfac);
			}

			double xph = x + h;
			//  starting values for Newton iteration
			if (first || startn) {
				for (int i = 0; i < n; i++)
					z1[i] = z2[i] = z3[i] = f1[i] = f2[i] = f3[i] = 0.0;
			}
			else {
				double c3q = h/hold;
				double c1q = c1*c3q;
				double c2q = c2*c3q;
				double ak1, ak2, ak3;
				for (int i = 0; i < n; i++) {
					ak1 = cont[i+n];
					ak2 = cont[i+2*n];
					ak3 = cont[i+3*n];
					z1[i] = c1q*(ak1 + (c1q - c2m1)*(ak2 + (c1q - c1m1)*ak3));
					z2[i] = c2q*(ak1 + (c2q - c2m1)*(ak2 + (c2q - c1m1)*ak3));
					z3[i] = c3q*(ak1 + (c3q - c2m1)*(ak2 + (c3q - c1m1)*ak3));
					f1[i] = ti11*z1[i] + ti12*z2[i] + ti13*z3[i];
					f2[i] = ti21*z1[i] + ti22*z2[i] + ti23*z3[i];
					f3[i] = ti31*z1[i] + ti32*z2[i] + ti33*z3[i];
				}
			}

			//  loop for the simplified Newton iteration
			int newt = 0;
			faccon = pow(max(faccon, uround), 0.8);
			double theta = fabs(thet);
			double dyno, dynold;

			while (true) {
				if (newt >= nit) {
					if (ier != 0) {
						nsing++;
						if (nsing >= 5) {
							cout << " exit of RADAU5 at x = " << x << endl;
							cout << " matrix is repeatedly singular, ier= " << ier << endl;
							return -4;
						}
					}
					h *= 0.5;
					hhfac = 0.5;
					reject = true;
					last = false;
					if (!caljac) ComputeJacobian(eps);
					loop = true;
					break;
				}
				// compute the right-hand side
				for (int i = 0; i < n; i++) cont[i] = y[i] + z1[i];
					Function(x+c1*h, cont, z1,eps);

				for (int i = 0; i < n; i++) cont[i] = y[i] + z2[i];
					Function(x+c2*h, cont, z2, eps);

				for (int i = 0; i < n; i++) cont[i] = y[i] + z3[i];
					Function(xph, cont, z3, eps);

				nfcn += 3;

				// solve the linear systems
				for (int i = 0; i < n; i++) {
					double a1 = z1[i];
					double a2 = z2[i];
					double a3 = z3[i];
					z1[i] = ti11*a1 + ti12*a2 + ti13*a3;
					z2[i] = ti21*a1 + ti22*a2 + ti23*a3;
					z3[i] = ti31*a1 + ti32*a2 + ti33*a3;
				}
				ier = LinearSolve();
				if (ier == -1) return -1;
				nsol++;
				newt++;
				dyno = 0.0;
				double denom;
				for (int i = 0; i < n; i++) {
					denom = scal[i];
					dyno = dyno + pow(z1[i]/denom, 2) + pow(z2[i]/denom, 2) +
						pow(z3[i]/denom, 2);
				}
				dyno = sqrt(dyno/(3*n));
				// bad convergence or number of iterations to large
				if ((newt > 1) && (newt < nit)) {
					double thq = dyno/dynold;
					if (newt == 2) theta = thq;
					else theta = sqrt(thq*thqold);
					thqold = thq;
					if (theta < 0.99) {
						faccon = theta/(1.0 - theta);
						double dyth = faccon*dyno*pow(theta, nit-1-newt)/fnewt;
						if (dyth >= 1.0) {
							double qnewt = max(1.0e-4, min(20.0, dyth));
							hhfac = 0.8*pow(qnewt, -1.0/(4.0+nit-1-newt));
							h *= hhfac;
							reject = true;
							last = false;
							if (caljac) ComputeJacobian(eps);
							loop = true;
							break;
						}
					}
					else {
						if (ier != 0) {
							nsing++;
							if (nsing >= 5) {
								cout << " exit of RADAU5 at x = " << x << endl;
								cout << " matrix is repeatedly singular, ier= " << ier << endl;
								return -4;
							}
						}
						h *= 0.5;
						hhfac = 0.5;
						reject = true;
						last = false;
						if (!caljac) ComputeJacobian(eps);
						loop = true;
						break;
					}
				}
				dynold = max(dyno, uround);
				for (int i = 0; i < n; i++) {
					f1[i] = f1[i] + z1[i];
					f2[i] = f2[i] + z2[i];
					f3[i] = f3[i] + z3[i];
					z1[i] = t11*f1[i] + t12*f2[i] + t13*f3[i];
					z2[i] = t21*f1[i] + t22*f2[i] + t23*f3[i];
					z3[i] = t31*f1[i] + f2[i];
				}
				if (faccon*dyno <= fnewt) break;
			}

			if (loop) break;

			// error estimation
			err = 0.0;
			ier = ErrorEstimate(eps);
			if (ier == -1) return -1;

			// computation of hnew -- require 0.2 <= hnew/h <= 8.
			double fac = min(safe, cfac/(newt+2*nit));
			double quot = max(facr, min(facl, pow(err, 0.25)/fac));
			double hnew = h/quot;

			//  is the error small enough ?
			if (err < 1.0) {
			// step is accepted
				first = false;
				naccpt++;
				if (pred) {
					// predictive controller of Gustafsson
					if (naccpt > 1) {
						double facgus = (hacc/(h))*pow(err*err/erracc, 0.25)/safe;
						facgus = max(facr, min(facl, facgus));
						quot = max(quot,facgus);
						hnew = h/quot;
					}
					hacc = h;
					erracc = max(1.0e-2, err);
				}
				xold = x;
				hold = h;
				x = xph;
				double ak, acont3;
				for (int i = 0; i < n; i++) {
					y[i] = y[i] + z3[i];
					cont[i+n] = (z2[i] - z3[i])/c2m1;
					ak = (z1[i] - z2[i])/c1mc2;
					acont3 = z1[i]/c1;
					acont3 = (ak - acont3)/c2;
					cont[i+2*n] = (ak - cont[i+n])/c1m1;
					cont[i+3*n] = cont[i+2*n] - acont3;
				}
				if (itoler == 0) {
					for (int i = 0; i < n; i++) scal[i] = atoler[0] + rtoler[0]*fabs(y[i]);
				}
				else {
					for (int i = 0; i < n; i++) scal[i] = atoler[i] + rtoler[i]*fabs(y[i]);
				}
				if (iout != 0) {
					for (int i = 0; i < n; i++) cont[i] = y[i];
					int irtrn = 1;
					irtrn = SolutionOutput();
					if (irtrn < 0) {
						cout << " exit of RADAU5 at x = " << x << endl;
						return 2;
					}
				}
				caljac = false;
				if (last) {
					h = hopt;
					return 1;
				}
				Function(x, y, y0,eps);
				nfcn++;
				hnew = posneg*min(fabs(hnew), hmaxn);
				hopt = min(h, hnew);
				if (reject) hnew = posneg*min(fabs(hnew), fabs(h));
				reject = false;
				if ((x + hnew/quot1 - xend)*posneg >= 0.0) {
					h = xend - x;
					last = true;
				}
				else {
					double qt = hnew/(h);
					hhfac = h;
					if ((theta <= thet) && (qt >= quot1) && (qt <= quot2)) continue;
					h = hnew;
				}
				hhfac = h;
				if (theta > thet) ComputeJacobian(eps);
				loop = true;
			}
			else {
			// step is rejected
				reject = true;
				last = false;
				if (first) {
					h *= 0.1;
					hhfac = 0.1;
				}
				else {
					hhfac = hnew/(h);
					h = hnew;
				}
				if (naccpt >= 1) nrejct++;
				if (!caljac) ComputeJacobian(eps);
				loop = true;
			}
			break;
		}
	}
	
	return 1;
}  // CoreIntegrator

void StiffIntegratorT::ComputeJacobian(double eps)
{
	njac++;
	if (ijac == 0) {
	// compute Jacobian matrix numerically
		if (jband) {
		// Jacobian is banded
			int mujacp = mujac + 1;
			int md = min(mbjac, m2);
			for (int mm1 = 0; mm1 < m1/m2 + 1; mm1++) {
				for (int k = 0; k < md; k++) {
					int j = k + mm1*m2;
					while (true) {
						f1[j] = y[j];
						f2[j] = sqrt(uround*max(1.0e-5, fabs(y[j])));
						y[j] = y[j] + f2[j];
						j += md;
						if (j > (mm1+1)*m2 - 1) break;
					}
					Function(x, y, cont, eps);
					j = k + mm1*m2;
					int j1 = k;
					int lbeg = max(0, j1 - mujac) + m1;
					int lend, mujacj;
					while (true) {
						lend = min(m2, j1 + mljac) + m1;
						y[j] = f1[j];
						mujacj = mujacp - j1 - m1 - 1;
						for (int l = lbeg; l <= lend; l++) {
							fjac[l+mujacj][j] = (cont[l] - y0[l])/f2[j];
						}
						j += md;
						j1 += md;
						lbeg = lend + 1;
						if (j > (mm1+1)*m2 - 1) break;
					}
				}
			}
		}
		else {
		// Jacobian is full
			double delt, ysafe;
			for (int i = 0; i < n; i++) {
				ysafe = y[i];
				delt = sqrt(uround*max(1.0e-5, fabs(ysafe)));
				y[i] = ysafe + delt;
				Function(x, y, cont, eps);
				for (int j = m1; j < n; j++)
					fjac[j-m1][i] = (cont[j] - y0[j])/delt;
				y[i] = ysafe;
			}
		}
	}
	else {
	// compute Jacobian matrix analytically
		Jacobian(x, y, fjac,eps);
	}

	caljac = true;
	calhes = true;

	return;
}

double StiffIntegratorT::ContinuousOutput(unsigned i)
{

// ----------------------------------------------------------
//     This function can be used for continuous output. It provides an
//     approximation to the i-th component of the solution at xd.
//     It gives the value of the collocation polynomial, defined for
//     the last successfully computed step.
// ----------------------------------------------------------

	double s, sq6, c1, c2, c2m1, c1m1;

	sq6 = sqrt(6.0);
	c1 = (4.0 - sq6)/10.0;
	c2 = (4.0 + sq6)/10.0;
	c1m1 = c1 - 1.0;
	c2m1 = c2 - 1.0;

	s = (xd - x)/hold;

	return (cont[i] + s*(cont[i+n] + (s - c2m1)*(cont[i+2*n] +
		(s - c1m1)*cont[i+3*n])));

}


int StiffIntegratorT::DecompReal()
{
	int mm, ier;
	double sum;

	switch (ijob)
	{
		case (1):
			// mass = identity, Jacobian a full matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < n; i++) {
					e1[i][j] = -fjac[i][j];
				}
				e1[j][j] += fac1;
			}
			ier = dec(n, e1, ip1);
			break;

		case (2):
			// mass = identity, Jacobian a banded matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < mbjac; i++) {
					e1[i+mle][j] = -fjac[i][j];
				}
				e1[mdiag][j] += fac1;
			}
			ier = decb(n, e1, mle, mue, ip1);
			break;

		case (3):
			// mass is a banded matrix, Jacobian a full matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < n; i++)
					e1[i][j] = -fjac[i][j];
				for (int i = max(0, j-mumas); i < min(n, j+mlmas+1); i++)
					e1[i][j] += fac1*fmas[i-j+mbdiag-1][j];
			}
			ier = dec(n, e1, ip1);
			break;

		case (4):
			// mass is a banded matrix, Jacobian a banded matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < mbjac; i++)
					e1[i+mle][j] = -fjac[i][j];
				for (int i = 0; i < mbb; i++)
					e1[i+mdiff][j] += fac1*fmas[i][j];
			}
			ier = decb(n, e1, mle, mue, ip1);
			break;

		case (5):
			// mass is a full matrix, Jacobian a full matrix
			for (int j = 0; j < n; j++)
				for (int i = 0; i < n; i++)
					e1[i][j] = fmas[i][j]*fac1 - fjac[i][j];
			ier = dec(n, e1, ip1);
			break;

		case (6):
			// mass is a full matrix, Jacobian a banded matrix
			// This option is not provided
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);

		case (7):
			// mass = identity, Jacobian a full matrix, Hessenberg-option
			if (calhes) elmhes(n, 0, n, fjac, iphes);
			calhes = false;
			for (int j = 0; j < n - 1; j++) e1[j+1][j] = -fjac[j+1][j];
			for (int j = 0; j < n; j++) {
				for (int i = 0; i <= j; i++) e1[i][j] = -fjac[i][j];
				e1[j][j] += fac1;
			}
			ier = dech(n, e1, 1, ip1);
			break;

		case (11):
			// mass = identity, Jacobian a full matrix, second order
			for (int j = 0; j < nm1; j++) {
				for (int i = 0; i < nm1; i++) {
					e1[i][j] = -fjac[i][j+m1];
				}
				e1[j][j] += fac1;
			}
			break;

		case (12):
			// mass = identity, Jacobian a banded matrix, second order
			for (int j = 0; j < nm1; j++) {
				for (int i = 0; i < mbjac; i++)
					e1[i+mle][j] = -fjac[i][j+m1];
				e1[mdiag][j] += fac1;
			}
			break;

		case (13):
			// mass is a banded matrix, Jacobian a full matrix, second order
			for (int j = 0; j < nm1; j++) {
				for (int i = 0; i < nm1; i++)
					e1[i][j] = -fjac[i][j+m1];
				for (int i = max(0, j-mumas); i < min(n, j+mlmas+1); i++)
					e1[i][j] += fac1*fmas[i-j+mbdiag-1][j];
			}
			break;

		case (14):
			// mass is a banded matrix, Jacobian a banded matrix, second order
			for (int j = 0; j < nm1; j++) {
		 		for (int i = 0; i < mbjac; i++)
					e1[i+mle][j] = -fjac[i][j+m1];
				for (int i = 0; i < mbb; i++)
					e1[i+mdiff][j] += fac1*fmas[i][j];
			}
			break;

		case (15):
			// mass is a full matrix, Jacobian a full matrix, second order
			for (int j = 0; j < nm1; j++)
				for (int i = 0; i < nm1; i++)
					e1[i][j] = fmas[i][j]*fac1 - fjac[i][j+m1];
			break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);
	}

	switch (ijob)
	{
		case (1):
		case (2):
		case (3):
		case (4):
		case (5):
		case (7):
			break;

		case (11):
		case (13):
		case (15):
			mm = m1/m2;
			for (int j = 0; j < m2; j++) {
				for (int i = 0; i < nm1; i++) {
					sum = 0.0;
					for (int k = 0; k < mm; k++)
						sum = (sum + fjac[i][j+k*m2])/fac1;
					e1[i][j] -= sum;
				}
			}
			ier = dec(nm1, e1, ip1);
			break;

		case (12):
		case (14):
			mm = m1/m2;
			for (int j = 0; j < m2; j++) {
				for (int i = 0; i < mbjac; i++) {
					sum = 0.0;
					for (int k = 0; k < mm; k++)
						sum = (sum + fjac[i][j+k*m2])/fac1;
					e1[i+mle][j] -= sum;
				}
			}
			ier = decb(nm1, e1, mle, mue, ip1);
			break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);
	}

	return (ier);

} // DecompReal


int StiffIntegratorT::DecompComplex()
{

	int mm, ier;
	double bb, ffma, abno, alp, bet, sumr, sumi, sums;

	switch (ijob) {
		case (1):
			// mass = identity, Jacobian a full matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < n; i++) {
					e2r[i][j] = -fjac[i][j];
					e2i[i][j] = 0.0;
				}
				e2r[j][j] += alphn;
				e2i[j][j] = betan;
			}
			ier = decc(n, e2r, e2i, ip2);
			break;

		case (2):
			// mass = identiy, Jacobian a banded matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < mbjac; i++) {
					e2r[i+mle][j] = -fjac[i][j];
					e2i[i+mle][j] = 0.0;
				}
				e2r[mdiag][j] += alphn;
				e2i[mdiag][j] = betan;
			}
			ier = decbc(n, e2r, e2i, mle, mue, ip2);
			break;

		case (3):
			// mass is a banded matrix, Jacobian a full matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < n; i++) {
					e2r[i][j] = -fjac[i][j];
					e2i[i][j] = 0.0;
				}
			}
			for (int j = 0; j < n; j++) {
				for (int i = max(0, j-mumas); i < min(n, j+mlmas+1); i++) {
					bb = fmas[i-j+mbdiag-1][j];
					e2r[i][j] += alphn*bb;
					e2i[i][j] = betan*bb;
				}
			}
			ier = decc(n, e2r, e2i, ip2);
			break;

		case (4):
			// mass is a banded matrix, Jacobian a banded matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < mbjac; i++) {
					e2r[i+mle][j] = -fjac[i][j];
					e2i[i+mle][j] = 0.0;
				}
				for (int i = max(0, mumas-j); i < min(mbb, mumas-j+n); i++) {
					bb = fmas[i][j];
					e2r[i+mdiff][j] += alphn*bb;
					e2i[i+mdiff][j] = betan*bb;
				}
			}
			ier = decbc(n, e2r, e2i, mle, mue, ip2);
			break;

		case (5):
			// mass is a full matrix, Jacobian a full matrix
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < n; i++) {
					bb = fmas[i][j];
					e2r[i][j] = alphn*bb - fjac[i][j];
					e2i[i][j] = betan*bb;
				}
			}
			ier = decc(n, e2r, e2i, ip2);
			break;

		case (6):
			// mass is a full matrix, Jacobian a banded matrix
			// This option is not provided
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);

		case (7):
			// mass = identity, Jacobian a full matrix, Hessenberg-option
			for (int j = 0; j < n - 1; j++) {
				e2r[j+1][j] = -fjac[j+1][j];
				e2i[j+1][j] = 0.0;
			}
			for (int j = 0; j < n; j++) {
				for (int i = 0; i <= j; i++) {
					e2i[i][j] = 0.0;
					e2r[i][j] = -fjac[i][j];
				}
				e2r[j][j] += alphn;
				e2i[j][j] = betan;
			}
			ier = dechc(n, e2r, e2i, 1, ip2);
			break;

		case (11):
			// mass = identity, Jacobian a full matrix, second order
			for (int j = 0; j < nm1; j++) {
				for (int i = 0; i < nm1; i++) {
					e2r[i][j] = -fjac[i][j+m1];
					e2i[i][j] = 0.0;
				}
				e2r[j][j] += alphn;
				e2i[j][j] = betan;
			}
			break;

		case (12):
			// mass = identity, Jacobian a banded matrix, second order
	  		for (int j = 0; j < nm1; j++) {
				for (int i = 0; i < mbjac; i++) {
					e2r[i+mle][j] = -fjac[i][j+m1];
					e2i[i+mle][j] = 0.0;
				}
				e2r[mdiag][j] += alphn;
				e2i[mdiag][j] += betan;
			}
			break;

		case (13):
			// mass is a banded matrix, Jacobian a full matrix, second order
			for (int j = 0; j < nm1; j++) {
				for (int i = 0; i < nm1; i++) {
					e2r[i][j] = -fjac[i][j+m1];
					e2i[i][j] = 0.0;
				}
				for (int i = max(0, j-mumas); i < min(nm1, j+mlmas+1); i++) {
					ffma = fmas[i-j+mbdiag-1][j];
					e2r[j][j] += alphn*ffma;
					e2i[j][j] += betan*ffma;
				}
			}
			break;

		case (14):
			// mass is a banded matrix, Jacobian a banded matrix, second order
	  		for (int j = 0; j < nm1; j++) {
				for (int i = 0; i < mbjac; i++) {
					e2r[i+mle][j] = -fjac[i][j+m1];
					e2i[i+mle][j] = 0.0;
				}
				for (int i = 0; i < mbb; i++) {
					ffma = fmas[i][j];
					e2r[i+mdiff][j] += alphn*ffma;
					e2i[i+mdiff][j] += betan*ffma;
				}
			}
			break;

		case (15):
			// mass is a full matrix, Jacobian a full matrix, second order
			for (int j = 0; j < nm1; j++) {
				for (int i = 0; i < nm1; i++) {
					e2r[i][j] = alphn*fmas[i][j] - fjac[i][j+m1];
					e2i[i][j] = betan*fmas[i][j];
				}
			}
			break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);

	}

	switch (ijob) {
		case (1):
		case (2):
		case (3):
		case (4):
		case (5):
		case (7):
			break;

		case (11):
		case (13):
		case (15):
			mm = m1/m2;
			abno = alphn*alphn + betan*betan;
			alp = alphn/abno;
			bet = betan/abno;
			for (int j = 0; j < m2; j++) {
				for (int i = 0; i < nm1; i++) {
					sumr = sumi = 0.0;
					for (int k = 0; k < mm; k++) {
						sums = sumr + fjac[i][j+k*m2];
						sumr = sums*alp + sumi*bet;
						sumi = sumi*alp - sums*bet;
					}
					e2r[i][j] -= sumr;
					e2i[i][j] -= sumi;
				}
			}
			ier = decc(nm1, e2r, e2i, ip2);
			break;

		case (12):
		case (14):
			mm = m1/m2;
			abno = alphn*alphn + betan*betan;
			alp = alphn/abno;
			bet = betan/abno;
			for (int j = 0; j < m2; j++) {
				for (int i = 0; i < mbjac; i++) {
					sumr = sumi = 0.0;
					for (int k = 0; k < mm; k++) {
						sums = sumr + fjac[i][j+k*m2];
						sumr = sums*alp + sumi*bet;
						sumi = sumi*alp - sums*bet;
					}
					e2r[i+mle][j] -= sumr;
					e2i[i+mle][j] -= sumi;
				}
			}
			ier = decbc(nm1, e2r, e2i, mle, mue, ip2);
			break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);
	}

	return (ier);

} // DecompComplex


int StiffIntegratorT::LinearSolve()
{

	int mm, mp, mp1, ii, jkm, mpi, ier = 0;
	double abno, bb, e1imp, s1, s2, s3, sum1, sum2, sum3, sumh;
	double ffja, z2i, z3i, zsafe;

	switch (ijob) {
		case (1):
			// mass = identity, Jacobian a full matrix
			for (int i = 0; i < n; i++) {
				s2 = -f2[i];
				s3 = -f3[i];
				z1[i] -= f1[i]*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			sol(n, e1, z1, ip1);
			solc(n, e2r, e2i, z2, z3, ip2);
			break;

		case (2):
			// mass = identity, Jacobian a banded matrix
			for (int i = 0; i < n; i++) {
				s2 = -f2[i];
				s3 = -f3[i];
				z1[i] -= f1[i]*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			solb(n, e1, mle, mue, z1, ip1);
			solbc(n, e2r, e2i, mle, mue, z2, z3, ip2);
			break;

		case (3):
			// mass is a banded matrix, Jacobian a full matrix
			for (int i = 0; i < n; i++) {
				s1 = s2 = s3 = 0.0;
				for (int j = max(0, i-mlmas); j < min(n, i+mumas+1); j++) {
					bb = fmas[i-j+mbdiag-1][j];
					s1 -= bb*f1[j];
					s2 -= bb*f2[j];
					s3 -= bb*f3[j];
				}
				z1[i] += s1*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			sol(n, e1, z1, ip1);
			solc(n, e2r, e2i, z2, z3, ip2);
			break;

		case (4):
			// mass is a banded matrix, Jacobian a banded matrix
			for (int i = 0; i < n; i++) {
				s1 = s2 = s3 = 0.0;
				for (int j = max(0, i-mlmas); j < min(n, i+mumas+1); j++) {
					bb = fmas[i-j+mbdiag-1][j];
					s1 -= bb*f1[j];
					s2 -= bb*f2[j];
					s3 -= bb*f3[j];
				}
				z1[i] += s1*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			solb(n, e1, mle, mue, z1, ip1);
			solbc(n, e2r, e2i, mle, mue, z2, z3, ip2);
			break;

		case (5):
			// mass is a full matrix, Jacobian a full matrix
			for (int i = 0; i < n; i++) {
				s1 = s2 = s3 = 0.0;
				for (int j = 0; j < n; j++) {
					bb = fmas[i][j];
					s1 -= bb*f1[j];
					s2 -= bb*f2[j];
					s3 -= bb*f3[j];
				}
				z1[i] += s1*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			sol(n, e1, z1, ip1);
			solc(n, e2r, e2i, z2, z3, ip2);
			break;

		case (6):
			// mass is a full matrix, Jacobian a banded matrix
			// This option is not provided
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);

		case (7):
			// mass = identity, Jacobian a full matrix, Hessenberg-option
			for (int i = 0; i < n; i++) {
				s2 = -f2[i];
				s3 = -f3[i];
				z1[i] -= f1[i]*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			for (int mm1 = n - 3; mm1 >= 0; mm1--) {
				mp = n - mm1 - 2;
				mp1 = mp - 1;
				ii = iphes[mp];
				if (ii != mp) {
					zsafe = z1[mp];
					z1[mp] = z1[ii];
					z1[ii] = zsafe;
					zsafe = z2[mp];
					z2[mp] = z2[ii];
					z2[ii] = zsafe;
					zsafe = z3[mp];
					z3[mp] = z3[ii];
					z3[ii] = zsafe;
				}
				for (int i = mp + 1; i < n ; i++) {
					e1imp = fjac[i][mp1];
					z1[i] -= e1imp*z1[mp];
					z2[i] -= e1imp*z2[mp];
					z3[i] -= e1imp*z3[mp];
				}
			}
			solh(n, e1, 1, z1, ip1);
			solhc(n, e2r, e2i, 1, z2, z3, ip2);
			for (int mm1 = 0; mm1 < n - 2; mm1++) {
				mp = n - mm1 - 2;
				mp1 = mp - 1;
				for (int i = mp; i < n; i++) {
					e1imp = fjac[i][mp1];
					z1[i] += e1imp*z1[mp];
					z2[i] += e1imp*z2[mp];
					z3[i] += e1imp*z3[mp];
				}
				ii = iphes[mp];
				if (ii != mp) {
					zsafe = z1[mp];
					z1[mp] = z1[ii];
					z1[ii] = zsafe;
					zsafe = z2[mp];
					z2[mp] = z2[ii];
					z2[ii] = zsafe;
					zsafe = z3[mp];
					z3[mp] = z3[ii];
					z3[ii] = zsafe;
				}
			}
			break;

		case (11):
			// mass = identity, Jacobian a full matrix, second order
		case (12):
			// ---  b = identity, Jacobian a banded matrix, second order
			for (int i = 0; i < n; i++) {
				s2 = -f2[i];
				s3 = -f3[i];
				z1[i] -= f1[i]*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			break;

		case (13):
			// mass is a banded matrix, Jacobian a full matrix, second order
		case (14):
			// mass is a banded matrix, Jacobian a banded matrix, second order
			for (int i = 0; i < m1; i++) {
				s2 = -f2[i];
				s3 = -f3[i];
				z1[i] -= f1[i]*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			for (int i = 0; i < nm1; i++) {
				s1 = s2 = s3 = 0.0;
				for (int j = max(0, i-mlmas); j < min(nm1, i+mumas+1); j++) {
					bb = fmas[i-j+mbdiag-1][j];
					s1 -= bb*f1[j+m1];
					s2 -= bb*f2[j+m1];
					s3 -= bb*f3[j+m1];
				}
				z1[i+m1] += s1*fac1;
				z2[i+m1] = z2[i+m1] + s2*alphn - s3*betan;
				z3[i+m1] = z3[i+m1] + s3*alphn + s2*betan;
			}
			break;

		case (15):
			// mass is a full matrix, Jacobian a full matrix, second order
			for (int i = 0; i < m1; i++) {
				s2 = -f2[i];
				s3 = -f3[i];
				z1[i] -= f1[i]*fac1;
				z2[i] = z2[i] + s2*alphn - s3*betan;
				z3[i] = z3[i] + s3*alphn + s2*betan;
			}
			for (int i = 0; i < nm1; i++) {
				s1 = s2 = s3 = 0.0;
				for (int j = 0; j < nm1; j++) {
					bb = fmas[i][j];
					s1 -= bb*f1[j+m1];
					s2 -= bb*f2[j+m1];
					s3 -= bb*f3[j+m1];
				}
				z1[i+m1] += s1*fac1;
				z2[i+m1] = z2[i+m1] + s2*alphn - s3*betan;
				z3[i+m1] = z3[i+m1] + s3*alphn + s2*betan;
			}
			break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);
	}

	switch (ijob) {
		case (1):
		case (2):
		case (3):
		case (4):
		case (5):
		case (7):
			break;

		case (11):
		case (13):
		case (15):
			abno = alphn*alphn + betan*betan;
			mm = m1/m2;
			for (int j = 0; j < m2; j++) {
				sum1 = sum2 = sum3 = 0.0;
				for (int k = mm - 1; k >= 0; k--) {
					jkm = j + k*m2;
					sum1 = (z1[jkm] + sum1)/fac1;
					sumh = (z2[jkm] + sum2)/abno;
					sum3 = (z3[jkm] + sum3)/abno;
					sum2 = sumh*alphn + sum3*betan;
					sum3 = sum3*alphn - sumh*betan;
					for (int i = 0; i < nm1; i++) {
						z1[i+m1] += fjac[i][jkm]*sum1;
						z2[i+m1] += fjac[i][jkm]*sum2;
						z3[i+m1] += fjac[i][jkm]*sum3;
					}
				}
			}
			sol(nm1, e1, &z1[m1], ip1);
			solc(nm1, e2r, e2i, &z2[m1], &z3[m1], ip2);
			break;

		case (12):
		case (14):
			abno = alphn*alphn + betan*betan;
			mm = m1/m2;
			for (int j = 0; j < m2; j++) {
				sum1 = sum2 = sum3 = 0.0;
				for (int k = mm - 1; k >= 0; k--) {
					jkm = j + k*m2;
					sum1 = (z1[jkm] + sum1)/fac1;
					sumh = (z2[jkm] + sum2)/abno;
					sum3 = (z3[jkm] + sum3)/abno;
					sum2 = sumh*alphn + sum3*betan;
					sum3 = sum3*alphn - sumh*betan;
					for (int i = max(0, j-mujac); i < min(nm1, j+mljac+1); i++) {
						ffja = fjac[i+mujac-j][jkm];
						z1[i+m1] += ffja*sum1;
						z2[i+m1] += ffja*sum2;
						z3[i+m1] += ffja*sum3;
					}
				}
			}
			solb(nm1, e1, mle, mue, &z1[m1], ip1);
			solbc(nm1, e2r, e2i, mle, mue, &z2[m1], &z3[m1], ip2);
			break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);
	}

	switch (ijob) {
		case (1):
		case (2):
		case (3):
		case (4):
		case (5):
		case (7):
			break;

		case (11):
		case (12):
		case (13):
		case (14):
		case (15):
		  	for (int i = m1 - 1; i >= 0; i--) {
				mpi = m2 + i;
				z1[i] = (z1[i] + z1[mpi])/fac1;
				z2i = z2[i] + z2[mpi];
				z3i = z3[i] + z3[mpi];
				z3[i] = (z3i*alphn - z2i*betan)/abno;
				z2[i] = (z2i*alphn + z3i*betan)/abno;
			}
		 	break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);

	}

	return (ier);

} // LinearSolve


int StiffIntegratorT::ErrorEstimate(double eps)
{

	int mm, ii, mp, ier = 0;
	double sum, zsafe;

	double hee1 = -(13.0 + 7.0*sqrt(6.0))/(3.0*h);
	double hee2 = (-13.0 + 7.0*sqrt(6.0))/(3.0*h);
	double hee3 = -1.0/(3.0*h);

	switch (ijob)
	{
		case (1):
			// mass = identity, Jacobian a full matrix
			for (int i = 0; i < n; i++) {
				f2[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
				cont[i] = f2[i] + y0[i];
			}
			sol(n, e1, cont, ip1);
			break;

		case (2):
			// mass = identity, Jacobian a banded matrix
			for (int i = 0; i < n; i++) {
				f2[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
				cont[i] = f2[i] + y0[i];
			}
			solb(n, e1, mle, mue, cont, ip1);
			break;

		case (3):
			// mass is a banded matrix, Jacobian a full matrix
			for (int i = 0; i < n; i++)
				f1[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
			for (int i = 0; i < n; i++) {
				sum = 0.0;
				for (int j = max(0, i-mlmas); j < min(n, i+mumas+1); j++)
					sum += fmas[i-j+mbdiag-1][j]*f1[j];
				f2[i] = sum;
				cont[i] = sum + y0[i];
			}
			sol(n, e1, cont, ip1);
			break;

		case (4):
			// mass is a banded matrix, Jacobian a banded matrix
			for (int i = 0; i < n; i++)
				f1[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
			for (int i = 0; i < n; i++) {
				sum = 0.0;
				for (int j = max(0, i-mlmas); j < min(n, i+mumas+1); j++)
					sum = sum + fmas[i-j+mbdiag-1][j]*f1[j];
				f2[i] = sum;
				cont[i] = sum + y0[i];
			}
			solb(n, e1, mle, mue, cont, ip1);
			break;

		case (5):
			// mass is a full matrix, Jacobian a full matrix
			for (int i = 0; i < n; i++)
				f1[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
			for (int i = 0; i < n; i++) {
				sum = 0.0;
				for (int j = 0; j < n; j++)
					sum += fmas[j][i]*f1[j];
				f2[i] = sum;
				cont[i] = sum + y0[i];
			}
			sol(n, e1, cont, ip1);
			break;

		case (6):
			// mass is a full matrix, Jacobian a banded matrix
			// this option is not provided
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);

		case (7):
			// mass = identity, Jacobian a full matrix, Hessenberg-option
			for (int i = 0; i < n; i++) {
				f2[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
				cont[i] = f2[i] + y0[i];
			}
			for (int mm1 = n - 3; mm1 >= 0; mm1--) {
				mp = n - mm1 - 2;
				ii = iphes[mp];
				if (ii != mp) {
					zsafe = cont[mp];
					cont[mp] = cont[ii];
					cont[ii] = zsafe;
				}
				for (int i = mp; i < n; i++)
					cont[i] -= fjac[i][mp-1]*cont[mp];
			}
			solh(n, e1, 1, cont, ip1);
			for (int mm1 = 0; mm1 < n - 2; mm1++) {
				mp = n - mm1 - 2;
				for (int i = mp; i < n; i++)
					cont[i] += fjac[i][mp-1]*cont[mp];
				ii = iphes[mp];
				if (ii != mp) {
					zsafe = cont[mp];
					cont[mp] = cont[ii];
					cont[ii] = zsafe;
				}
			}
			break;

		case (11):
			// mass = identity, Jacobian a full matrix, second order
			for (int i = 0; i < n; i++) {
				f2[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
				cont[i] = f2[i] + y0[i];
			}
			break;

		case (12):
			// mass = identity, Jacobian a banded matrix, second order
			for (int i = 0; i < n; i++) {
				f2[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
				cont[i] = f2[i] + y0[i];
			}
			break;

		case (13):
			// mass is a banded matrix, Jacobian a full matrix, second order
			for (int i = 0; i < m1; i++) {
				f1[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
				cont[i] = f2[i] + y0[i];
			}
			for (int i = m1; i < n; i++)
				f1[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
			for (int i = 0; i < nm1; i++) {
				sum = 0.0;
				for (int j = max(0, i-mlmas); j < min(nm1, i+mumas+1); j++)
					sum += fmas[i-j+mbdiag-1][j]*f1[j+m1];
				f2[i+m1] = sum;
				cont[i+m1] = sum + y0[i+m1];
			}
			break;

		case (14):
			// mass is a banded matrix, Jacobian a banded matrix, second order
			for (int i = 0; i < m1; i++) {
				f2[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
				cont[i] = f2[i] + y0[i];
			}
			for (int i = m1; i < n; i++)
				f1[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
			for (int i = 0; i < nm1; i++) {
				sum = 0.0;
				for (int j = max(0, i-mlmas); j < min(nm1, i+mumas+1); j++)
					sum += fmas[i-j+mbdiag-1][j]*f1[j+m1];
				f2[i+m1] = sum;
				cont[i+m1] = sum + y0[i+m1];
			}
			break;

		case (15):
			// mass is a banded matrix, Jacobian a full matrix, second order
			for (int i = 0; i < m1; i++) {
				f2[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
				cont[i] = f2[i] + y0[i];
			}
			for (int i = m1; i < n; i++)
				f1[i] = hee1*z1[i] + hee2*z2[i] + hee3*z3[i];
			for (int i = 0; i < nm1; i++) {
				sum = 0.0;
				for (int j = 0; j < nm1; j++)
					sum += fmas[j][i]*f1[j+m1];
				f2[i+m1] = sum;
				cont[i+m1] = sum + y0[i+m1];
			}
			break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);

	}

	switch (ijob)
	{
		case (1):
		case (2):
		case (3):
		case (4):
		case (5):
		case (7):
			break;

		case (11):
		case (13):
		case (15):
			mm = m1/m2;
			for (int j = 0; j < m2; j++) {
				sum = 0.0;
				for (int k = mm - 1; k >= 0; k--) {
					sum = (cont[j+k*m2] + sum)/fac1;
					for (int i = 0; i < nm1; i++)
						cont[i+m1] += fjac[i][j+k*m2]*sum;
				}
			}
			sol(nm1, e1, &cont[m1], ip1);
			for (int i = m1 - 1; i >= 0; i--)
				cont[i] = (cont[i] + cont[m2+i])/fac1;
			break;

		case (12):
		case (14):
			mm = m1/m2;
			for (int j = 0; j < m2; j++) {
				sum = 0.0;
				for (int k = mm - 1; k >= 0; k--) {
					sum = (cont[j+k*m2] + sum)/fac1;
					for (int i = max(0, j - mujac); i < min(nm1, j+mljac); i++)
						cont[i+m1] += fjac[i+mujac-j][j+k*m2]*sum;
				}
			}
			solb(nm1, e1, mle, mue, &cont[m1], ip1);
			for (int i = m1 - 1; i >= 0; i--)
					cont[i] = (cont[i] + cont[m2+i])/fac1;
			break;
		default:
			cout << "Not a value ijob. \tijob = " << ijob << endl;
			ier = -1;
			return (ier);

	}

	err = 0.0;
	for (int i = 0; i < n; i++)
		err += pow(cont[i]/scal[i], 2);
	err = max(sqrt(err/n), 1.0e-10);

	if (err < 1.0) return (ier);

	if (first || reject) {
		for (int i = 0; i < n; i++) cont[i] = y[i] + cont[i];
		Function(x, cont, f1,eps);
		nfcn++;
		for (int i = 0; i < n; i++) cont[i] = f1[i] + f2[i];

		switch (ijob)
		{
			case (1):
			case (3):
			case (5):
				// full matrix option
				sol(n, e1, cont, ip1);
				break;

			case (2):
			case (4):
				// banded matrix option
				solb(n, e1, mle, mue, cont, ip1);
				break;

			case (7):
				//Hessenberg matrix option
				// mass = identity, Jacobian a full matrix, Hessenberg-option
				for (int mm1 = n - 3; mm1 >= 0; mm1--) {
					mp = n - mm1 - 2;
					ii = iphes[mp];
					if (ii != mp) {
						zsafe = cont[mp];
						cont[mp] = cont[ii];
						cont[ii] = zsafe;
					}
					for (int i = mp; i < n; i++)
						cont[i] -= fjac[i][mp-1]*cont[mp];
				}
				solh(n, e1, 1, cont, ip1);
				for (int mm1 = 0; mm1 < n - 2; mm1++) {
					mp = n - mm1 - 2;
					for (int i = mp; i < n; i++)
						cont[i] += fjac[i][mp-1]*cont[mp];
					ii = iphes[mp];
					if (ii != mp) {
						zsafe = cont[mp];
						cont[mp] = cont[ii];
						cont[ii] = zsafe;
					}
				}
				break;

			case (11):
			case (13):
			case (15):
				// Full matrix option, second order
				for (int j = 0; j < m2; j++) {
					sum = 0.0;
					for (int k = mm - 1; k >= 0; k--) {
						sum = (cont[j+k*m2] + sum)/fac1;
						for (int i = 0; i < nm1; i++)
							cont[i+m1] += fjac[i][j+k*m2]*sum;
					}
				}
				sol(nm1, e1, &cont[m1], ip1);
				for (int i = m1 - 1; i >= 0; i--)
					cont[i] = (cont[i] + cont[m2+i])/fac1;
				break;

			case (12):
			case (14):
				// Banded matrix option, second order
				for (int j = 0; j < m2; j++) {
					sum = 0.0;
					for (int k = mm - 1; k >= 0; k--) {
						sum = (cont[j+k*m2] + sum)/fac1;
						for (int i = max(0, j-mujac); i < min(nm1, j+mljac); i++)
							cont[i+m1] += fjac[i+mujac-j][j+k*m2]*sum;
					}
				}
				solb(nm1, e1, mle, mue, &cont[m1], ip1);
				for (int i = m1 - 1; i >= 0; i--)
					cont[i] = (cont[i] + cont[m2+i])/fac1;
				break;
			default:
				cout << "Not a value ijob. \tijob = " << ijob << endl;
				ier = -1;
				return (ier);
		}

		err = 0.0;
		for (int i = 0; i < n; i++)
			err += pow(cont[i]/scal[i], 2);
		err = max(sqrt(err/n), 1.0e-10);
	}

	return (ier);

} // ErrorEstimate

