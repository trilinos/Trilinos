/* Try to test all combinations of types and operations */

#ifdef SACADO_NAMESPACE
#define ADT_RAD Sacado::Rad::
#else
#define ADT_RAD /*nothing*/
#endif

#include "Sacado_trad.hpp"
#include <stdio.h>

typedef ADT_RAD IndepADvar<double> AI;
typedef ADT_RAD ADvar<double> A;
typedef ADT_RAD ConstADvar<double> C;
typedef ADT_RAD ADvari<double> Ai;
typedef const ADT_RAD IndepADvar<double> cAI;
typedef const ADT_RAD ADvar<double> cA;
typedef const ADT_RAD ConstADvar<double> cC;
typedef const ADT_RAD ADvari<double> cAi;
static int rc;

/* This is to be run through an awk program that changes lines */
/* with "BINTEST" or "UNOPTEST" at the beginning of the line into */
/* a the desired C++ (which we can then inspect). */

 void
botch(const char *what, double wanted, double got)
{
	printf("%s: expected %g, got %g, diff = %.2g\n", what, wanted, got, wanted-got);
	rc = 1;
	}

 const double tol = 5e-16;

 int
differ(double a, double b)
{
	double d = a - b;
	if (d < 0.)
		d = -d;
	if (a < 0.)
		a = -a;
	if (b < 0.)
		b = -b;
	if (a < b)
		a = b;
	if (a > 0.)
		d /= a;
	return d > tol;
	}

#ifdef RAD_NO_EQ_ALIAS
#define Plus_dx 1.
#else
#ifdef RAD_AUTO_AD_Const
#define Plus_dx 1.
#else
#define Plus_dx 0.
#endif
#endif

 int
main(void)
{
	AI xAI, yAI;
	A fA, xA, yA;
	C xC, yC;
	double dx, dy, f, xd, yd;
	long xL, yL;
	int xi, yi;

	rc = 0;

	/**** Test of operator+ ****/

	xd = 28.; f = 28.; dx = Plus_dx;
	xAI = xd;
	fA = operator+(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator+(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = operator+(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator+(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = operator+(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator+(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = operator+(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator+(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of operator+ ****/

	xd = 28.; f = 28.; dx = 1.;
	xAI = xd;
	fA = operator+(copy(xAI));
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(copy(xAI))", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(copy(xAI))/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator+(copy(xcAI));
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(copy(xcAI))", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(copy(xcAI))/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = operator+(copy(xA));
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(copy(xA))", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(copy(xA))/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator+(copy(xcA));
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(copy(xcA))", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(copy(xcA))/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = operator+(copy(xC));
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(copy(xC))", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(copy(xC))/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator+(copy(xcC));
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(copy(xcC))", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(copy(xcC))/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = operator+(copy(xcAi));
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(copy(xcAi))", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(copy(xcAi))/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator+(copy(xcAi));
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(copy(xcAi))", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(copy(xcAi))/dx", dx, xcAi.aval);
	}


	/**** Test of operator- ****/

	xd = 42.; f = -42.; dx = -1.;
	xAI = xd;
	fA = operator-(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator-(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = operator-(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator-(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = operator-(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator-(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = operator-(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator-(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of acos ****/

	xd = .7; f = acos(.7); dx = -1.4002800840280099;
	xAI = xd;
	fA = acos(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acos(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d acos(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = acos(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acos(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d acos(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = acos(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acos(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d acos(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = acos(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acos(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d acos(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = acos(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acos(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d acos(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = acos(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acos(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d acos(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = acos(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acos(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d acos(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = acos(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acos(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d acos(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of acosh ****/

	xd = 1.25; f = acosh(1.25); dx = 1.3333333333333333;
	xAI = xd;
	fA = acosh(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acosh(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d acosh(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = acosh(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acosh(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d acosh(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = acosh(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acosh(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d acosh(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = acosh(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acosh(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d acosh(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = acosh(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acosh(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d acosh(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = acosh(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acosh(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d acosh(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = acosh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acosh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d acosh(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = acosh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = acosh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d acosh(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of asin ****/

	xd = .7; f = asin(.7); dx = 1.4002800840280099;
	xAI = xd;
	fA = asin(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asin(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d asin(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = asin(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asin(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d asin(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = asin(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asin(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d asin(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = asin(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asin(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d asin(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = asin(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asin(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d asin(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = asin(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asin(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d asin(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = asin(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asin(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d asin(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = asin(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asin(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d asin(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of asinh ****/

	xd = -2.; f = asinh(-2.); dx = 0.4472135954999579;
	xAI = xd;
	fA = asinh(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asinh(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d asinh(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = asinh(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asinh(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d asinh(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = asinh(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asinh(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d asinh(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = asinh(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asinh(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d asinh(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = asinh(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asinh(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d asinh(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = asinh(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asinh(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d asinh(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = asinh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asinh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d asinh(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = asinh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = asinh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d asinh(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of atan ****/

	xd = 2.; f = atan(2.); dx = 0.2;
	xAI = xd;
	fA = atan(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = atan(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = atan(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = atan(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = atan(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = atan(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = atan(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = atan(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of atanh ****/

	xd = .6; f = atanh(.6); dx = 1.5625;
	xAI = xd;
	fA = atanh(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atanh(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atanh(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = atanh(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atanh(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atanh(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = atanh(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atanh(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atanh(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = atanh(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atanh(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atanh(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = atanh(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atanh(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atanh(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = atanh(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atanh(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atanh(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = atanh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atanh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atanh(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = atanh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atanh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atanh(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of cos ****/

	xd = 3.; f = cos(3.); dx = -sin(3.);
	xAI = xd;
	fA = cos(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cos(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d cos(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = cos(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cos(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d cos(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = cos(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cos(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d cos(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = cos(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cos(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d cos(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = cos(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cos(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d cos(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = cos(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cos(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d cos(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = cos(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cos(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d cos(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = cos(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cos(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d cos(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of cosh ****/

	xd = 2.; f = cosh(2.); dx = sinh(2.);
	xAI = xd;
	fA = cosh(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cosh(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d cosh(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = cosh(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cosh(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d cosh(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = cosh(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cosh(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d cosh(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = cosh(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cosh(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d cosh(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = cosh(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cosh(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d cosh(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = cosh(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cosh(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d cosh(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = cosh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cosh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d cosh(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = cosh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = cosh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d cosh(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of exp ****/

	xd = 2.; f = exp(2.); dx = exp(2.);
	xAI = xd;
	fA = exp(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = exp(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d exp(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = exp(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = exp(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d exp(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = exp(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = exp(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d exp(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = exp(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = exp(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d exp(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = exp(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = exp(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d exp(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = exp(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = exp(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d exp(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = exp(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = exp(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d exp(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = exp(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = exp(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d exp(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of log ****/

	xd = 4.; f = log(4.); dx = 0.25;
	xAI = xd;
	fA = log(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d log(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = log(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d log(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = log(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d log(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = log(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d log(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = log(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d log(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = log(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d log(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = log(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d log(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = log(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d log(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of log10 ****/

	xd = 100.; f = 2.; dx = .01/log(10.);
	xAI = xd;
	fA = log10(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log10(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d log10(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = log10(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log10(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d log10(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = log10(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log10(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d log10(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = log10(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log10(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d log10(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = log10(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log10(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d log10(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = log10(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log10(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d log10(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = log10(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log10(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d log10(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = log10(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = log10(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d log10(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of sin ****/

	xd = 5.; f = sin(5.); dx = cos(5.);
	xAI = xd;
	fA = sin(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sin(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d sin(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = sin(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sin(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d sin(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = sin(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sin(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d sin(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = sin(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sin(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d sin(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = sin(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sin(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d sin(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = sin(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sin(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d sin(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = sin(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sin(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d sin(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = sin(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sin(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d sin(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of sinh ****/

	xd = 3.; f = sinh(3.); dx = cosh(3.);
	xAI = xd;
	fA = sinh(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sinh(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d sinh(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = sinh(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sinh(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d sinh(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = sinh(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sinh(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d sinh(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = sinh(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sinh(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d sinh(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = sinh(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sinh(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d sinh(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = sinh(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sinh(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d sinh(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = sinh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sinh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d sinh(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = sinh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sinh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d sinh(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of sqrt ****/

	xd = 25.; f = 5.; dx = 0.1;
	xAI = xd;
	fA = sqrt(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sqrt(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d sqrt(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = sqrt(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sqrt(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d sqrt(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = sqrt(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sqrt(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d sqrt(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = sqrt(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sqrt(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d sqrt(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = sqrt(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sqrt(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d sqrt(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = sqrt(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sqrt(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d sqrt(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = sqrt(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sqrt(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d sqrt(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = sqrt(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = sqrt(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d sqrt(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of tan ****/

	xd = 3.; f = tan(3.); dx = 1./(cos(3.)*cos(3.));
	xAI = xd;
	fA = tan(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tan(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d tan(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = tan(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tan(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d tan(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = tan(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tan(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d tan(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = tan(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tan(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d tan(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = tan(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tan(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d tan(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = tan(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tan(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d tan(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = tan(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tan(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d tan(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = tan(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tan(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d tan(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of tanh ****/

	xd = 1.; f = tanh(1.); dx = (1./cosh(1.))*(1./cosh(1.));
	xAI = xd;
	fA = tanh(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tanh(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d tanh(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = tanh(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tanh(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d tanh(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = tanh(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tanh(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d tanh(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = tanh(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tanh(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d tanh(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = tanh(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tanh(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d tanh(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = tanh(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tanh(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d tanh(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = tanh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tanh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d tanh(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = tanh(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = tanh(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d tanh(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of fabs ****/

	xd = -12.; f = 12.; dx = -1.;
	xAI = xd;
	fA = fabs(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d fabs(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = fabs(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d fabs(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = fabs(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d fabs(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = fabs(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d fabs(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = fabs(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d fabs(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = fabs(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d fabs(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = fabs(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d fabs(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = fabs(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d fabs(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of fabs ****/

	xd = 37.; f = 37.; dx = 1.;
	xAI = xd;
	fA = fabs(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d fabs(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = fabs(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d fabs(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = fabs(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d fabs(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = fabs(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d fabs(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = fabs(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d fabs(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = fabs(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d fabs(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = fabs(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d fabs(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = fabs(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = fabs(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d fabs(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of abs ****/

	xd = -12.; f = 12.; dx = -1.;
	xAI = xd;
	fA = abs(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d abs(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = abs(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d abs(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = abs(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d abs(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = abs(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d abs(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = abs(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d abs(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = abs(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d abs(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = abs(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d abs(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = abs(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d abs(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of abs ****/

	xd = 37.; f = 37.; dx = 1.;
	xAI = xd;
	fA = abs(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d abs(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = abs(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d abs(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = abs(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d abs(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = abs(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d abs(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = abs(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d abs(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = abs(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d abs(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = abs(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d abs(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = abs(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = abs(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d abs(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of + ****/

	xd = 17.; f = 17.; dx = Plus_dx;
	xAI = xd;
	fA = +xAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +xAI", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d +xAI/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = +xcAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +xcAI", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d +xcAI/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = +xA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +xA", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d +xA/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = +xcA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +xcA", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d +xcA/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = +xC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +xC", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d +xC/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = +xcC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +xcC", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d +xcC/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = +xcAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +xcAi", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d +xcAi/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = +xcAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +xcAi", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d +xcAi/dx", dx, xcAi.aval);
	}


	/**** Test of + ****/

	xd = 17.; f = 17.; dx = 1.;
	xAI = xd;
	fA = +copy(xAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +copy(xAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d +copy(xAI)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = +copy(xcAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +copy(xcAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d +copy(xcAI)/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = +copy(xA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +copy(xA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d +copy(xA)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = +copy(xcA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +copy(xcA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d +copy(xcA)/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = +copy(xC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +copy(xC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d +copy(xC)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = +copy(xcC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +copy(xcC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d +copy(xcC)/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = +copy(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +copy(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d +copy(xcAi)/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = +copy(xcAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = +copy(xcAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d +copy(xcAi)/dx", dx, xcAi.aval);
	}


	/**** Test of - ****/

	xd = 19.; f = -19.; dx = -1.;
	xAI = xd;
	fA = -xAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = -xAI", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d -xAI/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = -xcAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = -xcAI", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d -xcAI/dx", dx, xcAI.adj());
	}

	xA = xd;
	fA = -xA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = -xA", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d -xA/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = -xcA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = -xcA", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d -xcA/dx", dx, xcA.adj());
	}

	xC = xd;
	fA = -xC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = -xC", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d -xC/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = -xcC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = -xcC", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d -xcC/dx", dx, xcC.adj());
	}

	{
	cAi xcAi(xd);
	fA = -xcAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = -xcAi", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d -xcAi/dx", dx, xcAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = -xcAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = -xcAi", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d -xcAi/dx", dx, xcAi.aval);
	}


	if (!rc) // chatter for cppunit test, which cannot tolerate silence
		printf("OK\n");
	return rc;
	}
