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

	xd = 1.; yd = 2.; f = 3.; dx = 1.; dy = 1.;
	xAI = xd;
	yAI = yd;
	fA = operator+(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = operator+(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = operator+(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = operator+(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = operator+(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator+(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = operator+(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator+(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = operator+(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = operator+(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = operator+(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator+(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = operator+(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator+(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = operator+(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = operator+(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = operator+(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator+(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = operator+(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator+(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = operator+(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = operator+(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = operator+(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator+(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = operator+(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = operator+(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = operator+(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator+(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = operator+(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator+(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = operator+(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = operator+(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = operator+(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = operator+(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = operator+(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator+(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = operator+(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator+(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = operator+(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = operator+(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = operator+(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator+(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = operator+(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator+(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = operator+(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = operator+(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = operator+(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator+(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = operator+(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator+(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = operator+(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = operator+(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = operator+(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator+(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = operator+(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = operator+(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = operator+(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator+(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = operator+(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator+(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = operator+(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = operator+(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = operator+(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = operator+(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = operator+(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator+(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = operator+(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator+(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = operator+(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = operator+(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = operator+(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator+(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = operator+(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator+(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = operator+(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = operator+(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = operator+(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator+(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = operator+(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator+(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = operator+(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = operator+(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = operator+(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator+(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = operator+(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = operator+(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = operator+(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator+(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = operator+(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator+(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = operator+(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator+(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = operator+(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator+(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = operator+(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = operator+(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator+(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = operator+(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator+(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = operator+(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d operator+(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = operator+(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator+(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = operator+(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator+(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = operator+(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d operator+(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = operator+(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator+(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = operator+(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator+(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = operator+(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d operator+(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = operator+(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator+(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = operator+(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = operator+(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = operator+(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator+(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = operator+(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator+(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = operator+(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = operator+(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = operator+(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator+(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = operator+(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = operator+(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator+(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = operator+(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = operator+(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator+(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = operator+(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = operator+(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = operator+(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = operator+(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator+(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = operator+(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = operator+(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator+(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = operator+(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = operator+(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator+(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = operator+(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = operator+(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator+(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = operator+(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator+(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = operator+(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator+(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = operator+(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator+(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = operator+(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator+(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = operator+(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator+(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = operator+(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator+(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = operator+(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator+(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator+(xi,ycAi)/dy", dy, ycAi.aval);
	}


	/**** Test of operator- ****/

	xd = 7.; yd = 4; f = 3.; dx = 1.; dy = -1.;
	xAI = xd;
	yAI = yd;
	fA = operator-(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = operator-(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = operator-(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = operator-(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = operator-(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator-(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = operator-(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator-(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = operator-(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = operator-(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = operator-(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator-(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = operator-(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator-(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = operator-(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = operator-(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = operator-(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator-(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = operator-(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator-(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = operator-(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = operator-(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = operator-(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator-(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = operator-(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = operator-(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = operator-(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator-(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = operator-(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator-(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = operator-(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = operator-(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = operator-(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = operator-(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = operator-(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator-(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = operator-(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator-(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = operator-(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = operator-(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = operator-(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator-(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = operator-(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator-(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = operator-(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = operator-(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = operator-(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator-(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = operator-(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator-(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = operator-(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = operator-(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = operator-(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator-(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = operator-(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = operator-(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = operator-(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator-(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = operator-(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator-(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = operator-(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = operator-(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = operator-(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = operator-(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = operator-(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator-(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = operator-(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator-(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = operator-(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = operator-(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = operator-(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator-(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = operator-(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator-(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = operator-(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = operator-(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = operator-(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator-(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = operator-(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator-(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = operator-(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = operator-(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = operator-(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator-(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = operator-(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = operator-(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = operator-(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator-(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = operator-(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator-(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = operator-(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator-(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = operator-(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator-(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = operator-(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = operator-(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator-(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = operator-(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator-(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = operator-(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d operator-(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = operator-(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator-(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = operator-(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator-(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = operator-(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d operator-(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = operator-(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator-(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = operator-(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator-(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = operator-(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d operator-(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = operator-(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator-(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = operator-(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = operator-(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = operator-(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator-(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = operator-(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator-(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = operator-(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = operator-(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = operator-(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator-(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = operator-(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = operator-(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator-(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = operator-(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = operator-(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator-(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = operator-(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = operator-(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = operator-(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = operator-(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator-(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = operator-(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = operator-(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator-(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = operator-(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = operator-(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator-(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = operator-(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = operator-(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator-(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = operator-(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator-(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = operator-(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator-(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = operator-(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator-(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = operator-(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator-(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = operator-(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator-(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = operator-(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator-(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = operator-(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator-(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator-(xi,ycAi)/dy", dy, ycAi.aval);
	}


	/**** Test of operator* ****/

	xd = 6.; yd = 7.; f = 42.; dx = 7.; dy = 6.;
	xAI = xd;
	yAI = yd;
	fA = operator*(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = operator*(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = operator*(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = operator*(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = operator*(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator*(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = operator*(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator*(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = operator*(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = operator*(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = operator*(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator*(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = operator*(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator*(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = operator*(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = operator*(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = operator*(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator*(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = operator*(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator*(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = operator*(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = operator*(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = operator*(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator*(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = operator*(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = operator*(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = operator*(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator*(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = operator*(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator*(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = operator*(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = operator*(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = operator*(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = operator*(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = operator*(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator*(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = operator*(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator*(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = operator*(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = operator*(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = operator*(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator*(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = operator*(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator*(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = operator*(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = operator*(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = operator*(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator*(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = operator*(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator*(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = operator*(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = operator*(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = operator*(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator*(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = operator*(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = operator*(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = operator*(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator*(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = operator*(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator*(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = operator*(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = operator*(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = operator*(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = operator*(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = operator*(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator*(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = operator*(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator*(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = operator*(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = operator*(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = operator*(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator*(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = operator*(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator*(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = operator*(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = operator*(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = operator*(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator*(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = operator*(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator*(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = operator*(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = operator*(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = operator*(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator*(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = operator*(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = operator*(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = operator*(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator*(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = operator*(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator*(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = operator*(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator*(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = operator*(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator*(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator*(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = operator*(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = operator*(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator*(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = operator*(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator*(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator*(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = operator*(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d operator*(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = operator*(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator*(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = operator*(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator*(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator*(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = operator*(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d operator*(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = operator*(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator*(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = operator*(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator*(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator*(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = operator*(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d operator*(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = operator*(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator*(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator*(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = operator*(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = operator*(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator*(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = operator*(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator*(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = operator*(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator*(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = operator*(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = operator*(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = operator*(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator*(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = operator*(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = operator*(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator*(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = operator*(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = operator*(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator*(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = operator*(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = operator*(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = operator*(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = operator*(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator*(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = operator*(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = operator*(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator*(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = operator*(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = operator*(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator*(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = operator*(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = operator*(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator*(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = operator*(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator*(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = operator*(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator*(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = operator*(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator*(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = operator*(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator*(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = operator*(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator*(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = operator*(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator*(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = operator*(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator*(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator*(xi,ycAi)/dy", dy, ycAi.aval);
	}


	/**** Test of operator/ ****/

	xd = 24.; yd = 6.; f = 4.; dx = 1/6.; dy = -2./3.;
	xAI = xd;
	yAI = yd;
	fA = operator/(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = operator/(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = operator/(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = operator/(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = operator/(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator/(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = operator/(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator/(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = operator/(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = operator/(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = operator/(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator/(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = operator/(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator/(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = operator/(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = operator/(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = operator/(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator/(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = operator/(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator/(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = operator/(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = operator/(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = operator/(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator/(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = operator/(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = operator/(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = operator/(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator/(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = operator/(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator/(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = operator/(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = operator/(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = operator/(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = operator/(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = operator/(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator/(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = operator/(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator/(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = operator/(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = operator/(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = operator/(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator/(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = operator/(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator/(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = operator/(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = operator/(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = operator/(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator/(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = operator/(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator/(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = operator/(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = operator/(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = operator/(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator/(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = operator/(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = operator/(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = operator/(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator/(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = operator/(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator/(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = operator/(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = operator/(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = operator/(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = operator/(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = operator/(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator/(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = operator/(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator/(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = operator/(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = operator/(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = operator/(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator/(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = operator/(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator/(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = operator/(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = operator/(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = operator/(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator/(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = operator/(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator/(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = operator/(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = operator/(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = operator/(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator/(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = operator/(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = operator/(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = operator/(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator/(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = operator/(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator/(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = operator/(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator/(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = operator/(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator/(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator/(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = operator/(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = operator/(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator/(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = operator/(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator/(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator/(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = operator/(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d operator/(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = operator/(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator/(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = operator/(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator/(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator/(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = operator/(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d operator/(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = operator/(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator/(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = operator/(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator/(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator/(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = operator/(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d operator/(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = operator/(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator/(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator/(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = operator/(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = operator/(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator/(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = operator/(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator/(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = operator/(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator/(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = operator/(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = operator/(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = operator/(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator/(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = operator/(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = operator/(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator/(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = operator/(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = operator/(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator/(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = operator/(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = operator/(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = operator/(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = operator/(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator/(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = operator/(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = operator/(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator/(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = operator/(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = operator/(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator/(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = operator/(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = operator/(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator/(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = operator/(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator/(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = operator/(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator/(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = operator/(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator/(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = operator/(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator/(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = operator/(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator/(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = operator/(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator/(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = operator/(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator/(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator/(xi,ycAi)/dy", dy, ycAi.aval);
	}


	/**** Test of atan2 ****/

	xd = 3.; yd = 4.; f = atan2(3.,4.); dx = 0.16; dy = -0.12;
	xAI = xd;
	yAI = yd;
	fA = atan2(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = atan2(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = atan2(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = atan2(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = atan2(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d atan2(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = atan2(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d atan2(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = atan2(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = atan2(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = atan2(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d atan2(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = atan2(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d atan2(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = atan2(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = atan2(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = atan2(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d atan2(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = atan2(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d atan2(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = atan2(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = atan2(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = atan2(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = atan2(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = atan2(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = atan2(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = atan2(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d atan2(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = atan2(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d atan2(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = atan2(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = atan2(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = atan2(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = atan2(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = atan2(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d atan2(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = atan2(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d atan2(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = atan2(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = atan2(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = atan2(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d atan2(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = atan2(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d atan2(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = atan2(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = atan2(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = atan2(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d atan2(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = atan2(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d atan2(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = atan2(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = atan2(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = atan2(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = atan2(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = atan2(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = atan2(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = atan2(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d atan2(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = atan2(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d atan2(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = atan2(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = atan2(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = atan2(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = atan2(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = atan2(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d atan2(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = atan2(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d atan2(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = atan2(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = atan2(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = atan2(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d atan2(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = atan2(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d atan2(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = atan2(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = atan2(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = atan2(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d atan2(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = atan2(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d atan2(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = atan2(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = atan2(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = atan2(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = atan2(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = atan2(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = atan2(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = atan2(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d atan2(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = atan2(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d atan2(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = atan2(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d atan2(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = atan2(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan2(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d atan2(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = atan2(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = atan2(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d atan2(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = atan2(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan2(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d atan2(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = atan2(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d atan2(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = atan2(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d atan2(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = atan2(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan2(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d atan2(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = atan2(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d atan2(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = atan2(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d atan2(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = atan2(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan2(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d atan2(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = atan2(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d atan2(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = atan2(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = atan2(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan2(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = atan2(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = atan2(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan2(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = atan2(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d atan2(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = atan2(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d atan2(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = atan2(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = atan2(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = atan2(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d atan2(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = atan2(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = atan2(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d atan2(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = atan2(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = atan2(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d atan2(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = atan2(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = atan2(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = atan2(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = atan2(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d atan2(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = atan2(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = atan2(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d atan2(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = atan2(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = atan2(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d atan2(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = atan2(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = atan2(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d atan2(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = atan2(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d atan2(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = atan2(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d atan2(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = atan2(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d atan2(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = atan2(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d atan2(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = atan2(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d atan2(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = atan2(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d atan2(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = atan2(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = atan2(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d atan2(xi,ycAi)/dy", dy, ycAi.aval);
	}


	/**** Test of pow ****/

	xd = 2.; yd = 3.; f = 8.; dx = 12.; dy = 8.*log(2.);
	xAI = xd;
	yAI = yd;
	fA = pow(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d pow(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = pow(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d pow(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = pow(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = pow(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = pow(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d pow(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = pow(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d pow(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = pow(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d pow(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = pow(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d pow(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = pow(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d pow(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = pow(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d pow(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = pow(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d pow(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = pow(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d pow(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = pow(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d pow(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = pow(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d pow(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = pow(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d pow(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = pow(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d pow(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = pow(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = pow(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = pow(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = pow(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = pow(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d pow(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = pow(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d pow(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = pow(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d pow(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = pow(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d pow(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = pow(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = pow(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = pow(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d pow(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = pow(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d pow(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = pow(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d pow(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = pow(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d pow(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = pow(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d pow(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = pow(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d pow(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = pow(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d pow(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = pow(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d pow(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = pow(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d pow(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = pow(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d pow(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = pow(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d pow(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = pow(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d pow(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = pow(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = pow(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = pow(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = pow(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = pow(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d pow(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = pow(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d pow(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = pow(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d pow(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = pow(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d pow(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = pow(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = pow(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = pow(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d pow(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = pow(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d pow(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = pow(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d pow(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = pow(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d pow(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = pow(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d pow(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = pow(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d pow(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = pow(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d pow(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = pow(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d pow(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = pow(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d pow(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = pow(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d pow(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = pow(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d pow(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = pow(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d pow(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = pow(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = pow(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = pow(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = pow(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = pow(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d pow(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = pow(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d pow(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = pow(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d pow(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = pow(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d pow(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d pow(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = pow(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d pow(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = pow(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d pow(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = pow(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d pow(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d pow(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = pow(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d pow(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = pow(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d pow(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = pow(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d pow(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d pow(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = pow(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d pow(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = pow(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d pow(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = pow(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d pow(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d pow(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = pow(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d pow(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = pow(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = pow(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d pow(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = pow(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = pow(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d pow(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = pow(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d pow(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = pow(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d pow(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = pow(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d pow(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = pow(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = pow(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d pow(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = pow(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d pow(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = pow(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d pow(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = pow(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d pow(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = pow(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d pow(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = pow(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d pow(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = pow(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d pow(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = pow(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = pow(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d pow(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = pow(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d pow(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = pow(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d pow(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = pow(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d pow(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = pow(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d pow(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = pow(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d pow(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = pow(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d pow(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = pow(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d pow(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = pow(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d pow(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = pow(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d pow(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = pow(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d pow(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = pow(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d pow(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = pow(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d pow(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = pow(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = pow(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d pow(xi,ycAi)/dy", dy, ycAi.aval);
	}


	/**** Test of operator< ****/

	xd = 2.; yd = 3.; f = 1.; dx = 0.; dy = 0.;
	xAI = xd;
	yAI = yd;
	fA = operator<(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = operator<(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = operator<(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = operator<(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = operator<(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = operator<(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = operator<(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = operator<(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = operator<(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = operator<(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = operator<(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = operator<(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = operator<(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = operator<(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = operator<(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = operator<(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = operator<(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator<(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = operator<(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = operator<(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = operator<(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = operator<(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = operator<(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = operator<(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = operator<(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = operator<(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = operator<(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = operator<(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = operator<(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = operator<(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = operator<(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = operator<(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = operator<(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = operator<(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = operator<(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = operator<(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = operator<(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = operator<(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = operator<(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator<(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = operator<(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = operator<(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = operator<(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = operator<(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = operator<(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = operator<(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = operator<(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = operator<(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = operator<(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = operator<(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = operator<(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = operator<(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = operator<(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = operator<(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = operator<(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = operator<(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = operator<(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = operator<(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = operator<(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = operator<(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = operator<(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator<(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = operator<(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = operator<(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = operator<(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = operator<(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = operator<(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator<(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = operator<(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator<(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = operator<(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = operator<(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator<(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = operator<(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator<(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = operator<(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d operator<(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = operator<(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator<(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = operator<(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator<(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = operator<(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d operator<(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = operator<(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator<(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = operator<(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator<(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = operator<(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d operator<(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = operator<(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator<(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = operator<(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = operator<(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = operator<(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = operator<(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = operator<(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = operator<(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = operator<(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = operator<(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = operator<(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = operator<(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = operator<(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = operator<(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = operator<(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = operator<(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = operator<(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = operator<(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = operator<(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = operator<(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = operator<(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = operator<(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = operator<(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = operator<(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = operator<(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = operator<(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = operator<(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = operator<(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = operator<(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = operator<(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xi,ycAi)/dy", dy, ycAi.aval);
	}


	/**** Test of operator< ****/

	xd = 3.; yd = 3.; f = 0.; dx = 0.; dy = 0.;
	xAI = xd;
	yAI = yd;
	fA = operator<(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = operator<(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = operator<(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = operator<(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = operator<(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = operator<(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = operator<(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = operator<(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = operator<(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = operator<(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = operator<(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = operator<(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = operator<(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = operator<(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = operator<(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = operator<(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = operator<(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator<(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = operator<(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = operator<(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = operator<(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = operator<(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = operator<(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = operator<(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = operator<(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = operator<(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = operator<(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = operator<(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = operator<(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = operator<(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = operator<(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = operator<(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = operator<(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = operator<(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = operator<(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = operator<(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = operator<(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = operator<(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = operator<(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator<(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = operator<(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = operator<(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = operator<(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = operator<(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = operator<(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = operator<(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = operator<(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = operator<(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = operator<(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = operator<(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = operator<(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = operator<(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = operator<(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = operator<(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = operator<(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = operator<(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = operator<(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = operator<(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = operator<(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = operator<(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = operator<(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator<(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = operator<(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = operator<(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = operator<(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = operator<(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = operator<(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator<(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = operator<(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator<(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = operator<(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = operator<(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator<(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = operator<(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator<(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = operator<(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d operator<(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = operator<(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator<(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = operator<(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator<(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = operator<(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d operator<(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = operator<(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator<(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = operator<(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator<(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = operator<(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d operator<(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = operator<(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator<(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = operator<(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = operator<(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = operator<(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = operator<(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = operator<(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = operator<(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = operator<(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = operator<(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = operator<(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = operator<(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = operator<(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = operator<(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = operator<(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = operator<(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = operator<(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = operator<(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = operator<(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = operator<(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = operator<(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = operator<(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = operator<(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = operator<(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = operator<(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = operator<(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = operator<(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = operator<(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = operator<(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = operator<(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<(xi,ycAi)/dy", dy, ycAi.aval);
	}


	/**** Test of operator<= ****/

	xd = 2.; yd = 3.; f = 1.; dx = 0.; dy = 0.;
	xAI = xd;
	yAI = yd;
	fA = operator<=(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xAI,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = operator<=(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xcAI,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = operator<=(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xAI,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = operator<=(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xcAI,ycAI)/dy", dy, ycAI.adj());
	}

	xAI = xd;
	yA = yd;
	fA = operator<=(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xAI,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = operator<=(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xcAI,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = operator<=(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xAI,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = operator<=(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xcAI,ycA)/dy", dy, ycA.adj());
	}

	xAI = xd;
	yC = yd;
	fA = operator<=(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xAI,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = operator<=(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xcAI,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = operator<=(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xAI,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = operator<=(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xcAI,ycC)/dy", dy, ycC.adj());
	}

	{
	xAI = xd;
	Ai yAi(yd);
	fA = operator<=(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = operator<=(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xcAI,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = operator<=(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xAI,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = operator<=(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xcAI,ycAi)/dy", dy, ycAi.aval);
	}

	xAI = xd;
	fA = operator<=(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,yd)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator<=(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,yd)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yL = (long)yd;
	fA = operator<=(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,yL)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = operator<=(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,yL)/dx", dx, xcAI.adj());
	}

	xAI = xd;
	yi = (int)yd;
	fA = operator<=(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator<=(xAI,yi)/dx", dx, xAI.adj());

	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = operator<=(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator<=(xcAI,yi)/dx", dx, xcAI.adj());
	}

	xA = xd;
	yAI = yd;
	fA = operator<=(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xA,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = operator<=(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xcA,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = operator<=(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xA,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = operator<=(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xcA,ycAI)/dy", dy, ycAI.adj());
	}

	xA = xd;
	yA = yd;
	fA = operator<=(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xA,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = operator<=(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xcA,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = operator<=(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xA,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = operator<=(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xcA,ycA)/dy", dy, ycA.adj());
	}

	xA = xd;
	yC = yd;
	fA = operator<=(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xA,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = operator<=(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xcA,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = operator<=(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xA,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = operator<=(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xcA,ycC)/dy", dy, ycC.adj());
	}

	{
	xA = xd;
	Ai yAi(yd);
	fA = operator<=(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = operator<=(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xcA,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = operator<=(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xA,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = operator<=(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xcA,ycAi)/dy", dy, ycAi.aval);
	}

	xA = xd;
	fA = operator<=(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,yd)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator<=(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,yd)/dx", dx, xcA.adj());
	}

	xA = xd;
	yL = (long)yd;
	fA = operator<=(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,yL)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = operator<=(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,yL)/dx", dx, xcA.adj());
	}

	xA = xd;
	yi = (int)yd;
	fA = operator<=(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator<=(xA,yi)/dx", dx, xA.adj());

	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = operator<=(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator<=(xcA,yi)/dx", dx, xcA.adj());
	}

	xC = xd;
	yAI = yd;
	fA = operator<=(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xC,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = operator<=(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xcC,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = operator<=(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xC,ycAI)/dy", dy, ycAI.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = operator<=(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xcC,ycAI)/dy", dy, ycAI.adj());
	}

	xC = xd;
	yA = yd;
	fA = operator<=(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xC,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = operator<=(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xcC,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = operator<=(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xC,ycA)/dy", dy, ycA.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = operator<=(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xcC,ycA)/dy", dy, ycA.adj());
	}

	xC = xd;
	yC = yd;
	fA = operator<=(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xC,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = operator<=(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xcC,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = operator<=(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xC,ycC)/dy", dy, ycC.adj());
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = operator<=(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xcC,ycC)/dy", dy, ycC.adj());
	}

	{
	xC = xd;
	Ai yAi(yd);
	fA = operator<=(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = operator<=(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xcC,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = operator<=(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xC,ycAi)/dy", dy, ycAi.aval);
	}

	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = operator<=(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xcC,ycAi)/dy", dy, ycAi.aval);
	}

	xC = xd;
	fA = operator<=(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,yd)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator<=(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,yd)/dx", dx, xcC.adj());
	}

	xC = xd;
	yL = (long)yd;
	fA = operator<=(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,yL)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = operator<=(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,yL)/dx", dx, xcC.adj());
	}

	xC = xd;
	yi = (int)yd;
	fA = operator<=(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator<=(xC,yi)/dx", dx, xC.adj());

	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = operator<=(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator<=(xcC,yi)/dx", dx, xcC.adj());
	}

	{
	Ai xAi(xd);
	yAI = yd;
	fA = operator<=(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = operator<=(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<=(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xcAi,yAI)/dy", dy, yAI.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = operator<=(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xAi,ycAI)/dy", dy, ycAI.adj());
	}

	{
	Ai xAi(xd);
	yA = yd;
	fA = operator<=(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator<=(xAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = operator<=(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<=(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator<=(xcAi,yA)/dy", dy, yA.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = operator<=(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xAi,ycA)/dy", dy, ycA.adj());
	}

	{
	Ai xAi(xd);
	yC = yd;
	fA = operator<=(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator<=(xAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = operator<=(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<=(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator<=(xcAi,yC)/dy", dy, yC.adj());
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = operator<=(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xAi,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = operator<=(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator<=(xAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = operator<=(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<=(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator<=(xcAi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = operator<=(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xAi,ycAi)/dy", dy, ycAi.aval);
	}

	{
	Ai xAi(xd);
	fA = operator<=(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,yd)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator<=(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<=(xcAi,yd)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = operator<=(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,yL)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = operator<=(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<=(xcAi,yL)/dx", dx, xcAi.aval);
	}

	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = operator<=(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator<=(xAi,yi)/dx", dx, xAi.aval);
	}

	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = operator<=(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator<=(xcAi,yi)/dx", dx, xcAi.aval);
	}

	yAI = yd;
	fA = operator<=(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xd,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = operator<=(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xd,ycAI)/dy", dy, ycAI.adj());
	}

	yA = yd;
	fA = operator<=(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xd,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	cA ycA(yd);
	fA = operator<=(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xd,ycA)/dy", dy, ycA.adj());
	}

	yC = yd;
	fA = operator<=(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xd,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	cC ycC(yd);
	fA = operator<=(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xd,ycC)/dy", dy, ycC.adj());
	}

	{
	Ai yAi(yd);
	fA = operator<=(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xd,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = operator<=(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xd,ycAi)/dy", dy, ycAi.aval);
	}

	xL = (long)xd;
	yAI = yd;
	fA = operator<=(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xL,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = operator<=(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xL,ycAI)/dy", dy, ycAI.adj());
	}

	xL = (long)xd;
	yA = yd;
	fA = operator<=(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xL,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = operator<=(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xL,ycA)/dy", dy, ycA.adj());
	}

	xL = (long)xd;
	yC = yd;
	fA = operator<=(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xL,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = operator<=(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xL,ycC)/dy", dy, ycC.adj());
	}

	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = operator<=(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xL,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = operator<=(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xL,ycAi)/dy", dy, ycAi.aval);
	}

	xi = (int)xd;
	yAI = yd;
	fA = operator<=(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator<=(xi,yAI)/dy", dy, yAI.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = operator<=(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator<=(xi,ycAI)/dy", dy, ycAI.adj());
	}

	xi = (int)xd;
	yA = yd;
	fA = operator<=(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator<=(xi,yA)/dy", dy, yA.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = operator<=(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator<=(xi,ycA)/dy", dy, ycA.adj());
	}

	xi = (int)xd;
	yC = yd;
	fA = operator<=(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator<=(xi,yC)/dy", dy, yC.adj());

	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = operator<=(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator<=(xi,ycC)/dy", dy, ycC.adj());
	}

	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = operator<=(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator<=(xi,yAi)/dy", dy, yAi.aval);
	}

	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = operator<=(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator<=(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator<=(xi,ycAi)/dy", dy, ycAi.aval);
	}


	if (!rc) // chatter for cppunit test, which cannot tolerate silence
		printf("OK\n");
	return rc;
	}
