// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/* Try to test all combinations of types and operations */



#define ADT_RAD Sacado::Rad::



#include "Sacado_trad.hpp"

#include <cstdio>

using std::printf;



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



#ifndef RAD_EQ_ALIAS

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


	/**** Test of operator!= ****/

	xd = 2.; yd = 3.; f = 1.; dx = 0.; dy = 0.;
	xAI = xd;
	yAI = yd;
	fA = operator!=(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xAI,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = operator!=(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xcAI,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = operator!=(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xAI,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = operator!=(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xcAI,ycAI)/dy", dy, ycAI.adj());
	}
	xAI = xd;
	yA = yd;
	fA = operator!=(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xAI,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = operator!=(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xcAI,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = operator!=(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xAI,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = operator!=(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xcAI,ycA)/dy", dy, ycA.adj());
	}
	xAI = xd;
	yC = yd;
	fA = operator!=(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xAI,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = operator!=(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xcAI,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = operator!=(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xAI,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = operator!=(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xcAI,ycC)/dy", dy, ycC.adj());
	}
	{
	xAI = xd;
	Ai yAi(yd);
	fA = operator!=(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xAI,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = operator!=(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xcAI,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = operator!=(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xAI,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = operator!=(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xcAI,ycAi)/dy", dy, ycAi.aval);
	}
	xAI = xd;
	fA = operator!=(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,yd)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = operator!=(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,yd)/dx", dx, xcAI.adj());
	}
	xAI = xd;
	yL = (long)yd;
	fA = operator!=(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,yL)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = operator!=(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,yL)/dx", dx, xcAI.adj());
	}
	xAI = xd;
	yi = (int)yd;
	fA = operator!=(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d operator!=(xAI,yi)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = operator!=(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d operator!=(xcAI,yi)/dx", dx, xcAI.adj());
	}
	xA = xd;
	yAI = yd;
	fA = operator!=(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xA,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = operator!=(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xcA,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = operator!=(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xA,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = operator!=(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xcA,ycAI)/dy", dy, ycAI.adj());
	}
	xA = xd;
	yA = yd;
	fA = operator!=(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xA,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = operator!=(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xcA,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = operator!=(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xA,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = operator!=(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xcA,ycA)/dy", dy, ycA.adj());
	}
	xA = xd;
	yC = yd;
	fA = operator!=(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xA,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = operator!=(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xcA,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = operator!=(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xA,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = operator!=(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xcA,ycC)/dy", dy, ycC.adj());
	}
	{
	xA = xd;
	Ai yAi(yd);
	fA = operator!=(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xA,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = operator!=(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xcA,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = operator!=(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xA,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = operator!=(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xcA,ycAi)/dy", dy, ycAi.aval);
	}
	xA = xd;
	fA = operator!=(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,yd)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	fA = operator!=(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,yd)/dx", dx, xcA.adj());
	}
	xA = xd;
	yL = (long)yd;
	fA = operator!=(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,yL)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = operator!=(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,yL)/dx", dx, xcA.adj());
	}
	xA = xd;
	yi = (int)yd;
	fA = operator!=(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d operator!=(xA,yi)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = operator!=(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d operator!=(xcA,yi)/dx", dx, xcA.adj());
	}
	xC = xd;
	yAI = yd;
	fA = operator!=(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xC,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = operator!=(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xcC,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = operator!=(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xC,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = operator!=(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xcC,ycAI)/dy", dy, ycAI.adj());
	}
	xC = xd;
	yA = yd;
	fA = operator!=(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xC,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = operator!=(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xcC,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = operator!=(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xC,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = operator!=(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xcC,ycA)/dy", dy, ycA.adj());
	}
	xC = xd;
	yC = yd;
	fA = operator!=(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xC,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = operator!=(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xcC,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = operator!=(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xC,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = operator!=(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xcC,ycC)/dy", dy, ycC.adj());
	}
	{
	xC = xd;
	Ai yAi(yd);
	fA = operator!=(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xC,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = operator!=(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xcC,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = operator!=(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xC,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = operator!=(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xcC,ycAi)/dy", dy, ycAi.aval);
	}
	xC = xd;
	fA = operator!=(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,yd)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	fA = operator!=(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,yd)/dx", dx, xcC.adj());
	}
	xC = xd;
	yL = (long)yd;
	fA = operator!=(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,yL)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = operator!=(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,yL)/dx", dx, xcC.adj());
	}
	xC = xd;
	yi = (int)yd;
	fA = operator!=(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d operator!=(xC,yi)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = operator!=(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d operator!=(xcC,yi)/dx", dx, xcC.adj());
	}
	{
	Ai xAi(xd);
	yAI = yd;
	fA = operator!=(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xAi,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = operator!=(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator!=(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xcAi,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = operator!=(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xAi,ycAI)/dy", dy, ycAI.adj());
	}
	{
	Ai xAi(xd);
	yA = yd;
	fA = operator!=(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator!=(xAi,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = operator!=(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator!=(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d operator!=(xcAi,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = operator!=(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xAi,ycA)/dy", dy, ycA.adj());
	}
	{
	Ai xAi(xd);
	yC = yd;
	fA = operator!=(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator!=(xAi,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = operator!=(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator!=(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d operator!=(xcAi,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = operator!=(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xAi,ycC)/dy", dy, ycC.adj());
	}
	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = operator!=(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator!=(xAi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = operator!=(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator!=(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d operator!=(xcAi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = operator!=(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xAi,ycAi)/dy", dy, ycAi.aval);
	}
	{
	Ai xAi(xd);
	fA = operator!=(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,yd)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = operator!=(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator!=(xcAi,yd)/dx", dx, xcAi.aval);
	}
	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = operator!=(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,yL)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = operator!=(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator!=(xcAi,yL)/dx", dx, xcAi.aval);
	}
	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = operator!=(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d operator!=(xAi,yi)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = operator!=(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d operator!=(xcAi,yi)/dx", dx, xcAi.aval);
	}
	yAI = yd;
	fA = operator!=(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xd,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = operator!=(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xd,ycAI)/dy", dy, ycAI.adj());
	}
	yA = yd;
	fA = operator!=(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xd,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cA ycA(yd);
	fA = operator!=(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xd,ycA)/dy", dy, ycA.adj());
	}
	yC = yd;
	fA = operator!=(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xd,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cC ycC(yd);
	fA = operator!=(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xd,ycC)/dy", dy, ycC.adj());
	}
	{
	Ai yAi(yd);
	fA = operator!=(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xd,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = operator!=(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xd,ycAi)/dy", dy, ycAi.aval);
	}
	xL = (long)xd;
	yAI = yd;
	fA = operator!=(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xL,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = operator!=(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xL,ycAI)/dy", dy, ycAI.adj());
	}
	xL = (long)xd;
	yA = yd;
	fA = operator!=(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xL,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = operator!=(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xL,ycA)/dy", dy, ycA.adj());
	}
	xL = (long)xd;
	yC = yd;
	fA = operator!=(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xL,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = operator!=(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xL,ycC)/dy", dy, ycC.adj());
	}
	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = operator!=(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xL,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = operator!=(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xL,ycAi)/dy", dy, ycAi.aval);
	}
	xi = (int)xd;
	yAI = yd;
	fA = operator!=(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d operator!=(xi,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = operator!=(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d operator!=(xi,ycAI)/dy", dy, ycAI.adj());
	}
	xi = (int)xd;
	yA = yd;
	fA = operator!=(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d operator!=(xi,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = operator!=(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d operator!=(xi,ycA)/dy", dy, ycA.adj());
	}
	xi = (int)xd;
	yC = yd;
	fA = operator!=(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d operator!=(xi,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = operator!=(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d operator!=(xi,ycC)/dy", dy, ycC.adj());
	}
	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = operator!=(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d operator!=(xi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = operator!=(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = operator!=(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d operator!=(xi,ycAi)/dy", dy, ycAi.aval);
	}


	if (!rc) // chatter for cppunit test, which cannot tolerate silence

		printf("OK\n");

	return rc;

	}
