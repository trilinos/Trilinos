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


	/**** Test of min ****/

	xd = 4.; yd = 3.; f = 3.; dx = 0.; dy = 1.;
	xAI = xd;
	yAI = yd;
	fA = min(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d min(xAI,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = min(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d min(xcAI,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = min(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d min(xAI,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = min(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d min(xcAI,ycAI)/dy", dy, ycAI.adj());
	}
	xAI = xd;
	yA = yd;
	fA = min(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d min(xAI,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = min(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d min(xcAI,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = min(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d min(xAI,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = min(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d min(xcAI,ycA)/dy", dy, ycA.adj());
	}
	xAI = xd;
	yC = yd;
	fA = min(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d min(xAI,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = min(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d min(xcAI,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = min(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d min(xAI,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = min(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d min(xcAI,ycC)/dy", dy, ycC.adj());
	}
	{
	xAI = xd;
	Ai yAi(yd);
	fA = min(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d min(xAI,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = min(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d min(xcAI,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = min(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d min(xAI,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = min(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d min(xcAI,ycAi)/dy", dy, ycAi.aval);
	}
	xAI = xd;
	fA = min(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,yd)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = min(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,yd)/dx", dx, xcAI.adj());
	}
	xAI = xd;
	yL = (long)yd;
	fA = min(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,yL)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = min(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,yL)/dx", dx, xcAI.adj());
	}
	xAI = xd;
	yi = (int)yd;
	fA = min(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d min(xAI,yi)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = min(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d min(xcAI,yi)/dx", dx, xcAI.adj());
	}
	xA = xd;
	yAI = yd;
	fA = min(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d min(xA,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = min(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d min(xcA,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = min(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d min(xA,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = min(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d min(xcA,ycAI)/dy", dy, ycAI.adj());
	}
	xA = xd;
	yA = yd;
	fA = min(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d min(xA,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = min(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d min(xcA,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = min(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d min(xA,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = min(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d min(xcA,ycA)/dy", dy, ycA.adj());
	}
	xA = xd;
	yC = yd;
	fA = min(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d min(xA,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = min(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d min(xcA,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = min(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d min(xA,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = min(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d min(xcA,ycC)/dy", dy, ycC.adj());
	}
	{
	xA = xd;
	Ai yAi(yd);
	fA = min(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d min(xA,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = min(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d min(xcA,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = min(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d min(xA,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = min(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d min(xcA,ycAi)/dy", dy, ycAi.aval);
	}
	xA = xd;
	fA = min(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,yd)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	fA = min(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,yd)/dx", dx, xcA.adj());
	}
	xA = xd;
	yL = (long)yd;
	fA = min(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,yL)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = min(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,yL)/dx", dx, xcA.adj());
	}
	xA = xd;
	yi = (int)yd;
	fA = min(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d min(xA,yi)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = min(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d min(xcA,yi)/dx", dx, xcA.adj());
	}
	xC = xd;
	yAI = yd;
	fA = min(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d min(xC,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = min(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d min(xcC,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = min(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d min(xC,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = min(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d min(xcC,ycAI)/dy", dy, ycAI.adj());
	}
	xC = xd;
	yA = yd;
	fA = min(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d min(xC,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = min(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d min(xcC,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = min(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d min(xC,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = min(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d min(xcC,ycA)/dy", dy, ycA.adj());
	}
	xC = xd;
	yC = yd;
	fA = min(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d min(xC,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = min(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d min(xcC,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = min(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d min(xC,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = min(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d min(xcC,ycC)/dy", dy, ycC.adj());
	}
	{
	xC = xd;
	Ai yAi(yd);
	fA = min(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d min(xC,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = min(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d min(xcC,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = min(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d min(xC,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = min(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d min(xcC,ycAi)/dy", dy, ycAi.aval);
	}
	xC = xd;
	fA = min(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,yd)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	fA = min(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,yd)/dx", dx, xcC.adj());
	}
	xC = xd;
	yL = (long)yd;
	fA = min(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,yL)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = min(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,yL)/dx", dx, xcC.adj());
	}
	xC = xd;
	yi = (int)yd;
	fA = min(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d min(xC,yi)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = min(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d min(xcC,yi)/dx", dx, xcC.adj());
	}
	{
	Ai xAi(xd);
	yAI = yd;
	fA = min(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d min(xAi,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = min(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d min(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d min(xcAi,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = min(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d min(xAi,ycAI)/dy", dy, ycAI.adj());
	}
	{
	Ai xAi(xd);
	yA = yd;
	fA = min(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d min(xAi,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = min(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d min(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d min(xcAi,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = min(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d min(xAi,ycA)/dy", dy, ycA.adj());
	}
	{
	Ai xAi(xd);
	yC = yd;
	fA = min(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d min(xAi,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = min(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d min(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d min(xcAi,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = min(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d min(xAi,ycC)/dy", dy, ycC.adj());
	}
	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = min(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d min(xAi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = min(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d min(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d min(xcAi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = min(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d min(xAi,ycAi)/dy", dy, ycAi.aval);
	}
	{
	Ai xAi(xd);
	fA = min(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,yd)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = min(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d min(xcAi,yd)/dx", dx, xcAi.aval);
	}
	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = min(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,yL)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = min(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d min(xcAi,yL)/dx", dx, xcAi.aval);
	}
	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = min(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d min(xAi,yi)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = min(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d min(xcAi,yi)/dx", dx, xcAi.aval);
	}
	yAI = yd;
	fA = min(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d min(xd,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = min(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d min(xd,ycAI)/dy", dy, ycAI.adj());
	}
	yA = yd;
	fA = min(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d min(xd,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cA ycA(yd);
	fA = min(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d min(xd,ycA)/dy", dy, ycA.adj());
	}
	yC = yd;
	fA = min(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d min(xd,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cC ycC(yd);
	fA = min(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d min(xd,ycC)/dy", dy, ycC.adj());
	}
	{
	Ai yAi(yd);
	fA = min(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d min(xd,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = min(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d min(xd,ycAi)/dy", dy, ycAi.aval);
	}
	xL = (long)xd;
	yAI = yd;
	fA = min(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d min(xL,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = min(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d min(xL,ycAI)/dy", dy, ycAI.adj());
	}
	xL = (long)xd;
	yA = yd;
	fA = min(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d min(xL,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = min(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d min(xL,ycA)/dy", dy, ycA.adj());
	}
	xL = (long)xd;
	yC = yd;
	fA = min(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d min(xL,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = min(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d min(xL,ycC)/dy", dy, ycC.adj());
	}
	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = min(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d min(xL,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = min(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d min(xL,ycAi)/dy", dy, ycAi.aval);
	}
	xi = (int)xd;
	yAI = yd;
	fA = min(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d min(xi,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = min(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d min(xi,ycAI)/dy", dy, ycAI.adj());
	}
	xi = (int)xd;
	yA = yd;
	fA = min(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d min(xi,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = min(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d min(xi,ycA)/dy", dy, ycA.adj());
	}
	xi = (int)xd;
	yC = yd;
	fA = min(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d min(xi,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = min(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d min(xi,ycC)/dy", dy, ycC.adj());
	}
	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = min(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d min(xi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = min(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = min(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d min(xi,ycAi)/dy", dy, ycAi.aval);
	}


	if (!rc) // chatter for cppunit test, which cannot tolerate silence

		printf("OK\n");

	return rc;

	}
