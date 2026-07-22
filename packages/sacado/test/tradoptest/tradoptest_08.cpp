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


	/**** Test of max ****/

	xd = 4.; yd = 5.; f = 5.; dx = 0.; dy = 1.;
	xAI = xd;
	yAI = yd;
	fA = max(xAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,yAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,yAI)/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d max(xAI,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = max(xcAI,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,yAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,yAI)/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d max(xcAI,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = max(xAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,ycAI)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,ycAI)/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d max(xAI,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = max(xcAI,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,ycAI)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,ycAI)/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d max(xcAI,ycAI)/dy", dy, ycAI.adj());
	}
	xAI = xd;
	yA = yd;
	fA = max(xAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,yA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,yA)/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d max(xAI,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = max(xcAI,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,yA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,yA)/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d max(xcAI,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = max(xAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,ycA)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,ycA)/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d max(xAI,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = max(xcAI,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,ycA)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,ycA)/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d max(xcAI,ycA)/dy", dy, ycA.adj());
	}
	xAI = xd;
	yC = yd;
	fA = max(xAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,yC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,yC)/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d max(xAI,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = max(xcAI,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,yC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,yC)/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d max(xcAI,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = max(xAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,ycC)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,ycC)/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d max(xAI,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = max(xcAI,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,ycC)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,ycC)/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d max(xcAI,ycC)/dy", dy, ycC.adj());
	}
	{
	xAI = xd;
	Ai yAi(yd);
	fA = max(xAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,yAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,yAi)/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d max(xAI,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = max(xcAI,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,yAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,yAi)/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d max(xcAI,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = max(xAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,ycAi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,ycAi)/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d max(xAI,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = max(xcAI,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,ycAi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,ycAi)/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d max(xcAI,ycAi)/dy", dy, ycAi.aval);
	}
	xAI = xd;
	fA = max(xAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,yd)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,yd)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = max(xcAI,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,yd)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,yd)/dx", dx, xcAI.adj());
	}
	xAI = xd;
	yL = (long)yd;
	fA = max(xAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,yL)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,yL)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = max(xcAI,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,yL)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,yL)/dx", dx, xcAI.adj());
	}
	xAI = xd;
	yi = (int)yd;
	fA = max(xAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAI,yi)", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d max(xAI,yi)/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = max(xcAI,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAI,yi)", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d max(xcAI,yi)/dx", dx, xcAI.adj());
	}
	xA = xd;
	yAI = yd;
	fA = max(xA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,yAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,yAI)/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d max(xA,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = max(xcA,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,yAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,yAI)/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d max(xcA,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = max(xA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,ycAI)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,ycAI)/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d max(xA,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = max(xcA,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,ycAI)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,ycAI)/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d max(xcA,ycAI)/dy", dy, ycAI.adj());
	}
	xA = xd;
	yA = yd;
	fA = max(xA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,yA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,yA)/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d max(xA,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = max(xcA,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,yA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,yA)/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d max(xcA,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = max(xA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,ycA)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,ycA)/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d max(xA,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = max(xcA,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,ycA)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,ycA)/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d max(xcA,ycA)/dy", dy, ycA.adj());
	}
	xA = xd;
	yC = yd;
	fA = max(xA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,yC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,yC)/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d max(xA,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = max(xcA,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,yC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,yC)/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d max(xcA,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = max(xA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,ycC)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,ycC)/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d max(xA,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = max(xcA,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,ycC)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,ycC)/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d max(xcA,ycC)/dy", dy, ycC.adj());
	}
	{
	xA = xd;
	Ai yAi(yd);
	fA = max(xA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,yAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,yAi)/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d max(xA,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = max(xcA,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,yAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,yAi)/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d max(xcA,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = max(xA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,ycAi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,ycAi)/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d max(xA,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = max(xcA,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,ycAi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,ycAi)/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d max(xcA,ycAi)/dy", dy, ycAi.aval);
	}
	xA = xd;
	fA = max(xA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,yd)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,yd)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	fA = max(xcA,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,yd)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,yd)/dx", dx, xcA.adj());
	}
	xA = xd;
	yL = (long)yd;
	fA = max(xA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,yL)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,yL)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = max(xcA,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,yL)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,yL)/dx", dx, xcA.adj());
	}
	xA = xd;
	yi = (int)yd;
	fA = max(xA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xA,yi)", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d max(xA,yi)/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = max(xcA,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcA,yi)", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d max(xcA,yi)/dx", dx, xcA.adj());
	}
	xC = xd;
	yAI = yd;
	fA = max(xC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,yAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,yAI)/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d max(xC,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = max(xcC,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,yAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,yAI)/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d max(xcC,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = max(xC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,ycAI)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,ycAI)/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d max(xC,ycAI)/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = max(xcC,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,ycAI)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,ycAI)/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d max(xcC,ycAI)/dy", dy, ycAI.adj());
	}
	xC = xd;
	yA = yd;
	fA = max(xC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,yA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,yA)/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d max(xC,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = max(xcC,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,yA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,yA)/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d max(xcC,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = max(xC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,ycA)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,ycA)/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d max(xC,ycA)/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = max(xcC,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,ycA)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,ycA)/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d max(xcC,ycA)/dy", dy, ycA.adj());
	}
	xC = xd;
	yC = yd;
	fA = max(xC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,yC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,yC)/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d max(xC,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = max(xcC,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,yC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,yC)/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d max(xcC,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = max(xC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,ycC)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,ycC)/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d max(xC,ycC)/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = max(xcC,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,ycC)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,ycC)/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d max(xcC,ycC)/dy", dy, ycC.adj());
	}
	{
	xC = xd;
	Ai yAi(yd);
	fA = max(xC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,yAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,yAi)/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d max(xC,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = max(xcC,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,yAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,yAi)/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d max(xcC,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = max(xC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,ycAi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,ycAi)/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d max(xC,ycAi)/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = max(xcC,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,ycAi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,ycAi)/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d max(xcC,ycAi)/dy", dy, ycAi.aval);
	}
	xC = xd;
	fA = max(xC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,yd)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,yd)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	fA = max(xcC,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,yd)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,yd)/dx", dx, xcC.adj());
	}
	xC = xd;
	yL = (long)yd;
	fA = max(xC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,yL)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,yL)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = max(xcC,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,yL)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,yL)/dx", dx, xcC.adj());
	}
	xC = xd;
	yi = (int)yd;
	fA = max(xC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xC,yi)", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d max(xC,yi)/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = max(xcC,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcC,yi)", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d max(xcC,yi)/dx", dx, xcC.adj());
	}
	{
	Ai xAi(xd);
	yAI = yd;
	fA = max(xAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,yAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,yAI)/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d max(xAi,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = max(xcAi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAi,yAI)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d max(xcAi,yAI)/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d max(xcAi,yAI)/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = max(xAi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,ycAI)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,ycAI)/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d max(xAi,ycAI)/dy", dy, ycAI.adj());
	}
	{
	Ai xAi(xd);
	yA = yd;
	fA = max(xAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,yA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,yA)/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d max(xAi,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = max(xcAi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAi,yA)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d max(xcAi,yA)/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d max(xcAi,yA)/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = max(xAi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,ycA)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,ycA)/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d max(xAi,ycA)/dy", dy, ycA.adj());
	}
	{
	Ai xAi(xd);
	yC = yd;
	fA = max(xAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,yC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,yC)/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d max(xAi,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = max(xcAi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAi,yC)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d max(xcAi,yC)/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d max(xcAi,yC)/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = max(xAi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,ycC)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,ycC)/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d max(xAi,ycC)/dy", dy, ycC.adj());
	}
	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = max(xAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,yAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,yAi)/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d max(xAi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = max(xcAi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAi,yAi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d max(xcAi,yAi)/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d max(xcAi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = max(xAi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,ycAi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,ycAi)/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d max(xAi,ycAi)/dy", dy, ycAi.aval);
	}
	{
	Ai xAi(xd);
	fA = max(xAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,yd)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,yd)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = max(xcAi,yd);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAi,yd)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d max(xcAi,yd)/dx", dx, xcAi.aval);
	}
	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = max(xAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,yL)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,yL)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = max(xcAi,yL);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAi,yL)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d max(xcAi,yL)/dx", dx, xcAi.aval);
	}
	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = max(xAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xAi,yi)", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d max(xAi,yi)/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = max(xcAi,yi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xcAi,yi)", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d max(xcAi,yi)/dx", dx, xcAi.aval);
	}
	yAI = yd;
	fA = max(xd,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xd,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d max(xd,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = max(xd,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xd,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d max(xd,ycAI)/dy", dy, ycAI.adj());
	}
	yA = yd;
	fA = max(xd,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xd,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d max(xd,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	cA ycA(yd);
	fA = max(xd,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xd,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d max(xd,ycA)/dy", dy, ycA.adj());
	}
	yC = yd;
	fA = max(xd,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xd,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d max(xd,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	cC ycC(yd);
	fA = max(xd,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xd,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d max(xd,ycC)/dy", dy, ycC.adj());
	}
	{
	Ai yAi(yd);
	fA = max(xd,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xd,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d max(xd,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = max(xd,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xd,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d max(xd,ycAi)/dy", dy, ycAi.aval);
	}
	xL = (long)xd;
	yAI = yd;
	fA = max(xL,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xL,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d max(xL,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = max(xL,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xL,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d max(xL,ycAI)/dy", dy, ycAI.adj());
	}
	xL = (long)xd;
	yA = yd;
	fA = max(xL,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xL,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d max(xL,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = max(xL,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xL,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d max(xL,ycA)/dy", dy, ycA.adj());
	}
	xL = (long)xd;
	yC = yd;
	fA = max(xL,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xL,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d max(xL,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = max(xL,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xL,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d max(xL,ycC)/dy", dy, ycC.adj());
	}
	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = max(xL,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xL,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d max(xL,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = max(xL,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xL,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d max(xL,ycAi)/dy", dy, ycAi.aval);
	}
	xi = (int)xd;
	yAI = yd;
	fA = max(xi,yAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xi,yAI)", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d max(xi,yAI)/dy", dy, yAI.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = max(xi,ycAI);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xi,ycAI)", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d max(xi,ycAI)/dy", dy, ycAI.adj());
	}
	xi = (int)xd;
	yA = yd;
	fA = max(xi,yA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xi,yA)", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d max(xi,yA)/dy", dy, yA.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = max(xi,ycA);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xi,ycA)", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d max(xi,ycA)/dy", dy, ycA.adj());
	}
	xi = (int)xd;
	yC = yd;
	fA = max(xi,yC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xi,yC)", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d max(xi,yC)/dy", dy, yC.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = max(xi,ycC);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xi,ycC)", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d max(xi,ycC)/dy", dy, ycC.adj());
	}
	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = max(xi,yAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xi,yAi)", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d max(xi,yAi)/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = max(xi,ycAi);
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = max(xi,ycAi)", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d max(xi,ycAi)/dy", dy, ycAi.aval);
	}


	if (!rc) // chatter for cppunit test, which cannot tolerate silence

		printf("OK\n");

	return rc;

	}
