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


	if (!rc) // chatter for cppunit test, which cannot tolerate silence

		printf("OK\n");

	return rc;

	}
