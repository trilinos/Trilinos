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


	if (!rc) // chatter for cppunit test, which cannot tolerate silence

		printf("OK\n");

	return rc;

	}
