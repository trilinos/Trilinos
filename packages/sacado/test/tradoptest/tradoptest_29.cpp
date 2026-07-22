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


	/**** Test of <= ****/

	xd = 2.; yd = 3.; f = 1.; dx = 0.; dy = 0.;
	xAI = xd;
	yAI = yd;
	fA = xAI <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= yAI", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= yAI/dx", dx, xAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d xAI <= yAI/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yAI = yd;
	fA = xcAI <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= yAI", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= yAI/dx", dx, xcAI.adj());
	else if (differ(yAI.adj(), dy)) botch("d xcAI <= yAI/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cAI ycAI(yd);
	fA = xAI <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= ycAI", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= ycAI/dx", dx, xAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d xAI <= ycAI/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cAI ycAI(yd);
	fA = xcAI <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= ycAI", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= ycAI/dx", dx, xcAI.adj());
	else if (differ(ycAI.adj(), dy)) botch("d xcAI <= ycAI/dy", dy, ycAI.adj());
	}
	xAI = xd;
	yA = yd;
	fA = xAI <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= yA", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= yA/dx", dx, xAI.adj());
	else if (differ(yA.adj(), dy)) botch("d xAI <= yA/dy", dy, yA.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yA = yd;
	fA = xcAI <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= yA", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= yA/dx", dx, xcAI.adj());
	else if (differ(yA.adj(), dy)) botch("d xcAI <= yA/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cA ycA(yd);
	fA = xAI <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= ycA", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= ycA/dx", dx, xAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d xAI <= ycA/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cA ycA(yd);
	fA = xcAI <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= ycA", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= ycA/dx", dx, xcAI.adj());
	else if (differ(ycA.adj(), dy)) botch("d xcAI <= ycA/dy", dy, ycA.adj());
	}
	xAI = xd;
	yC = yd;
	fA = xAI <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= yC", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= yC/dx", dx, xAI.adj());
	else if (differ(yC.adj(), dy)) botch("d xAI <= yC/dy", dy, yC.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yC = yd;
	fA = xcAI <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= yC", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= yC/dx", dx, xcAI.adj());
	else if (differ(yC.adj(), dy)) botch("d xcAI <= yC/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xAI = xd;
	cC ycC(yd);
	fA = xAI <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= ycC", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= ycC/dx", dx, xAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d xAI <= ycC/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cC ycC(yd);
	fA = xcAI <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= ycC", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= ycC/dx", dx, xcAI.adj());
	else if (differ(ycC.adj(), dy)) botch("d xcAI <= ycC/dy", dy, ycC.adj());
	}
	{
	xAI = xd;
	Ai yAi(yd);
	fA = xAI <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= yAi", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= yAi/dx", dx, xAI.adj());
	else if (differ(yAi.aval, dy)) botch("d xAI <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	Ai yAi(yd);
	fA = xcAI <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= yAi", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= yAi/dx", dx, xcAI.adj());
	else if (differ(yAi.aval, dy)) botch("d xcAI <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xAI = xd;
	cAi ycAi(yd);
	fA = xAI <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= ycAi", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= ycAi/dx", dx, xAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d xAI <= ycAi/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cAI xcAI(xd);
	cAi ycAi(yd);
	fA = xcAI <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= ycAi", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= ycAi/dx", dx, xcAI.adj());
	else if (differ(ycAi.aval, dy)) botch("d xcAI <= ycAi/dy", dy, ycAi.aval);
	}
	xAI = xd;
	fA = xAI <= yd;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= yd", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= yd/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	fA = xcAI <= yd;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= yd", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= yd/dx", dx, xcAI.adj());
	}
	xAI = xd;
	yL = (long)yd;
	fA = xAI <= yL;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= yL", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= yL/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yL = (long)yd;
	fA = xcAI <= yL;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= yL", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= yL/dx", dx, xcAI.adj());
	}
	xAI = xd;
	yi = (int)yd;
	fA = xAI <= yi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAI <= yi", f, fA.val());
	else if (differ(xAI.adj(), dx)) botch("d xAI <= yi/dx", dx, xAI.adj());
	{
	A::aval_reset();
	cAI xcAI(xd);
	yi = (int)yd;
	fA = xcAI <= yi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAI <= yi", f, fA.val());
	else if (differ(xcAI.adj(), dx)) botch("d xcAI <= yi/dx", dx, xcAI.adj());
	}
	xA = xd;
	yAI = yd;
	fA = xA <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= yAI", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= yAI/dx", dx, xA.adj());
	else if (differ(yAI.adj(), dy)) botch("d xA <= yAI/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yAI = yd;
	fA = xcA <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= yAI", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= yAI/dx", dx, xcA.adj());
	else if (differ(yAI.adj(), dy)) botch("d xcA <= yAI/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cAI ycAI(yd);
	fA = xA <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= ycAI", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= ycAI/dx", dx, xA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d xA <= ycAI/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cAI ycAI(yd);
	fA = xcA <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= ycAI", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= ycAI/dx", dx, xcA.adj());
	else if (differ(ycAI.adj(), dy)) botch("d xcA <= ycAI/dy", dy, ycAI.adj());
	}
	xA = xd;
	yA = yd;
	fA = xA <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= yA", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= yA/dx", dx, xA.adj());
	else if (differ(yA.adj(), dy)) botch("d xA <= yA/dy", dy, yA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yA = yd;
	fA = xcA <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= yA", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= yA/dx", dx, xcA.adj());
	else if (differ(yA.adj(), dy)) botch("d xcA <= yA/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cA ycA(yd);
	fA = xA <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= ycA", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= ycA/dx", dx, xA.adj());
	else if (differ(ycA.adj(), dy)) botch("d xA <= ycA/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cA ycA(yd);
	fA = xcA <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= ycA", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= ycA/dx", dx, xcA.adj());
	else if (differ(ycA.adj(), dy)) botch("d xcA <= ycA/dy", dy, ycA.adj());
	}
	xA = xd;
	yC = yd;
	fA = xA <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= yC", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= yC/dx", dx, xA.adj());
	else if (differ(yC.adj(), dy)) botch("d xA <= yC/dy", dy, yC.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yC = yd;
	fA = xcA <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= yC", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= yC/dx", dx, xcA.adj());
	else if (differ(yC.adj(), dy)) botch("d xcA <= yC/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xA = xd;
	cC ycC(yd);
	fA = xA <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= ycC", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= ycC/dx", dx, xA.adj());
	else if (differ(ycC.adj(), dy)) botch("d xA <= ycC/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cC ycC(yd);
	fA = xcA <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= ycC", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= ycC/dx", dx, xcA.adj());
	else if (differ(ycC.adj(), dy)) botch("d xcA <= ycC/dy", dy, ycC.adj());
	}
	{
	xA = xd;
	Ai yAi(yd);
	fA = xA <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= yAi", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= yAi/dx", dx, xA.adj());
	else if (differ(yAi.aval, dy)) botch("d xA <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cA xcA(xd);
	Ai yAi(yd);
	fA = xcA <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= yAi", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= yAi/dx", dx, xcA.adj());
	else if (differ(yAi.aval, dy)) botch("d xcA <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xA = xd;
	cAi ycAi(yd);
	fA = xA <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= ycAi", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= ycAi/dx", dx, xA.adj());
	else if (differ(ycAi.aval, dy)) botch("d xA <= ycAi/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cA xcA(xd);
	cAi ycAi(yd);
	fA = xcA <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= ycAi", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= ycAi/dx", dx, xcA.adj());
	else if (differ(ycAi.aval, dy)) botch("d xcA <= ycAi/dy", dy, ycAi.aval);
	}
	xA = xd;
	fA = xA <= yd;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= yd", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= yd/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	fA = xcA <= yd;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= yd", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= yd/dx", dx, xcA.adj());
	}
	xA = xd;
	yL = (long)yd;
	fA = xA <= yL;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= yL", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= yL/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yL = (long)yd;
	fA = xcA <= yL;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= yL", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= yL/dx", dx, xcA.adj());
	}
	xA = xd;
	yi = (int)yd;
	fA = xA <= yi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xA <= yi", f, fA.val());
	else if (differ(xA.adj(), dx)) botch("d xA <= yi/dx", dx, xA.adj());
	{
	A::aval_reset();
	cA xcA(xd);
	yi = (int)yd;
	fA = xcA <= yi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcA <= yi", f, fA.val());
	else if (differ(xcA.adj(), dx)) botch("d xcA <= yi/dx", dx, xcA.adj());
	}
	xC = xd;
	yAI = yd;
	fA = xC <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= yAI", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= yAI/dx", dx, xC.adj());
	else if (differ(yAI.adj(), dy)) botch("d xC <= yAI/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yAI = yd;
	fA = xcC <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= yAI", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= yAI/dx", dx, xcC.adj());
	else if (differ(yAI.adj(), dy)) botch("d xcC <= yAI/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cAI ycAI(yd);
	fA = xC <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= ycAI", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= ycAI/dx", dx, xC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d xC <= ycAI/dy", dy, ycAI.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cAI ycAI(yd);
	fA = xcC <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= ycAI", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= ycAI/dx", dx, xcC.adj());
	else if (differ(ycAI.adj(), dy)) botch("d xcC <= ycAI/dy", dy, ycAI.adj());
	}
	xC = xd;
	yA = yd;
	fA = xC <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= yA", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= yA/dx", dx, xC.adj());
	else if (differ(yA.adj(), dy)) botch("d xC <= yA/dy", dy, yA.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yA = yd;
	fA = xcC <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= yA", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= yA/dx", dx, xcC.adj());
	else if (differ(yA.adj(), dy)) botch("d xcC <= yA/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cA ycA(yd);
	fA = xC <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= ycA", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= ycA/dx", dx, xC.adj());
	else if (differ(ycA.adj(), dy)) botch("d xC <= ycA/dy", dy, ycA.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cA ycA(yd);
	fA = xcC <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= ycA", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= ycA/dx", dx, xcC.adj());
	else if (differ(ycA.adj(), dy)) botch("d xcC <= ycA/dy", dy, ycA.adj());
	}
	xC = xd;
	yC = yd;
	fA = xC <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= yC", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= yC/dx", dx, xC.adj());
	else if (differ(yC.adj(), dy)) botch("d xC <= yC/dy", dy, yC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yC = yd;
	fA = xcC <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= yC", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= yC/dx", dx, xcC.adj());
	else if (differ(yC.adj(), dy)) botch("d xcC <= yC/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	xC = xd;
	cC ycC(yd);
	fA = xC <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= ycC", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= ycC/dx", dx, xC.adj());
	else if (differ(ycC.adj(), dy)) botch("d xC <= ycC/dy", dy, ycC.adj());
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cC ycC(yd);
	fA = xcC <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= ycC", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= ycC/dx", dx, xcC.adj());
	else if (differ(ycC.adj(), dy)) botch("d xcC <= ycC/dy", dy, ycC.adj());
	}
	{
	xC = xd;
	Ai yAi(yd);
	fA = xC <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= yAi", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= yAi/dx", dx, xC.adj());
	else if (differ(yAi.aval, dy)) botch("d xC <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cC xcC(xd);
	Ai yAi(yd);
	fA = xcC <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= yAi", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= yAi/dx", dx, xcC.adj());
	else if (differ(yAi.aval, dy)) botch("d xcC <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xC = xd;
	cAi ycAi(yd);
	fA = xC <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= ycAi", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= ycAi/dx", dx, xC.adj());
	else if (differ(ycAi.aval, dy)) botch("d xC <= ycAi/dy", dy, ycAi.aval);
	}
	{
	A::aval_reset();
	cC xcC(xd);
	cAi ycAi(yd);
	fA = xcC <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= ycAi", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= ycAi/dx", dx, xcC.adj());
	else if (differ(ycAi.aval, dy)) botch("d xcC <= ycAi/dy", dy, ycAi.aval);
	}
	xC = xd;
	fA = xC <= yd;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= yd", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= yd/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	fA = xcC <= yd;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= yd", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= yd/dx", dx, xcC.adj());
	}
	xC = xd;
	yL = (long)yd;
	fA = xC <= yL;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= yL", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= yL/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yL = (long)yd;
	fA = xcC <= yL;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= yL", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= yL/dx", dx, xcC.adj());
	}
	xC = xd;
	yi = (int)yd;
	fA = xC <= yi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xC <= yi", f, fA.val());
	else if (differ(xC.adj(), dx)) botch("d xC <= yi/dx", dx, xC.adj());
	{
	A::aval_reset();
	cC xcC(xd);
	yi = (int)yd;
	fA = xcC <= yi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcC <= yi", f, fA.val());
	else if (differ(xcC.adj(), dx)) botch("d xcC <= yi/dx", dx, xcC.adj());
	}
	{
	Ai xAi(xd);
	yAI = yd;
	fA = xAi <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= yAI", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= yAI/dx", dx, xAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d xAi <= yAI/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yAI = yd;
	fA = xcAi <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAi <= yAI", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d xcAi <= yAI/dx", dx, xcAi.aval);
	else if (differ(yAI.adj(), dy)) botch("d xcAi <= yAI/dy", dy, yAI.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cAI ycAI(yd);
	fA = xAi <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= ycAI", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= ycAI/dx", dx, xAi.aval);
	else if (differ(ycAI.adj(), dy)) botch("d xAi <= ycAI/dy", dy, ycAI.adj());
	}
	{
	Ai xAi(xd);
	yA = yd;
	fA = xAi <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= yA", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= yA/dx", dx, xAi.aval);
	else if (differ(yA.adj(), dy)) botch("d xAi <= yA/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yA = yd;
	fA = xcAi <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAi <= yA", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d xcAi <= yA/dx", dx, xcAi.aval);
	else if (differ(yA.adj(), dy)) botch("d xcAi <= yA/dy", dy, yA.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cA ycA(yd);
	fA = xAi <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= ycA", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= ycA/dx", dx, xAi.aval);
	else if (differ(ycA.adj(), dy)) botch("d xAi <= ycA/dy", dy, ycA.adj());
	}
	{
	Ai xAi(xd);
	yC = yd;
	fA = xAi <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= yC", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= yC/dx", dx, xAi.aval);
	else if (differ(yC.adj(), dy)) botch("d xAi <= yC/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yC = yd;
	fA = xcAi <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAi <= yC", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d xcAi <= yC/dx", dx, xcAi.aval);
	else if (differ(yC.adj(), dy)) botch("d xcAi <= yC/dy", dy, yC.adj());
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cC ycC(yd);
	fA = xAi <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= ycC", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= ycC/dx", dx, xAi.aval);
	else if (differ(ycC.adj(), dy)) botch("d xAi <= ycC/dy", dy, ycC.adj());
	}
	{
	Ai xAi(xd);
	Ai yAi(yd);
	fA = xAi <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= yAi", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= yAi/dx", dx, xAi.aval);
	else if (differ(yAi.aval, dy)) botch("d xAi <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	Ai yAi(yd);
	fA = xcAi <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAi <= yAi", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d xcAi <= yAi/dx", dx, xcAi.aval);
	else if (differ(yAi.aval, dy)) botch("d xcAi <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	Ai xAi(xd);
	cAi ycAi(yd);
	fA = xAi <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= ycAi", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= ycAi/dx", dx, xAi.aval);
	else if (differ(ycAi.aval, dy)) botch("d xAi <= ycAi/dy", dy, ycAi.aval);
	}
	{
	Ai xAi(xd);
	fA = xAi <= yd;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= yd", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= yd/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	fA = xcAi <= yd;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAi <= yd", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d xcAi <= yd/dx", dx, xcAi.aval);
	}
	{
	Ai xAi(xd);
	yL = (long)yd;
	fA = xAi <= yL;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= yL", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= yL/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yL = (long)yd;
	fA = xcAi <= yL;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAi <= yL", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d xcAi <= yL/dx", dx, xcAi.aval);
	}
	{
	Ai xAi(xd);
	yi = (int)yd;
	fA = xAi <= yi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xAi <= yi", f, fA.val());
	else if (differ(xAi.aval, dx)) botch("d xAi <= yi/dx", dx, xAi.aval);
	}
	{
	A::aval_reset();
	cAi xcAi(xd);
	yi = (int)yd;
	fA = xcAi <= yi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xcAi <= yi", f, fA.val());
	else if (differ(xcAi.aval, dx)) botch("d xcAi <= yi/dx", dx, xcAi.aval);
	}
	yAI = yd;
	fA = xd <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xd <= yAI", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d xd <= yAI/dy", dy, yAI.adj());
	{
	A::aval_reset();
	cAI ycAI(yd);
	fA = xd <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xd <= ycAI", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d xd <= ycAI/dy", dy, ycAI.adj());
	}
	yA = yd;
	fA = xd <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xd <= yA", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d xd <= yA/dy", dy, yA.adj());
	{
	A::aval_reset();
	cA ycA(yd);
	fA = xd <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xd <= ycA", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d xd <= ycA/dy", dy, ycA.adj());
	}
	yC = yd;
	fA = xd <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xd <= yC", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d xd <= yC/dy", dy, yC.adj());
	{
	A::aval_reset();
	cC ycC(yd);
	fA = xd <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xd <= ycC", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d xd <= ycC/dy", dy, ycC.adj());
	}
	{
	Ai yAi(yd);
	fA = xd <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xd <= yAi", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d xd <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	cAi ycAi(yd);
	fA = xd <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xd <= ycAi", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d xd <= ycAi/dy", dy, ycAi.aval);
	}
	xL = (long)xd;
	yAI = yd;
	fA = xL <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xL <= yAI", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d xL <= yAI/dy", dy, yAI.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cAI ycAI(yd);
	fA = xL <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xL <= ycAI", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d xL <= ycAI/dy", dy, ycAI.adj());
	}
	xL = (long)xd;
	yA = yd;
	fA = xL <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xL <= yA", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d xL <= yA/dy", dy, yA.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cA ycA(yd);
	fA = xL <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xL <= ycA", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d xL <= ycA/dy", dy, ycA.adj());
	}
	xL = (long)xd;
	yC = yd;
	fA = xL <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xL <= yC", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d xL <= yC/dy", dy, yC.adj());
	{
	A::aval_reset();
	xL = (long)xd;
	cC ycC(yd);
	fA = xL <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xL <= ycC", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d xL <= ycC/dy", dy, ycC.adj());
	}
	{
	xL = (long)xd;
	Ai yAi(yd);
	fA = xL <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xL <= yAi", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d xL <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xL = (long)xd;
	cAi ycAi(yd);
	fA = xL <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xL <= ycAi", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d xL <= ycAi/dy", dy, ycAi.aval);
	}
	xi = (int)xd;
	yAI = yd;
	fA = xi <= yAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xi <= yAI", f, fA.val());
	else if (differ(yAI.adj(), dy)) botch("d xi <= yAI/dy", dy, yAI.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cAI ycAI(yd);
	fA = xi <= ycAI;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xi <= ycAI", f, fA.val());
	else if (differ(ycAI.adj(), dy)) botch("d xi <= ycAI/dy", dy, ycAI.adj());
	}
	xi = (int)xd;
	yA = yd;
	fA = xi <= yA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xi <= yA", f, fA.val());
	else if (differ(yA.adj(), dy)) botch("d xi <= yA/dy", dy, yA.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cA ycA(yd);
	fA = xi <= ycA;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xi <= ycA", f, fA.val());
	else if (differ(ycA.adj(), dy)) botch("d xi <= ycA/dy", dy, ycA.adj());
	}
	xi = (int)xd;
	yC = yd;
	fA = xi <= yC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xi <= yC", f, fA.val());
	else if (differ(yC.adj(), dy)) botch("d xi <= yC/dy", dy, yC.adj());
	{
	A::aval_reset();
	xi = (int)xd;
	cC ycC(yd);
	fA = xi <= ycC;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xi <= ycC", f, fA.val());
	else if (differ(ycC.adj(), dy)) botch("d xi <= ycC/dy", dy, ycC.adj());
	}
	{
	xi = (int)xd;
	Ai yAi(yd);
	fA = xi <= yAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xi <= yAi", f, fA.val());
	else if (differ(yAi.aval, dy)) botch("d xi <= yAi/dy", dy, yAi.aval);
	}
	{
	A::aval_reset();
	xi = (int)xd;
	cAi ycAi(yd);
	fA = xi <= ycAi;
	A::Gradcomp();
	if (differ(fA.val(), f)) botch("fA = xi <= ycAi", f, fA.val());
	else if (differ(ycAi.aval, dy)) botch("d xi <= ycAi/dy", dy, ycAi.aval);
	}


	if (!rc) // chatter for cppunit test, which cannot tolerate silence

		printf("OK\n");

	return rc;

	}
