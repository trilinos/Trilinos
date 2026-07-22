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

	double dx, f, xd;



	rc = 0;


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


	if (!rc) // chatter for cppunit test, which cannot tolerate silence

		printf("OK\n");

	return rc;

	}
