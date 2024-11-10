// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// tradvec_example
//
//  usage:
//     tradvec_example
//
//  output:
//     prints the results of differentiating two simple functions and their
//     sum with reverse mode AD using functions Outvar_Gradcomp and
//     Weighted_GradcompVec in the Sacado::RadVec::ADvar class.

/* Simple test of Outvar_Gradcomp and Weighted_GradcompVec. */

#include "Sacado_tradvec.hpp"
#include <stdio.h>

#ifdef _MSC_VER
# define snprintf _snprintf
#endif

typedef Sacado::RadVec::ADvar<double> ADVar;

 ADVar
foo(double d, ADVar x, ADVar y)
{ return d*x*y; }

 ADVar
goo(double d, ADVar x, ADVar y)
{ return d*(x*x + y*y); }

 struct
ExpectedAnswer {
	const char *name;
	double v;	// value of name
	double dvdx;	// partial w.r.t. x
	double dvdy;	// partial w.r.t. y
	};

 static ExpectedAnswer expected[4] = {
	{ "a", 6., 3., 2. },
	{ "b", 13., 4., 6. },
	{ "c", 19., 7., 8. },
	{ "(a + b + c)", 38., 14., 16. }};

 int
botch(ExpectedAnswer *e, const char *partial, double got, double wanted)
{
	char buf[32];
	const char *what;

	what = e->name;
	if (partial) {
		snprintf(buf, sizeof(buf), "d%s/d%s", what, partial);
		what = buf;
		}
	fprintf(stderr, "Expected %s = %g, but got %g\n", what, wanted, got);
	return 1;
	}

 int
acheck(int k, double d, double v, double dvdx, double dvdy)
{
	ExpectedAnswer *e = &expected[k];
	int nbad = 0;

	/* There should be no round-off error in this simple example, so we */
	/* use exact comparisons, rather than, say, relative differences. */

	if (v != d*e->v)
		nbad += botch(e, 0, v, d*e->v);
	if (dvdx != d*e->dvdx)
		nbad += botch(e, "x", dvdx, d*e->dvdx);
	if (dvdy != d*e->dvdy)
		nbad += botch(e, "y", dvdy, d*e->dvdy);
	return nbad;
	}

 int
main(void)
{
	double d, z[4];
	int i, nbad;

	static ADVar a, b, c, x, y, *v[3] = {&a, &b, &c};
	static ADVar **V[4] = {v, v+1, v+2, v};
	static size_t np[4] = {1, 1, 1, 3};
	static double w[3] = { 1., 1., 1. };
	static double *W[4] = {w, w, w, w};

	nbad = 0;
	for(d = 1.; d <= 2.; ++d) {
		printf("\nd = %g\n", d);
		x = 2;
		y = 3;
		a = foo(d,x,y);
		b = goo(d,x,y);
		c = a + b;

		ADVar::Outvar_Gradcomp(a);
		printf("a = %g\n", a.val());
		printf("da/dx = %g\n", x.adj());
		printf("da/dy = %g\n", y.adj());
		nbad += acheck(0, d, a.val(), x.adj(), y.adj());
		z[0] = a.val();

		ADVar::Outvar_Gradcomp(b);
		printf("b = %g\n", b.val());
		printf("db/dx = %g\n", x.adj());
		printf("db/dy = %g\n", y.adj());
		nbad += acheck(1, d, b.val(), x.adj(), y.adj());
		z[1] = b.val();

		ADVar::Outvar_Gradcomp(c);
		printf("c = %g (should be a + b)\n", c.val());
		printf("dc/dx = %g\n", x.adj());
		printf("dc/dy = %g\n", y.adj());
		nbad += acheck(2, d, c.val(), x.adj(), y.adj());
		z[2] = c.val();
		z[3] = z[0] + z[1] + z[2];

		ADVar::Weighted_GradcompVec(4,np,V,W);
		for(i = 0; i < 4; ++i) {
			printf("w %d:\td/dx = %g\td/dy = %g\n", i, x.adj(i), y.adj(i));
			nbad += acheck(i, d, z[i], x.adj(i), y.adj(i));
			}
		}
	if (nbad == 0)
		printf("\nExample passed!\n");
	else
		printf("\nSomething is wrong, example failed!\n");
	return nbad > 0;
	}
