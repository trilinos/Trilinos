// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// Individually check concistency among Hv products via rad2.h, trad2.hpp,
// Fad<Rad> and Rad<Fad> for all unary and binary ops.

#include <cstdio>
#include <float.h>	// for DBL_MAX
#include "Sacado_DisableKokkosCuda.hpp" // Disable Cuda stuff that fails
#include "Sacado_MathFunctions.hpp"
#include "Sacado_trad.hpp"
#include "Sacado_trad2.hpp"
#include "Sacado_rad2.hpp"
#include "Sacado_Fad_SFad.hpp"

typedef Sacado::Rad2d::ADvar ADVar;
typedef Sacado::Rad2::ADvar<double> DADVar;
typedef Sacado::Rad::ADvar<double> DReal;
typedef Sacado::Fad::SFad<DReal,1> FDReal;
typedef Sacado::Fad::SFad<double,1> FReal;
typedef Sacado::Rad::ADvar<FReal> AFReal;

#define UF(T,f,fn) T fn(T &x) { return f(x); }

#define UF4(f) UF(ADVar,f,f ## _ADV) UF(DADVar,f,f ## _DADV) UF(FDReal,f,f ## _FDR) UF(AFReal,f,f ## _AFR)

UF4(abs)
UF4(acos)
UF4(acosh)
UF4(asin)
UF4(asinh)
UF4(atan)
UF4(atanh)
UF4(cos)
UF4(cosh)
UF4(exp)
UF4(fabs)
UF4(log)
UF4(log10)
UF4(sin)
UF4(sinh)
UF4(sqrt)
UF4(tan)
UF4(tanh)

#undef UF
#define UF(T,f,fn) T fn(T &x, T &y) { return f(x,y); }

UF4(atan2)
UF4(pow)

typedef ADVar	(*ADVar_uf )(ADVar  &);
typedef DADVar	(*DADVar_uf)(DADVar &);
typedef FDReal	(*FDReal_uf)(FDReal &);
typedef AFReal	(*AFReal_uf)(AFReal &);

typedef ADVar	(*ADVar_uf2 )(ADVar  &, ADVar  &);
typedef DADVar	(*DADVar_uf2)(DADVar &, DADVar &);
typedef FDReal	(*FDReal_uf2)(FDReal &, FDReal &);
typedef AFReal	(*AFReal_uf2)(AFReal &, AFReal &);

static double
	dom_acosh[]	= { 1., DBL_MAX },
	dom_all[]	= { -DBL_MAX, DBL_MAX },
	dom_invtrig[]	= { -1., 1. },
	dom_log[]	= { DBL_MIN, DBL_MAX };

 struct
Func4 {
	const char *name;
	const double *dom;
	ADVar_uf  f1;
	DADVar_uf f2;
	FDReal_uf f3;
	AFReal_uf f4;
	};

#define U4f(f,d) { #f, d, f ## _ADV, f ## _DADV, f ## _FDR, f ## _AFR }

static Func4 UT[] = {
	U4f(abs,   dom_all),
	U4f(acos,  dom_invtrig),
	U4f(acosh, dom_acosh),
	U4f(asin,  dom_invtrig),
	U4f(asinh, dom_all),
	U4f(atan,  dom_all),
	U4f(atanh, dom_invtrig),
	U4f(cos,   dom_all),
	U4f(cosh,  dom_all),
	U4f(exp,   dom_all),
	U4f(fabs,  dom_all),
	U4f(log,   dom_log),
	U4f(log10, dom_log),
	U4f(sin,   dom_all),
	U4f(sinh,  dom_all),
	U4f(sqrt,  dom_log),
	U4f(tan,   dom_all),
	U4f(tanh,  dom_all)
	};

 struct
Func42 {
	const char *name;
	const double *xdom;
	ADVar_uf2  f1;
	DADVar_uf2 f2;
	FDReal_uf2 f3;
	AFReal_uf2 f4;
	};

static Func42 UT2[] = {
	U4f(atan2, dom_all),
	U4f(pow,   dom_log)
	};

#undef U4f
#undef UF4
#undef UF

 static int
differ(double a, double b)
{
	double d = a - b;;
	if (d == 0.)
		return 0;
	if (d < 0.)
		d = -d;
	if (a < 0.)
		a = -a;
	if (b < 0.)
		b = - b;
	return d > 5e-15 * (a + b);
	}

 static char *progname;

 static int
usage(int rc)
{
	fprintf(rc ? stderr : stdout, "Usage: %s [-v]\n\
	to compare consistency among several schemes for Hessian-vector computations.\n\
	-v ==> verbose (show results)\n\
	No -v ==> just print OK if all goes well; otherwise complain.\n", progname);
	return rc;
	}

 int
main(int argc, char **argv)
{
	Func4 *f, *fe;
	Func42 *f2, *f2e;
	char *s;
	double a, *ap, c, h[4], h1[4], v[3];
	int b0, b1, b2, k, verbose;
	static double unopargs[] = { .7, -.9, 2.5, -4.2, .1, .3, 1.2, -2.4, 0. };
	static double binargs[] = {	.1,.2,	-.3,.4,	-.5,-.6,	.7,-.8,
					1.1,2.2, -2.3,3.4, -2.5,-3.6,	2.7,-3.8, 0.};

#define Dcl(T,T1) T x ## T1, y ## T1, z ## T1
	Dcl(ADVar,  ADV);
	Dcl(DADVar, DADV);
	Dcl(FDReal, FDR);
	Dcl(AFReal, AFR);
#undef	Dcl
	ADVar *pADV[2];
	DADVar *pDADV[2];
	DReal tDR;
	FDReal tfFDR;
	FReal  tfFR;

	progname = argv[0];
	verbose = 0;
	if (argc > 1) {
		if (argc > 2 || *(s = argv[1]) != '-')
			return usage(1);
		switch(s[1]) {
		 case 'h':
		 case '?':
			return usage(s[2] != 0);
		 case '-':
			if (!s[2])
				break;
			return usage(strcmp(s,"--help") ? 1 : 0);
		 case 'v':
			if (!s[2]) {
				verbose = 1;
				break;
				}
		 default:
			return usage(1);
		 }
		}

	v[0] = v[2] = 1.;
	v[1] = 0.;
#define In(T) p ## T[0] = &x ## T; p ## T[1] = &z ## T
	In(ADV);
	In(DADV);
#undef	In
	b0 = b1 = b2 = 0;
	fe = UT + sizeof(UT)/sizeof(UT[0]);
	for(ap = unopargs; (a = *ap++) != 0.; )
	    for(f = UT; f < fe; f++) {
		if (a < f->dom[0] || a > f->dom[1])
			continue;
		xADV = a;
		yADV = f->f1(xADV);
		ADVar::Gradcomp();
		ADVar::Hvprod(1,pADV,v,h);
		if (verbose) {
			printf("%s(%g) = %g\n", f->name, xADV.val(), yADV.val());
			printf("%s' = %g\n", f->name, xADV.adj());
			printf("%s\" = %g\n", f->name, h[0]);
			}
		xDADV = a;
		yDADV = f->f2(xDADV);
		DADVar::Gradcomp();
		if (differ(yDADV.val(), yADV.val())) {
			++b0;
			printf("%s(%g): yADV = %g, yDADV = %g, diff = %.2g\n",
				f->name, a, yADV.val(), yDADV.val(),
				yADV.val() - yDADV.val());
			}
		if (differ(xADV.adj(), xDADV.adj())) {
			++b1;
			printf("%s'(%g): ADV ==> %g, DADV ==> %g, diff = %.2g\n",
				f->name, a, xADV.adj(), xDADV.adj(),
				xADV.adj() - xDADV.adj());
			}
		DADVar::Hvprod(1,pDADV,v,h1);
		if (differ(h[0], h1[0])) {
			++b2;
			printf("%s\"(%g): ADV ==> %g, DADV ==> %g, diff = %.2g\n",
				f->name, a, h[0], h1[0], h[0] - h1[0]);
			}

		tfFDR = 0.;
		tfFDR.diff(0,1);
		xFDR = a + tfFDR*v[0];
		yFDR = f->f3(xFDR);
		tDR = yFDR.dx(0);
		DReal::Gradcomp();

		if (differ(yFDR.val().val(), yADV.val())) {
			++b0;
			printf("%s(%g): yFDR = %g, yADV = %g, diff = %.2g\n",
				f->name, a, yFDR.val().val(), yADV.val(),
				yFDR.val().val() - yADV.val());
			}
		if (differ(xFDR.val().adj(), h[0])) {
			++b2;
			printf("%s\"(%g): FDR ==> %g, ADV ==> %g, diff = %.2g\n",
				f->name, a, xFDR.val().adj(), h[0],
				xFDR.val().adj() - h[0]);
			}

		tfFR = 0.;
		tfFR.diff(0,1);
		xAFR = a + tfFR*v[0];
		yAFR = f->f4(xAFR);
		AFReal::Gradcomp();
		if (differ(yAFR.val().val(), yADV.val())) {
			++b0;
			printf("%s(%g): yAFR = %g, yADV = %g, diff = %.2g\n",
				f->name, a, yAFR.val().val(), yADV.val(),
				yAFR.val().val() - yADV.val());
			}
		if (differ(xAFR.adj().val(), xADV.adj())) {
			++b1;
			printf("%s'(%g): AFR ==> %g, ADV ==> %g, diff = %.2g\n",
				f->name, a, xAFR.adj().val(), xADV.adj(),
				xAFR.adj().val() - xADV.adj());
			}
		if (differ(xAFR.adj().dx(0), h[0])) {
			++b2;
			printf("%s\"(%g): AFR ==> %g, ADV ==> %g, diff = %.2g\n",
				f->name, a, xAFR.adj().dx(0), h[0],
				xAFR.adj().dx(0) - h[0]);
			}
		}
	f2e = UT2 + sizeof(UT2)/sizeof(UT2[0]);
	for(ap = binargs; (a = *ap) != 0.; ap += 2) {
	    c = ap[1];
	    for(f2 = UT2; f2 < f2e; f2++) {
	    	if (a < f2->xdom[0] || a > f2->xdom[1])
			continue;
		xADV = a;
		zADV = c;
		yADV = f2->f1(xADV, zADV);
		ADVar::Gradcomp();
		ADVar::Hvprod(2,pADV,v,h);
		ADVar::Hvprod(2,pADV,v+1,h+2);
		if (verbose) {
			printf("%s(%g,%g) = %g\n", f2->name, xADV.val(), zADV.val(), yADV.val());
			printf("%s' = (%g, %g)\n", f2->name, xADV.adj(), zADV.adj());
			printf("%s\" = (%8g\t%g)\n(%8g\t%g)", f2->name, h[0],h[1],h[2],h[3]);
			}
		xDADV = a;
		zDADV = c;
		yDADV = f2->f2(xDADV, zDADV);
		DADVar::Gradcomp();
		if (differ(yDADV.val(), yADV.val())) {
			++b0;
			printf("%s(%g,%g): yADV = %g, yDADV = %g, diff = %.2g\n",
				f2->name, a, c, yADV.val(), yDADV.val(),
				yADV.val() - yDADV.val());
			}
		if (differ(xADV.adj(), xDADV.adj()) || differ(zADV.adj(), zDADV.adj())) {
			++b1;
			printf("%s'(%g,%g): ADV ==> (%g,%g) DADV ==> (%g,%g), diff = (%.2g,%.2g)\n",
				f2->name, a, c, xADV.adj(), zADV.adj(),
				xDADV.adj(), zDADV.adj(),
				xADV.adj() - xDADV.adj(), zADV.adj() - zDADV.adj());
			}
		DADVar::Hvprod(2,pDADV,v,h1);
		DADVar::Hvprod(2,pDADV,v+1,h1+2);
		for(k = 0; k < 4; k++)
			if (differ(h[k], h1[k])) {
				++b2;
				printf("%s\"(%g,%g):\n ADV ==> (%8g\t%8g\t%8g\t%g)\n",
					f2->name, a, c, h[0], h[1], h[2], h[3]);
				printf("DADV ==> (%8g\t%8g\t%8g\t%g)\n",
					h1[0], h1[1], h1[2], h1[3]);
				printf("diff = (%.2g  %.2g  %.2g  %.2g)\n\n",
					h[0] - h1[0], h[1] - h1[1],
					h[2] - h1[2], h[3] - h1[3]);
				break;
				}
		tfFDR = 0.;
		tfFDR.diff(0,1);
		xFDR = a + tfFDR*v[0];
		zFDR = c + tfFDR*v[1];
		yFDR = f2->f3(xFDR, zFDR);
		tDR = yFDR.dx(0);
		DReal::Gradcomp();

		if (differ(yFDR.val().val(), yADV.val())) {
			++b0;
			printf("%s(%g,%g): yFDR = %g, yADV = %g, diff = %.2g\n",
				f2->name, a, c, yFDR.val().val(), yADV.val(),
				yFDR.val().val() - yADV.val());
			}
		if (differ(xFDR.val().adj(), h[0]) || differ(zFDR.val().adj(), h[1])) {
			++b2;
			printf("%s\"(%g,%g) row 1:\nFDR ==> (%g,%g)\nADV ==> (%g,%g)\n",
				f2->name, a, c, xFDR.val().adj(), zFDR.val().adj(),
				h[0], h[1]);
			printf("diff = %.2g  %g\n\n",
				xFDR.val().adj() - h[0], zFDR.val().adj() - h[1]);
			}
		tfFDR.diff(0,1);
		xFDR = a + tfFDR*v[1];
		zFDR = c + tfFDR*v[2];
		yFDR = f2->f3(xFDR, zFDR);
		tDR = yFDR.dx(0);
		DReal::Gradcomp();
		if (differ(xFDR.val().adj(), h[2]) || differ(zFDR.val().adj(), h[3])) {
			++b2;
			printf("%s\"(%g,%g) row 2:\nFDR ==> (%g,%g)\nADV ==> (%g,%g)\n",
				f2->name, a, c, xFDR.val().adj(), zFDR.val().adj(),
				h[2], h[3]);
			printf("diff = %.2g  %g\n\n",
				xFDR.val().adj() - h[2], zFDR.val().adj() - h[3]);
			}

		tfFR = 0.;
		tfFR.diff(0,1);
		xAFR = a + tfFR*v[0];
		zAFR = c + tfFR*v[1];
		yAFR = f2->f4(xAFR, zAFR);
		AFReal::Gradcomp();
		if (differ(yAFR.val().val(), yADV.val())) {
			++b0;
			printf("%s(%g,%g): yAFR = %g, yADV = %g, diff = %.2g\n",
				f2->name, a, c, yAFR.val().val(), yADV.val(),
				yAFR.val().val() - yADV.val());
			}
		if (differ(xAFR.adj().val(), xADV.adj()) || differ(zAFR.adj().val(), zADV.adj())) {
			++b1;
			printf("%s'(%g,%g):\n AFR ==> (%g,%g)\n ADV ==> (%g,%g)\n",
				f2->name, a, c, xAFR.adj().val(), zAFR.adj().val(),
				xADV.adj(), zADV.adj());
			printf(" diff = %.2g  %.2g\n\n",
				xAFR.adj().val() - xADV.adj(),
				zAFR.adj().val() - zADV.adj());
			}
		if (differ(xAFR.adj().dx(0), h[0]) || differ(zAFR.adj().dx(0), h[1])) {
			++b2;
			printf("%s\"(%g,%g) row 1:\n AFR ==> (%g,%g)\n",
				f2->name, a, c, xAFR.adj().dx(0), zAFR.adj().dx(0));
			printf(" ADV = (%g,%g)\n diff = %.2g  %.2g\n\n", h[0], h[1],
				xAFR.adj().dx(0) - h[0],
				zAFR.adj().dx(0) - h[1]);
			}
		tfFR.diff(0,1);
		xAFR = a + tfFR*v[1];
		zAFR = c + tfFR*v[2];
		yAFR = f2->f4(xAFR, zAFR);
		AFReal::Gradcomp();
		if (differ(xAFR.adj().dx(0), h[2]) || differ(zAFR.adj().dx(0), h[3])) {
			++b2;
			printf("%s\"(%g,%g) row 2:\n AFR ==> (%g,%g)\n",
				f2->name, a, c, xAFR.adj().dx(0), zAFR.adj().dx(0));
			printf(" ADV = (%g,%g)\n diff = %.2g  %.2g\n\n", h[2], h[3],
				xAFR.adj().dx(0) - h[2],
				zAFR.adj().dx(0) - h[3]);
			}
		}
	    }
	k = b0 + b1 + b2;
	if (k || verbose) {
		s = progname;
		if (*s == '.' && s[1] == '/')
			for(s += 2; *s == '/'; ++s);
		printf("%s: %d errors seen.\n", s, k);
		k = 1;
		}
	else
		printf("OK\n");
	return k;
	}
