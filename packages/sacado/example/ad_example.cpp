// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// Brief demo of Fad and Rad for computing first derivatives

#include <Sacado_No_Kokkos.hpp>		// for FAD and RAD
#include <cstdio>		// nicer than streams in some respects
using std::printf;

using namespace std;

// Typedefs reduce gobbledygook below

typedef Sacado::Fad::DFad<double>   F;  // FAD with # of ind. vars given later
typedef Sacado::Fad::SFad<double,2> F2; // FAD with # of ind. vars fixed at 2
typedef Sacado::Rad::ADvar<double>  R;  // for RAD

template <typename T>
const T func2(T &a, T &b)	// sample function of 2 variables
{ return sqrt(a*a + 2*b*b); }

template <typename T>
const T func(int n, T *x)	// sample function of n variables
				// == func2 when n == 2
{
	int i;
	T t = 0;
	for(i = 1; i < n; i++)
		t += i*x[i]*x[i];
	return sqrt(t);
	}

// demo of FAD (forward-mode AD), with general number of ind. vars

 void
Fad_demo()
{
	F a, b, x[5], y;
	int i, n;

	printf("Fad_demo...\n\n");

	// first try n == 2
	a = 1.;
	b = 2.;
	// indicate the independent variables, and initialize their partials to 1:
	a.diff(0,2);	// 0 ==> this is the first independent var., of 2
	b.diff(1,2);	// 1 ==> this is the second ind. var.
	
	y = func2(a,b);

	printf("func2(%g,%g) = %g\n", a.val(), b.val(), y.val());

	printf("partials of func2 = %g, %g\n", y.dx(0), y.dx(1));

	// When we know the result was not constant (i.e., did involve ind. vars)
	// or when hasFastAccess() is true, we access partials more quickly
	// by using member function fastAccessDx rather than dx

	if (y.hasFastAccess())
		printf("Repeat with fastAccess: partials of func2 = %g, %g\n",
			y.fastAccessDx(0), y.fastAccessDx(1));

	// Similar exercise with general n, in this case n == 5
	n = 5;
	for(i = 0; i < n; i++) {
		x[i] = i;
		x[i].diff(i, n);
		}
	y = func(n, x);
	printf("\nfunc(5,x) for x = (0,1,2,3,4) = %g\n", y.val());
	for(i = 0; i < n; i++)
		printf("d func / d x[%d] = %g == %g\n", i, y.dx(i), y.fastAccessDx(i));
	}

// Fad_demo2 == repeat first part of Fad_Demo with type F2 instead of F
// i.e., with fixed-size allocations

 void
Fad2_demo()
{
	F2 a, b, y;

	printf("\n\nFad2_demo...\n\n");

	a = 1.;
	b = 2.;
	// indicate the independent variables, and initialize their partials to 1:
	a.diff(0,2);	// 0 ==> this is the first independent var., of 2
	b.diff(1,2);	// 1 ==> this is the second ind. var.
	
	y = func2(a,b);

	printf("func2(%g,%g) = %g\n", a.val(), b.val(), y.val());

	printf("partials of func2 = %g, %g\n", y.dx(0), y.dx(1));

	if (y.hasFastAccess())
		printf("Repeat with fastAccess: partials of func2 = %g, %g\n",
			y.fastAccessDx(0), y.fastAccessDx(1));
	}

// Fad_demo3 == repeat of Fad_Demo2 with a different constructor, one that
// indicates the independent variables and their initial values
// and removes the need to invoke .diff()

 void
Fad3_demo()
{
	F2 a(2,0,1.), b(2,1,2.), y;

	printf("\n\nFad3_demo...\n\n");

	y = func2(a,b);

	printf("func2(%g,%g) = %g\n", a.val(), b.val(), y.val());

	printf("partials of func2 = %g, %g\n", y.dx(0), y.dx(1));

	if (y.hasFastAccess())
		printf("Repeat with fastAccess: partials of func2 = %g, %g\n",
			y.fastAccessDx(0), y.fastAccessDx(1));
	}

// Rad_demo == repeat of Fad_Demo with type R instead of F,
// i.e., with reverse-mode rather than forward-mode AD

 void
Rad_demo()
{
	R a, b, x[5], y;
	int i, n;

	printf("\n\nRad_demo...\n\n");

	// first try n == 2
	a = 1.;
	b = 2.;
	
	y = func2(a,b);

	R::Gradcomp();	// do the reverse sweep

	printf("func2(%g,%g) = %g\n", a.val(), b.val(), y.val());

	printf("partials of func2 = %g, %g\n", a.adj(), b.adj());

	// Similar exercise with general n, in this case n == 5
	n = 5;
	for(i = 0; i < n; i++)
		x[i] = i;
	y = func(n, x);
	printf("\nfunc(5,x) for x = (0,1,2,3,4) = %g\n", y.val());

	// the .val() values are always available; we must call Gradcomp
	// before accessing the adjoints, i.e., the .adj() values...

	R::Gradcomp();
	
	for(i = 0; i < n; i++)
		printf("d func / d x[%d] = %g\n", i, x[i].adj());
	}

 int
main()
{
	Fad_demo();
	Fad2_demo();
	Fad3_demo();
	Rad_demo();
	return 0;
	}
