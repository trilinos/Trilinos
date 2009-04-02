#!/bin/sh
# Lines ending in \n\ work around a bug in gawk
awk 'BEGIN { ftop = "\n\
/* Try to test all combinations of types and operations */\n\
\n\
#define ADT_RAD Sacado::Rad::\n\
\n\
#include \"Sacado_trad.hpp\"\n\
#include <cstdio>\n\
using std::printf;\n\
\n\
typedef ADT_RAD IndepADvar<double> AI;\n\
typedef ADT_RAD ADvar<double> A;\n\
typedef ADT_RAD ConstADvar<double> C;\n\
typedef ADT_RAD ADvari<double> Ai;\n\
typedef const ADT_RAD IndepADvar<double> cAI;\n\
typedef const ADT_RAD ADvar<double> cA;\n\
typedef const ADT_RAD ConstADvar<double> cC;\n\
typedef const ADT_RAD ADvari<double> cAi;\n\
static int rc;\n\
\n\
/* This is to be run through an awk program that changes lines */\n\
/* with \"BINTEST\" or \"UNOPTEST\" at the beginning of the line into */\n\
/* a the desired C++ (which we can then inspect). */\n\
\n\
 void\n\
botch(const char *what, double wanted, double got)\n\
{\n\
	printf(\"%s: expected %g, got %g, diff = %.2g\\n\", what, wanted, got, wanted-got);\n\
	rc = 1;\n\
	}\n\
\n\
 const double tol = 5e-16;\n\
\n\
 int\n\
differ(double a, double b)\n\
{\n\
	double d = a - b;\n\
	if (d < 0.)\n\
		d = -d;\n\
	if (a < 0.)\n\
		a = -a;\n\
	if (b < 0.)\n\
		b = -b;\n\
	if (a < b)\n\
		a = b;\n\
	if (a > 0.)\n\
		d /= a;\n\
	return d > tol;\n\
	}\n\
\n\
#ifndef RAD_EQ_ALIAS\n\
#define Plus_dx 1.\n\
#else\n\
#ifdef RAD_AUTO_AD_Const\n\
#define Plus_dx 1.\n\
#else\n\
#define Plus_dx 0.\n\
#endif\n\
#endif\n\
\n\
 int\n\
main(void)\n\
{\n\
	AI xAI, yAI;\n\
	A fA, xA, yA;\n\
	C xC, yC;\n\
	double dx, dy, f, xd, yd;\n\
	long xL, yL;\n\
	int xi, yi;\n\
\n\
	rc = 0;\n\
"
fbot = "\n\
	if (!rc) // chatter for cppunit test, which cannot tolerate silence\n\
		printf(\"OK\\n\");\n\
	return rc;\n\
	}"
suf[0] = "AI"
suf[1] = "A"
suf[2] = "C"
suf[3] = "Ai"
suf[4] = "d"
suf[5] = "L"
suf[6] = "i"
cast[0] = cast[1] = cast[2] = cast[3] = cast[4] = ""
cast[5] = "(long)"
cast[6] = "(int)"
adj[0] = adj[1] = adj[2] = "adj()"
adj[3] = "aval"
}
function crunch(op,i,j,cL,cR)
{
	L = suf[i]
	R = suf[j]
	newcontext = 0
	if (cL + cR || i ==3 || j == 3) {
		newcontext = 1
		printf "\t{\n" >>fname
		if (cL + cR) printf "\tA::aval_reset();\n" >>fname
		}
	if (cL) {
		x = "xc" L
		printf "\tc%s %s(xd);\n", L, x >>fname
		}
	else {
		x = "x" L
		if (i == 3) printf "\t%s %s(xd);\n", L, x >>fname
		else if (i != 4) printf "\t%s = %sxd;\n", x, cast[i] >>fname
		}
	if (cR) {
		y = "yc" R
		printf "\tc%s %s(yd);\n", R, y >>fname
		}
	else {
		y = "y" R
		if (j == 3) printf "\t%s %s(yd);\n", R, y >>fname
		else if (j != 4) printf "\t%s = %syd;\n", y, cast[j] >>fname
		}
	if (infix) opname = sprintf("%s %s %s", x, op, y)
	else	   opname = sprintf("%s(%s,%s)", op, x, y)
	printf "\tfA = %s;\n\tA::Gradcomp();\n", opname >>fname
	printf "\tif (differ(fA.val(), f)) botch(\"fA = %s\", f, fA.val());\n",opname >>fname
	if (i < 4)
	 printf "\telse if (differ(%s.%s, dx)) botch(\"d %s/dx\", dx, %s.%s);\n",x,adj[i],opname,x,adj[i] >>fname
	if (j < 4)
	 printf "\telse if (differ(%s.%s, dy)) botch(\"d %s/dy\", dy, %s.%s);\n",y,adj[j],opname,y,adj[j] >>fname
	if (newcontext) printf "\t}\n" >>fname
	}

function crunch1(op, i, cL)
{
	L = suf[i]
	newcontext = 0
	if (cL || i == 3) {
		newcontext = 1
		printf "\t{\n" >>fname
		if (cL) printf "\tA::aval_reset();\n" >>fname
		x = "xc" L
		printf "\tc%s %s(xd);\n", L, x >>fname
		}
	else {
		x = "x" L
		if (i == 3) printf "\t%s %s(xd);\n", L, x >>fname
		else if (i != 4) printf "\t%s = %sxd;\n", x, cast[i] >>fname
		}
	xx = x
	if (copy) xx = sprintf("copy(%s)", x)
	if (infix) opname = sprintf("%s%s", op, xx)
	else	   opname = sprintf("%s(%s)", op, xx)
	printf "\tfA = %s;\n\tA::Gradcomp();\n", opname >>fname
	printf "\tif (differ(fA.val(), f)) botch(\"fA = %s\", f, fA.val());\n",opname >>fname
	printf "\telse if (differ(%s.%s, dx)) botch(\"d %s/dx\", dx, %s.%s);\n",x,adj[i],opname,x,adj[i] >>fname
	if (newcontext) printf "\t}\n" >>fname
	}

/^BINTEST/ {
	fname = sprintf("tradoptest_%02d.cpp",NR)
	op = $2
	infix = $8 + 0
	printf "%s\n\t/**** Test of %s ****/\n\n", ftop, op >fname
	printf "\txd = %s; yd = %s; f = %s; dx = %s; dy = %s;\n",$3,$4,$5,$6,$7 >>fname
	for(i = 0; i < 7; i++) {
		jx = 7; if (i > 3) jx = 4
		for(j = 0; j < jx; j++) {
			crunch(op,i,j,0,0)
			if (i < 4) crunch(op,i,j,1,0)
			if (j < 4) {
				crunch(op,i,j,0,1)
				if (i < 3) crunch(op,i,j,1,1)
				}
			}
		}
	printf "%s\n", fbot >>fname
	close(fname)
	next
	}
/^UNOPTEST/ {
	fname = sprintf("tradoptest_%02d.cpp",NR)
	op = $2
	infix = $6 + 0
	copy = $7 + 0
	printf "%s\n\t/**** Test of %s ****/\n\n", ftop, op >fname
	printf "\txd = %s; f = %s; dx = %s;\n",$3,$4,$5 >>fname
	for(i = 0; i < 4; i++) {
		crunch1(op,i,0)
		crunch1(op,i,1)
		}
	printf "%s\n", fbot >>fname
	close(fname)
	next
	}
{print}' "$@"
