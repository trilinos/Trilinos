/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

/* az_reorder.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
        on Microsoft Windows system, link with libf2c.lib;
        on Linux or Unix systems, link with .../path/to/libf2c.a -lm
        or, if you install libf2c.a in a standard place, with -lf2c -lm
        -- in that order, at the end of the command line, as in
                cc *.o -lf2c -lm
        Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

                http://www.netlib.org/f2c/libf2c.zip
*/

/* From f2c.h  --  Standard Fortran to C header file */

/* 1/15/2008 JW changed the first typedef from typedef long int integer; to
   typedef int integer; to prevent two tests from segfaulting */
typedef int integer;
typedef unsigned long int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
#ifdef INTEGER_STAR_8	/* Adjust for integer*8. */
typedef long long longint;		/* system-dependent */
typedef unsigned long long ulongint;	/* system-dependent */
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

/*typedef long int Long;*/	/* No longer used; formerly in Namelist */

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;	/* complex function */
typedef VOID H_f;	/* character function */
typedef VOID Z_f;	/* double complex function */
typedef doublereal E_f;	/* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif

/*     This stuff comes from George & Liu. It basically corresponds */
/*     to computing a Reverse Cuthill-McKee Ordering corresponding */
/*     to the graph of a matrix. */

/* Subroutine */ int az_rcm_c(integer *root, integer *xadj, integer *adjncy, 
	integer *mask, integer *perm, integer *ccsize, integer *deg)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l;
    extern /* Subroutine */ int az_degree_c(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer nbr, node, fnbr, lnbr, lperm, jstop, jstrt, lbegin, lvlend;

/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */

    /* Parameter adjustments */
    --deg;
    --perm;
    --mask;
    --adjncy;
    --xadj;

    /* Function Body */
    az_degree_c(root, &xadj[1], &adjncy[1], &mask[1], &deg[1], ccsize, &perm[
	    1]);

    mask[*root] = 0;
    if (*ccsize <= 1) {
	return 0;
    }
    lvlend = 0;
    lnbr = 1;

L100:
    lbegin = lvlend + 1;

    lvlend = lnbr;

    i__1 = lvlend;
    for (i__ = lbegin; i__ <= i__1; ++i__) {
	node = perm[i__];
	jstrt = xadj[node];
	jstop = xadj[node + 1] - 1;
	fnbr = lnbr + 1;

	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    nbr = adjncy[j];
	    if (mask[nbr] == 0) {
		goto L200;
	    }
	    ++lnbr;
	    mask[nbr] = 0;
	    perm[lnbr] = nbr;
L200:
	    ;
	}

	if (fnbr >= lnbr) {
	    goto L600;
	}
	k = fnbr;

L300:
	l = k;

	++k;
	nbr = perm[k];

L400:
	if (l < fnbr) {
	    goto L500;
	}

	lperm = perm[l];
	if (deg[lperm] <= deg[nbr]) {
	    goto L500;
	}
	perm[l + 1] = lperm;
	--l;
	goto L400;

L500:
	perm[l + 1] = nbr;

	if (k < lnbr) {
	    goto L300;
	}
L600:
	;
    }

    if (lnbr > lvlend) {
	goto L100;
    }
    k = *ccsize / 2;
    l = *ccsize;

    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lperm = perm[l];
	perm[l] = perm[i__];
	perm[i__] = lperm;
	--l;
/* L700: */
    }


    return 0;
} /* az_rcm*/

/* Subroutine */ int az_degree_c(integer *root, integer *xadj, integer *
	adjncy, integer *mask, integer *deg, integer *ccsize, integer *ls)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, nbr, ideg, node, jstop, jstrt, lbegin, lvlend, 
	    lvsize;

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */

    /* Parameter adjustments */
    --ls;
    --deg;
    --mask;
    --adjncy;
    --xadj;

    /* Function Body */
    ls[1] = *root;
    xadj[*root] = -xadj[*root];
    lvlend = 0;
    *ccsize = 1;

L100:
    lbegin = lvlend + 1;

    lvlend = *ccsize;
    i__1 = lvlend;
    for (i__ = lbegin; i__ <= i__1; ++i__) {
	node = ls[i__];
	jstrt = -xadj[node];
	jstop = (i__2 = xadj[node + 1], abs(i__2)) - 1;
	ideg = 0;
	if (jstop < jstrt) {
	    goto L300;
	}

	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    nbr = adjncy[j];
	    if (mask[nbr] == 0) {
		goto L200;
	    }
	    ++ideg;
	    if (xadj[nbr] < 0) {
		goto L200;
	    }
	    xadj[nbr] = -xadj[nbr];
	    ++(*ccsize);
	    ls[*ccsize] = nbr;
L200:
	    ;
	}

L300:
	deg[node] = ideg;

/* L400: */
    }

    lvsize = *ccsize - lvlend;
    if (lvsize > 0) {
	goto L100;
    }

    i__1 = *ccsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	node = ls[i__];
	xadj[node] = -xadj[node];
/* L500: */
    }


    return 0;
} /* az_degree_c */

/* Subroutine */ int az_fnroot_c(integer *root, integer *xadj, integer *
	adjncy, integer *mask, integer *nlvl, integer *xls, integer *ls)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, k;
    extern /* Subroutine */ int az_rootls_c(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer ndeg, node, nabor, kstop, jstrt, kstrt, mindeg, ccsize, 
	    nunlvl;

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ls;
    --xls;
    --mask;
    --adjncy;
    --xadj;

    /* Function Body */
    az_rootls_c(root, &xadj[1], &adjncy[1], &mask[1], nlvl, &xls[1], &ls[1]);

    ccsize = xls[*nlvl + 1] - 1;
    if (*nlvl == 1 || *nlvl == ccsize) {
	return 0;
    }

L100:
    jstrt = xls[*nlvl];

    mindeg = ccsize;
    *root = ls[jstrt];
    if (ccsize == jstrt) {
	goto L400;
    }

    i__1 = ccsize;
    for (j = jstrt; j <= i__1; ++j) {
	node = ls[j];
	ndeg = 0;
	kstrt = xadj[node];
	kstop = xadj[node + 1] - 1;

	i__2 = kstop;
	for (k = kstrt; k <= i__2; ++k) {
	    nabor = adjncy[k];
	    if (mask[nabor] > 0) {
		++ndeg;
	    }
/* L200: */
	}

	if (ndeg >= mindeg) {
	    goto L300;
	}
	*root = node;
	mindeg = ndeg;
L300:
	;
    }

L400:
    az_rootls_c(root, &xadj[1], &adjncy[1], &mask[1], &nunlvl, &xls[1], &ls[1]
	    );

    if (nunlvl <= *nlvl) {
	return 0;
    }
    *nlvl = nunlvl;
    if (*nlvl < ccsize) {
	goto L100;
    }


    return 0;
} /* az_fnroot*/

/* Subroutine */ int az_rootls_c(integer *root, integer *xadj, integer *
	adjncy, integer *mask, integer *nlvl, integer *xls, integer *ls)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, nbr, node, jstop, jstrt, lbegin, ccsize, lvlend, 
	    lvsize;

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ls;
    --xls;
    --mask;
    --adjncy;
    --xadj;

    /* Function Body */
    mask[*root] = 0;
    ls[1] = *root;
    *nlvl = 0;
    lvlend = 0;
    ccsize = 1;

L200:
    lbegin = lvlend + 1;

    lvlend = ccsize;
    ++(*nlvl);
    xls[*nlvl] = lbegin;

    i__1 = lvlend;
    for (i__ = lbegin; i__ <= i__1; ++i__) {
	node = ls[i__];
	jstrt = xadj[node];
	jstop = xadj[node + 1] - 1;
	if (jstop < jstrt) {
	    goto L400;
	}

	i__2 = jstop;
	for (j = jstrt; j <= i__2; ++j) {
	    nbr = adjncy[j];
	    if (mask[nbr] == 0) {
		goto L300;
	    }
	    ++ccsize;
	    ls[ccsize] = nbr;
	    mask[nbr] = 0;
L300:
	    ;
	}

L400:
	;
    }
    lvsize = ccsize - lvlend;
    if (lvsize > 0) {
	goto L200;
    }
    xls[*nlvl + 1] = lvlend + 1;

    i__1 = ccsize;
    for (i__ = 1; i__ <= i__1; ++i__) {
	node = ls[i__];
	mask[node] = 1;
/* L500: */
    }


    return 0;
} /* az_rootls */

