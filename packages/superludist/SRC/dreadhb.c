#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "superlu_ddefs.h"

/*
 * Prototypes
 */
static void ReadVector(FILE *, int_t, int_t *, int_t, int_t);
static void dReadValues(FILE *, int_t, double *, int_t, int_t);
static void FormFullA(int_t, int_t *, double **, int_t **, int_t **);
static int dDumpLine(FILE *);
static int dParseIntFormat(char *, int_t *, int_t *);
static int dParseFloatFormat(char *, int_t *, int_t *);


void
dreadhb_dist(int iam, FILE *fp, int_t *nrow, int_t *ncol, int_t *nonz,
	     double **nzval, int_t **rowind, int_t **colptr)
{
/* 
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 *
 * Purpose
 * =======
 * 
 * Read a DOUBLE PRECISION matrix stored in Harwell-Boeing format 
 * as described below.
 * 
 * Line 1 (A72,A8) 
 *  	Col. 1 - 72   Title (TITLE) 
 *	Col. 73 - 80  Key (KEY) 
 * 
 * Line 2 (5I14) 
 * 	Col. 1 - 14   Total number of lines excluding header (TOTCRD) 
 * 	Col. 15 - 28  Number of lines for pointers (PTRCRD) 
 * 	Col. 29 - 42  Number of lines for row (or variable) indices (INDCRD) 
 * 	Col. 43 - 56  Number of lines for numerical values (VALCRD) 
 *	Col. 57 - 70  Number of lines for right-hand sides (RHSCRD) 
 *                    (including starting guesses and solution vectors 
 *		       if present) 
 *           	      (zero indicates no right-hand side data is present) 
 *
 * Line 3 (A3, 11X, 4I14) 
 *   	Col. 1 - 3    Matrix type (see below) (MXTYPE) 
 * 	Col. 15 - 28  Number of rows (or variables) (NROW) 
 * 	Col. 29 - 42  Number of columns (or elements) (NCOL) 
 *	Col. 43 - 56  Number of row (or variable) indices (NNZERO) 
 *	              (equal to number of entries for assembled matrices) 
 * 	Col. 57 - 70  Number of elemental matrix entries (NELTVL) 
 *	              (zero in the case of assembled matrices) 
 * Line 4 (2A16, 2A20) 
 * 	Col. 1 - 16   Format for pointers (PTRFMT) 
 *	Col. 17 - 32  Format for row (or variable) indices (INDFMT) 
 *	Col. 33 - 52  Format for numerical values of coefficient matrix (VALFMT) 
 * 	Col. 53 - 72 Format for numerical values of right-hand sides (RHSFMT) 
 *
 * Line 5 (A3, 11X, 2I14) Only present if there are right-hand sides present 
 *    	Col. 1 	      Right-hand side type: 
 *	         	  F for full storage or M for same format as matrix 
 *    	Col. 2        G if a starting vector(s) (Guess) is supplied. (RHSTYP) 
 *    	Col. 3        X if an exact solution vector(s) is supplied. 
 *	Col. 15 - 28  Number of right-hand sides (NRHS) 
 *	Col. 29 - 42  Number of row indices (NRHSIX) 
 *          	      (ignored in case of unassembled matrices) 
 *
 * The three character type field on line 3 describes the matrix type. 
 * The following table lists the permitted values for each of the three 
 * characters. As an example of the type field, RSA denotes that the matrix 
 * is real, symmetric, and assembled. 
 *
 * First Character: 
 *	R Real matrix 
 *	C Complex matrix 
 *	P Pattern only (no numerical values supplied) 
 *
 * Second Character: 
 *	S Symmetric 
 *	U Unsymmetric 
 *	H Hermitian 
 *	Z Skew symmetric 
 *	R Rectangular 
 *
 * Third Character: 
 *	A Assembled 
 *	E Elemental matrices (unassembled) 
 *
 */

    register int_t i, numer_lines, rhscrd = 0;
    int_t tmp, colnum, colsize, rownum, rowsize, valnum, valsize;
    char buf[100], type[4];
    int_t sym;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(0, "Enter dreadhb_dist()");
#endif

    /* Line 1 */
    fgets(buf, 100, fp);
    /*dDumpLine(fp);*/
    /*if ( !iam ) fflush(stdout);*/

    /* Line 2 */
    for (i=0; i<5; i++) {
	fscanf(fp, "%14c", buf); buf[14] = 0;
	tmp = atoi(buf); /*sscanf(buf, "%d", &tmp);*/
	if (i == 3) numer_lines = tmp;
	if (i == 4 && tmp) rhscrd = tmp;
    }
    dDumpLine(fp);

    /* Line 3 */
    fscanf(fp, "%3c", type);
    fscanf(fp, "%11c", buf); /* pad */
    type[3] = 0;
#if ( DEBUGlevel>=1 )
    if ( !iam ) printf("Matrix type %s\n", type);
#endif
    
    fscanf(fp, "%14c", buf); *nrow = atoi(buf); /*sscanf(buf, "%d", nrow);*/
    fscanf(fp, "%14c", buf); *ncol = atoi(buf); /*sscanf(buf, "%d", ncol);*/
    fscanf(fp, "%14c", buf); *nonz = atoi(buf); /*sscanf(buf, "%d", nonz);*/
    fscanf(fp, "%14c", buf); tmp = atoi(buf);   /*sscanf(buf, "%d", &tmp);*/
    
    if (tmp != 0)
	if ( !iam ) printf("This is not an assembled matrix!\n");
    if (*nrow != *ncol)
	if ( !iam ) printf("Matrix is not square.\n");
    dDumpLine(fp);

    /* Allocate storage for the three arrays ( nzval, rowind, colptr ) */
    dallocateA_dist(*ncol, *nonz, nzval, rowind, colptr);

    /* Line 4: format statement */
    fscanf(fp, "%16c", buf);
    dParseIntFormat(buf, &colnum, &colsize);
    fscanf(fp, "%16c", buf);
    dParseIntFormat(buf, &rownum, &rowsize);
    fscanf(fp, "%20c", buf);
    dParseFloatFormat(buf, &valnum, &valsize);
    fscanf(fp, "%20c", buf);
    dDumpLine(fp);

    /* Line 5: right-hand side */    
    if ( rhscrd ) dDumpLine(fp); /* skip RHSFMT */

#if ( DEBUGlevel>=1 )
    if ( !iam ) {
	printf("%d rows, %d nonzeros\n", *nrow, *nonz);
	printf("colnum %d, colsize %d\n", colnum, colsize);
	printf("rownum %d, rowsize %d\n", rownum, rowsize);
	printf("valnum %d, valsize %d\n", valnum, valsize);
    }
#endif
    
    ReadVector(fp, *ncol+1, *colptr, colnum, colsize);
#if ( DEBUGlevel>=1 )
    if ( !iam )	printf("read colptr[%d] = %d\n", *ncol, (*colptr)[*ncol]);
#endif
    ReadVector(fp, *nonz, *rowind, rownum, rowsize);
#if ( DEBUGlevel>=1 )
    if ( !iam )	printf("read rowind[%d] = %d\n", *nonz-1, (*rowind)[*nonz-1]);
#endif
    if ( numer_lines ) {
        dReadValues(fp, *nonz, *nzval, valnum, valsize);
#if ( DEBUGlevel>=1 )
	if ( !iam ) printf("read nzval[%d] = %e\n", *nonz-1, (*nzval)[*nonz-1]);
#endif
    }

    sym = (type[1] == 'S' || type[1] == 's');
    if ( sym ) {
	FormFullA(*ncol, nonz, nzval, rowind, colptr);
    }

    /*if ( !iam ) fflush(stdout);*/
    fclose(fp);
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(0, "Exit dreadhb_dist()");
#endif
}

/* Eat up the rest of the current line */
static int dDumpLine(FILE *fp)
{
    register int c;
    while ((c = fgetc(fp)) != '\n') ;
    return 0;
}

static int dParseIntFormat(char *buf, int_t *num, int_t *size)
{
    char *tmp;

    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
    while (*tmp != 'I' && *tmp != 'i') ++tmp;
    ++tmp;
    *size = atoi(tmp); /*sscanf(tmp, "%d", size);*/
    return 0;
}

static int dParseFloatFormat(char *buf, int_t *num, int_t *size)
{
    char *tmp, *period;
    
    tmp = buf;
    while (*tmp++ != '(') ;
    *num = atoi(tmp); /*sscanf(tmp, "%d", num);*/
    while (*tmp != 'E' && *tmp != 'e' && *tmp != 'D' && *tmp != 'd'
	   && *tmp != 'F' && *tmp != 'f') {
       /* May find kP before nE/nD/nF, like (1P6F13.6). In this case the
           num picked up refers to P, which should be skipped. */
        if (*tmp=='p' || *tmp=='P') {
           ++tmp;
           sscanf(tmp, "%d", num);
        } else {
           ++tmp;
        }
    }
    ++tmp;
    period = tmp;
    while (*period != '.' && *period != ')') ++period ;
    *period = '\0';
    *size = atoi(tmp); /*sscanf(tmp, "%2d", size);*/

    return 0;
}

static void
ReadVector(FILE *fp, int_t n, int_t *where, int_t perline, int_t persize)
{
    register int_t i, j, item;
    char tmp, buf[100];
    
    i = 0;
    while (i < n) {
	fgets(buf, 100, fp);    /* read a line at a time */
	for (j=0; j<perline && i<n; j++) {
	    tmp = buf[(j+1)*persize];     /* save the char at that place */
	    buf[(j+1)*persize] = 0;       /* null terminate */
	    item = atoi(&buf[j*persize]); 
	    buf[(j+1)*persize] = tmp;     /* recover the char at that place */
	    where[i++] = item - 1;
	}
    }
}

static void
dReadValues(FILE *fp, int_t n, double *destination, int_t perline, int_t persize)
{
    register int_t i, j, k, s;
    char tmp, buf[100];
    
    i = 0;
    while (i < n) {
	fgets(buf, 100, fp);    /* read a line at a time */
	for (j=0; j<perline && i<n; j++) {
	    tmp = buf[(j+1)*persize];     /* save the char at that place */
	    buf[(j+1)*persize] = 0;       /* null terminate */
	    s = j*persize;
	    for (k = 0; k < persize; ++k) /* No D_ format in C */
		if ( buf[s+k] == 'D' || buf[s+k] == 'd' ) buf[s+k] = 'E';
	    destination[i++] = atof(&buf[s]);
	    buf[(j+1)*persize] = tmp;     /* recover the char at that place */
	}
    }
}

/*
 * On input, nonz/nzval/rowind/colptr represents lower part of a symmetric
 * matrix. On exit, it represents the full matrix with lower and upper parts.
 */
static void
FormFullA(int_t n, int_t *nonz, double **nzval, int_t **rowind, int_t **colptr)
{
    register int_t i, j, k, col, new_nnz;
    int_t *t_rowind, *t_colptr, *al_rowind, *al_colptr, *a_rowind, *a_colptr;
    int_t *marker;
    double *t_val, *al_val, *a_val;

    al_rowind = *rowind;
    al_colptr = *colptr;
    al_val = *nzval;

    if ( !(marker =(int_t *) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (int_t *) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC t_colptr[]");
    if ( !(t_rowind = (int_t *) SUPERLU_MALLOC( *nonz * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for t_rowind[]");
    if ( !(t_val = (double*) SUPERLU_MALLOC( *nonz * sizeof(double)) ) )
	ABORT("SUPERLU_MALLOC fails for t_val[]");

    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < n; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = al_colptr[j]; i < al_colptr[j+1]; ++i)
	    ++marker[al_rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < n; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose matrix A to T */
    for (j = 0; j < n; ++j)
	for (i = al_colptr[j]; i < al_colptr[j+1]; ++i) {
	    col = al_rowind[i];
	    t_rowind[marker[col]] = j;
	    t_val[marker[col]] = al_val[i];
	    ++marker[col];
	}

    new_nnz = *nonz * 2 - n;
    if ( !(a_colptr = (int_t *) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC a_colptr[]");
    if ( !(a_rowind = (int_t *) SUPERLU_MALLOC( new_nnz * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for a_rowind[]");
    if ( !(a_val = (double*) SUPERLU_MALLOC( new_nnz * sizeof(double)) ) )
	ABORT("SUPERLU_MALLOC fails for a_val[]");
    
    a_colptr[0] = 0;
    k = 0;
    for (j = 0; j < n; ++j) {
      for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	if ( t_rowind[i] != j ) { /* not diagonal */
	  a_rowind[k] = t_rowind[i];
	  a_val[k] = t_val[i];
#ifdef DEBUG
	  if ( fabs(a_val[k]) < 4.047e-300 )
	      printf("%5d: %e\n", k, a_val[k]);
#endif
	  ++k;
	}
      }

      for (i = al_colptr[j]; i < al_colptr[j+1]; ++i) {
	a_rowind[k] = al_rowind[i];
	a_val[k] = al_val[i];
#ifdef DEBUG
	if ( fabs(a_val[k]) < 4.047e-300 )
	    printf("%5d: %e\n", k, a_val[k]);
#endif
	++k;
      }
      
      a_colptr[j+1] = k;
    }

    printf("FormFullA: new_nnz = %d, k = %d\n", new_nnz, k);

    SUPERLU_FREE(al_val);
    SUPERLU_FREE(al_rowind);
    SUPERLU_FREE(al_colptr);
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_val);
    SUPERLU_FREE(t_rowind);
    SUPERLU_FREE(t_colptr);

    *nzval = a_val;
    *rowind = a_rowind;
    *colptr = a_colptr;
    *nonz = new_nnz;
}
