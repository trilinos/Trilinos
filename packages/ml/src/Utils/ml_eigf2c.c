/* ml_eigf2c.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/*#include "f2c.h" */

/* Table of constant values */
#define integer int
#define doublereal double
#define logical int
#define FALSE_ (0)
#define TRUE_ (1)
#include "ml_eigf2c.h"

int ml_pdmout__(USR_COMM *comm, int *lout, int *m, int *n, double *a, 
		  int *lda, int *idigit)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
#ifdef HAVE_ML_PARPACK
    pdmout_(comm,
	    lout, m, n, &a[a_offset], lda, idigit, "Ritz values (Real,Imag) a\
nd direct residuals", (ftnlen)44);

#else
#ifdef HAVE_ML_ARPACK    
    dmout_(lout, m, n, &a[a_offset], lda, idigit, "Ritz values (Real,Imag) a\
nd direct residuals", (ftnlen)44);
   
       
#endif
#endif
    
    return 0;
} /* ml_c_pdmout__ */



int ml_pdneupc__(int *comm,
		int *ivec, char *howmny, int *celect, double *d__, 
		double *v, int *ldv, double *workev,  char *bmat, int *n, 
		char *which, int *nev, double *tol, double *resid, int *ncv, 
		int *iparam, int *ipntr, double *workd, double *workl, 
		int *lworkl, int *ierr, ftnlen howmny_len, ftnlen bmat_len, 
		ftnlen which_len)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;

    /* Builtin functions */
    /*    integer s_wsle(), do_lio(), e_wsle(); */

    /* Local variables */
    static logical rvec;
    static integer i__;
    /* static doublereal sigma, mu; */
    static logical *select;

    /* Parameter adjustments */

    select = (logical *) ML_allocate(sizeof(logical)* (*ncv));
    --workd;
    --resid;
    --workev;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1 * 1;
    v -= v_offset;
    --d__;
    --celect;
    --iparam;
    --ipntr;
    --workl;

    /* Function Body */
    if (*ivec == 0) {
	rvec = FALSE_;
    } else {
	rvec = TRUE_;
    }
    if (*(unsigned char *)howmny == 'S' || iparam[6] == 1) {
	i__1 = *ncv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (celect[i__] == 0) {
		select[i__ - 1] = FALSE_;
	    } else {
		select[i__ - 1] = TRUE_;
	    }
/* L100: */
	}

    } else if (*(unsigned char *)howmny != 'A') {
	printf("Error in haim_pdneupc\n");
	printf("unknown value of howmny %c\n",*howmny);
	exit(1);
    }
    

#ifdef HAVE_ML_PARPACK
    
    pdneupd_(comm,
	     &rvec, "A", select, &d__[1], &d__[*ncv + 1], &v[v_offset], ldv, &
	     sigma, &mu, &workev[1], bmat, n, which, nev, tol, &resid[1], ncv, 
	     &v[v_offset], ldv, &iparam[1], &ipntr[1], &workd[1], &workl[1], 
	     lworkl, ierr, (ftnlen)1, (ftnlen)1, (ftnlen)2);
       
       
#else
#ifdef HAVE_ML_ARPACK    
    
    dneupd_(&rvec, "A", select, &d__[1], &d__[*ncv + 1], &v[v_offset], ldv, &
	    sigma, &mu, &workev[1], bmat, n, which, nev, tol, &resid[1], ncv, 
	    &v[v_offset], ldv, &iparam[1], &ipntr[1], &workd[1], &workl[1], 
	    lworkl, ierr, (ftnlen)1, (ftnlen)1, (ftnlen)2);
   
    printf("\n\t\t Serial arpack iterations\n"); 

#endif
#endif

    ML_free(select);
    return 0;
} /* ml_dneupc__ */


