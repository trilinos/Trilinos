/*
 * -- SuperLU routine (version 1.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include "dcomplex.h"

int
izmax1_(int *n, doublecomplex *cx, int *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    IZMAX1 finds the index of the element whose real part has maximum   
    absolute value.   

    Based on IZAMAX from Level 1 BLAS.   
    The change is to use the 'genuine' absolute value.   

    Contributed by Nick Higham for use with ZLACON.   

    Arguments   
    =========   

    N       (input) INT   
            The number of elements in the vector CX.   

    CX      (input) COMPLEX*16 array, dimension (N)   
            The vector whose elements will be summed.   

    INCX    (input) INT   
            The spacing between successive values of CX.  INCX >= 1.   

   ===================================================================== 
*/  

    /* System generated locals */
    /* JJH mod 8/10/01 */
    int ret_val/*, i__1, i__2*/;
    /* --JJH */
    double d__1;
    
    /* Local variables */
    double smax;
    int i, ix;

#define CX(I) cx[(I)-1]
/* JJH  8/10/01 Added this macro to eliminate a compiler warning on cplant. */
#define SLU_dabs(x) (((x) > 0.) ? x : (-(x)))
/* --JJH */

    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L30;
    }

/*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    smax = (d__1 = CX(1).r, SLU_abs(d__1));
    ix += *incx;
    /* JJH mod 8/10/01 */
    /*i__1 = *n;*/
    /* --JJH */
    for (i = 2; i <= *n; ++i) {
    /* JJH mod 8/10/01 */
	/*i__2 = ix;*/
    /* --JJH */
	if ((d__1 = CX(ix).r, SLU_dabs(d__1)) <= smax) {
	    goto L10;
	}
	ret_val = i;
    /* JJH mod 8/10/01 */
	/*i__2 = ix;*/
    /* --JJH */
	smax = (d__1 = CX(ix).r, SLU_dabs(d__1));
L10:
	ix += *incx;
/* L20: */
    }
    return ret_val;

/*     CODE FOR INCREMENT EQUAL TO 1 */

L30:
    smax = (d__1 = CX(1).r, SLU_dabs(d__1));
    /* JJH mod 8/10/01 */
    /*i__1 = *n;*/
    /* --JJH */
    for (i = 2; i <= *n; ++i) {
    /* JJH mod 8/10/01 */
	/*i__2 = i;*/
    /* --JJH */
	if ((d__1 = CX(i).r, SLU_dabs(d__1)) <= smax) {
	    goto L40;
	}
	ret_val = i;
    /* JJH mod 8/10/01 */
	/*i__2 = i;*/
    /* --JJH */
	smax = (d__1 = CX(i).r, SLU_dabs(d__1));
L40:
	;
    }
    return ret_val;

/*     End of IZMAX1 */

} /* izmax1_ */

