/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include "superlu_ddefs.h"

/* Prototypes */
int_t mc21a_(int_t *, int_t *, int_t *, int_t *, int_t *, int_t *,
	     int_t *, int_t *);


int
zfdperm(int_t n, int_t nnz, int_t *adjncy, int_t *colbeg, int_t *len,
	int_t *perm)
{
/*
 * Purpose
 * =======
 *
 *   ZFDPERM finds a row permutation so that the matrix has no zeros on
 *   its diagonal.
 *
 * Arguments
 * =========
 *
 * n      (input) int
 *        The order of the matrix.
 *
 * nnz    (input) int
 *        The number of nonzeros in the matrix.
 *
 * adjncy (input) int*, of size nnz
 *        The adjacency structure of the matrix, which contains the row
 *        indices of the nonzeros.
 *
 * colbeg (input) int*, of size n
 *        The pointers to the beginning of each column in ADJNCY.
 *
 * len    (input) int*, of size n
 *        The number of nonzeros in each column.
 *
 * perm   (output) int*, of size n
 *        The permutation vector. perm[i] = j means row i in the
 *        original matrix is in row j of the permuted matrix.
 *
 */
 
    int_t i, numnz;
    int_t *iw;

    iw = intMalloc(5*n);
	    
    /* Increment one to get 1-based indexing. */
    for (i = 0; i < nnz; ++i) ++adjncy[i];
    for (i = 0; i < n; ++i)   ++colbeg[i];
	
    /* 
     * NOTE: MC21A assumes:
     *     o matrix A is stored in compressed row format
     *     o perm[i] = j means row i of permuted A is in row j of original A
     * Since our matrix is in compressed column format, we may think
     * of MC21A being applied to the transpose of A, i.e., perm is interpreted
     * as a column permutation vector. Since a symmetric permutation
     * preserves the diagonal entries. Then by the following relation:
     *     P'(A*P')P = P'A
     * we can apply inverse(perm) to rows of A to get zero-free diagonal.
     * But, since perm defined in MC21A happens to be the reverse of our
     * definition of permutation vector, therefore, it is already an inverse
     * for our purpose. We will thus use it directly.
     */
    mc21a_(&n, adjncy, &nnz, colbeg, len, perm, &numnz, iw);
	
    /* Restore to 0-based indexing. */
    for (i = 0; i < nnz; ++i) --adjncy[i];
    for (i = 0; i < n; ++i)   --colbeg[i];
    for (i = 0; i < n; ++i)   --perm[i];
#if ( PRNTlevel>=2 )
    printf(".. Size of maximum matching: numnz %d\n", numnz);
#endif
    free (iw);
}


/*  -- translated by f2c (version 19940927).

    REFERENCE: I.S. Duff, "Algorithm 575 Permutations for a zero-free
               diagonal." ACM Trans. Math. Softw. 7, 387-390, 1981.

   Subroutine */ int_t mc21a_(int_t *n, int_t *icn, int_t *licn, int_t *ip,
			      int_t *lenr, int_t *iperm, int_t *numnz, int_t *iw)
{
/*
 * DESCRIPTION OF PARAMETERS.
 * INPUT VARIABLES  N,ICN,LICN,IP
 * OUTPUT VARIABLES IPERM,NUMNZ
 *
 * N   ORDER OF THE MATRIX.
 * ICN ARRAY CONTAINING THE COLUMN INDICES OF THE NON-ZEROS. THOSE
 *     BELONGING TO A SINGLE ROW MUST BE CONTIGUOUS BUT THE ORDERING
 *     OF COLUMN INDICES WITHIN EACH ROW IS UNIMPORTANT AND WASTED
 *     SPACE BETWEEN ROWS IS PERMITTED.
 * LICN  LENGTH OF ARRAY ICN.
 * IP  IP(I), I=1,2,...N, IS THE POSITION IN ARRAY ICN OF THE FIRST
 *     COLUMN INDEX OF A NON-ZERO IN ROW I.
 * LENR  LENR(I) IS THE NUMBER OF NON-ZEROS IN ROWI, I=1,2,...N.
 * IPERM CONTAINS PERMUTATION TO MAKE DIAGONAL HAVE THE SMALLEST 
 *     NUMBER OF ZEROS ON IT. ELEMENTS (IPERM(I),I) I=1,2,...N ARE
 *     NON-ZERO AT THE END OF THE ALGORITHM UNLESS MATRIX
 *     IS STRUCTURALLY SINGULAR.  IN THIS CASE, (IPERM(I),I) WILL
 *     BE ZERO FOR N-NUMNZ ENTRIES.
 * NUMNZ NUMBER OF NON-ZEROS ON DIAGONAL OF PERMUTED MATRIX.
 * IW  WORK ARRAY  .. SEE LATER COMMENTS.
 *
 */
    /* System generated locals */
    int_t iw_dim1, iw_offset;

    /* Local variables */
    extern /* Subroutine */ int_t mc21b_(int_t *, int_t *, int_t *, int_t *,
					 int_t *, int_t *, int_t *, int_t *,
					 int_t *, int_t *, int_t *);

    /* Parameter adjustments */
    iw_dim1 = *n;
    iw_offset = iw_dim1 + 1;
    iw -= iw_offset;
    --iperm;
    --lenr;
    --ip;
    --icn;

    /* Function Body */
    /*mc21b_(n, &icn[1], licn, &ip[1], &lenr[1], &iperm[1], numnz, &iw[iw_dim1 
	    + 1], &iw[(iw_dim1 << 1) + 1], &iw[iw_dim1 * 3 + 1], &iw[(iw_dim1 
	    << 2) + 1]);*/
    mc21b_(n, &icn[1], licn, &ip[1], &lenr[1], &iperm[1], numnz, &iw[iw_dim1 
	    + 1], &iw[(iw_dim1 * 2) + 1], &iw[iw_dim1 * 3 + 1], &iw[(iw_dim1 
	    * 4) + 1]);
    return 0;
} /* mc21a_ */

/* Subroutine */ int_t mc21b_(int_t *n, int_t *icn, int_t *licn, int_t *ip,
			      int_t *lenr, int_t *iperm, int_t *numnz,
			      int_t *pr, int_t *arp, int_t *cv, int_t *out)
{
    /* System generated locals */
    int_t i__1, i__2, i__3, i__4;

    /* Local variables */
    static int_t jord, i, j, k, j1, ioutk, ii, kk, in1, in2;

/*
     DIVISION OF WORK ARRAY IS NOW DESCRIBED.

     PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.   
           IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.   
           ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE   
           ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE   
           (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.   
     ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I   
         WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.   
     CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I   
         WAS VISITED.   
     OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I   
         WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.   

     INITIALIZATION OF ARRAYS.   
       Parameter adjustments */
    --out;
    --cv;
    --arp;
    --pr;
    --iperm;
    --lenr;
    --ip;
    --icn;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	arp[i] = lenr[i] - 1;
	cv[i] = 0;
/* L10: */
	iperm[i] = 0;
    }
    *numnz = 0;


/*   MAIN LOOP.   
     EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT   
   OR GIVES A ROW WITH NO ASSIGNMENT. */
    i__1 = *n;
    for (jord = 1; jord <= i__1; ++jord) {
	j = jord;
	pr[j] = -1;
	i__2 = jord;
	for (k = 1; k <= i__2; ++k) {
/* LOOK FOR A CHEAP ASSIGNMENT */
	    in1 = arp[j];
	    if (in1 < 0) {
		goto L60;
	    }
	    in2 = ip[j] + lenr[j] - 1;
	    in1 = in2 - in1;
	    i__3 = in2;
	    for (ii = in1; ii <= i__3; ++ii) {
		i = icn[ii];
		if (iperm[i] == 0) {
		    goto L110;
		}
/* L50: */
	    }
/*   NO CHEAP ASSIGNMENT IN ROW. */
	    arp[j] = -1;
/*   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J. */
L60:
	    out[j] = lenr[j] - 1;
/* INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS. */
	    i__3 = jord;
	    for (kk = 1; kk <= i__3; ++kk) {
		in1 = out[j];
		if (in1 < 0) {
		    goto L80;
		}
		in2 = ip[j] + lenr[j] - 1;
		in1 = in2 - in1;
/* FORWARD SCAN. */
		i__4 = in2;
		for (ii = in1; ii <= i__4; ++ii) {
		    i = icn[ii];
		    if (cv[i] == jord) {
			goto L70;
		    }
/*   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS. 
*/
		    j1 = j;
		    j = iperm[i];
		    cv[i] = jord;
		    pr[j] = j1;
		    out[j1] = in2 - ii - 1;
		    goto L100;
L70:
		    ;
		}

/*   BACKTRACKING STEP. */
L80:
		j = pr[j];
		if (j == -1) {
		    goto L130;
		}
/* L90: */
	    }

L100:
	    ;
	}

/*   NEW ASSIGNMENT IS MADE. */
L110:
	iperm[i] = j;
	arp[j] = in2 - ii - 1;
	++(*numnz);
	i__2 = jord;
	for (k = 1; k <= i__2; ++k) {
	    j = pr[j];
	    if (j == -1) {
		goto L130;
	    }
	    ii = ip[j] + lenr[j] - out[j] - 2;
	    i = icn[ii];
	    iperm[i] = j;
/* L120: */
	}

L130:
	;
    }

/*   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE   
   PERMUTATION IPERM. */
    if (*numnz == *n) {
	return 0;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* L140: */
	arp[i] = 0;
    }
    k = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (iperm[i] != 0) {
	    goto L150;
	}
	++k;
	out[k] = i;
	goto L160;
L150:
	j = iperm[i];
	arp[j] = i;
L160:
	;
    }
    k = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (arp[i] != 0) {
	    goto L170;
	}
	++k;
	ioutk = out[k];
	iperm[ioutk] = i;
L170:
	;
    }
    return 0;
} /* mc21b_ */

