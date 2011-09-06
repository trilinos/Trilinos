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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "az_aztec.h"



void AZ_fact_ilut(int *n, AZ_MATRIX *A_overlapped, double *a, 
		  int *ja, double tol, int fill, int shift, 
		  int *iu,double *nextrow, double *unorm, 
		  int *L_lowest_cols, int *nz_used, int *pattern,
            double rthresh, double athresh)
{

    /* Local variables */

    int lenl, lenu, aptr;
    int i, j, k, m;
    int count, *U_cols, orig_U_len;
    int *L_largest_vals, bindx_prev;
    double  multiplier;
    int *U_largest_vals, jak, col, juj, heaplen;
    double  tmp;
    int row;
    int heap2len;
    double  rownorm;
    int  *bindx;
    double *val;


/*  Incomplete LU factorization of a sparse matrix A with fill and 
 *  drop tolerance.  This variant works on one row of A at a time. 
 *
 *  Parameters
    ==========
 *    Input:
 *       val, bindx - MSR matrix arrays of matrix to be factored.
 *                    On output, matrix factors are stored in these 
 *                    arrays in MSR/Fortran format.
 *
 *             a    - double precision array of length nnzrow used to 
 *                    hold 'current' row during factorization. 
 *             n    - Order of linear system.
 *             nnzrow - Maximum number of nonzero elements in each row of A.
 *             ja   - array of length nnzrow contain column indices
 *                    corresponding to elements of a.
 *             tol  - double precision
 *                    Drop tolerance.  Multipliers L_ik are dropped 
 *                    when | L_ik | < tol*|| A(i,:) ||.  Elements from the
 *                    reduced rows of A are dropped by the same criterion.
 *             fill - Fill level.  Each row of the triangular factor has
 *                    at most (fill + number of nonzeros in A(i,1:i))
 *                    nonzero elements.  (Similarly for the elements in 
 *                    the reduced rows of A.)
 *             shift- The off-diagonals of the input matrix will be shifted
 *                    to the right by 'shift' places to make room for
 *                    the fill-in. NOTE: if tol == 0, shift should be equal
 *                    to fill.
 *         next_row - double precision array of length n + 1.
 *                    Workspace for dense copies of sparse rows of A. 
 *           unorm  - Workspace array of length n. 
 *          pattern - Workspace array of length n.
 *    L_lowest_cols - Workspace array of length n+2.
 *
 *    Output:
 *       val, bindx - matrix factors in MSR/fortran fromat.
 *             iu   - int array of length n+1 
 *                    Array of row pointers into u:  u(iu(i)) is the
 *                    first element of row i in u.
 *	 nz_used    - number of nonzeros required for matrix factorization.
 * 
 *     
 */

    val     = A_overlapped->val;
    bindx   = A_overlapped->bindx;

   /* Scale diagonal if necessary */
   if (rthresh==0.0) rthresh = 1.0; /* Default means no rthresh */
   if (rthresh != 1.0 || athresh != 0.0)
     for (i = 0 ; i < *n ; i++ ) {
        if (val[i]>=0)
           val[i] = val[i] * rthresh + athresh;
        else
           val[i] = val[i] * rthresh - athresh;
     }

    /* Parameter adjustments */

    --unorm;
    --pattern;
    --iu;
    --ja;
    --a;
    --bindx;
    --val;


    L_largest_vals = (int *) AZ_allocate( (fill+1)*sizeof(int));
    U_largest_vals = (int *) AZ_allocate( (fill+1) * sizeof(int));
    if (U_largest_vals == NULL) 
       AZ_perror("Not enough space inside ilut\n");

    /* stick in some dummy values so that the case 'fill=0' works. */
    /* Specifically, this makes the following statement always     */
    /* false:                                                      */
    /*   if (fabs(multiplier) < fabs(nextrow[L_largest_vals[0]-1]))*/
          
    L_largest_vals[0] = *n+1;
    U_largest_vals[0] = *n+1;
    nextrow[*n] = DBL_MAX;


    /* Function Body */
    for (i = 1; i <= *n; ++i) pattern[i] = 0;

    /*  Initialize pointers to next open location */
    /*  in triangular factor and reduced rows of A. */

    aptr = *n + 2;

    /*  shift storage  from C to Fortran */

    m = bindx[*n + 1] - 1;
    for (i = 1; i <= m+1; ++i) ++bindx[i];

    /*  shift off-diagonals in original matrix to accomodate */
    /*  fill-in that will occur during the factorization.    */

    for (i = m + 1; i >= *n + 2; --i) {
	val[i + shift] = val[i];
	bindx[i + shift] = bindx[i];
    }
    for (i = 1; i <= *n+1 ; ++i) {
	bindx[i] += shift;
    }
    bindx_prev = bindx[1];
    bindx[1] = *n + 2;

    for (i = 1; i <= *n ; ++i) {
	lenl = 0;
	lenu = 1;
        U_cols = &(L_lowest_cols[i]);

        /* put row i into 'nextrow' */

	nextrow[i - 1] = val[i];
	U_cols[1]      = i;
	pattern[i]     = 1;
	rownorm        = fabs(nextrow[i - 1]);

	for (k = bindx_prev; k <= bindx[i + 1] - 1; ++k) {
	    jak = bindx[k];
	    if (jak < i) {
		AZ_put_in_heap(L_lowest_cols, &jak, &lenl);
	    } else if (jak > i) {
		++lenu;
		U_cols[lenu] = jak;
	    }
	    nextrow[jak - 1] = val[k];
	    pattern[jak] = 1;
	    rownorm += fabs(nextrow[jak - 1]);
	}
	orig_U_len = lenu;
	rownorm = tol * rownorm / (double ) (lenu + lenl);

        /* get rid of small values of 'nextrow' */

	for (k = 1; k <= lenl; ++k) {
	    if (fabs(nextrow[L_lowest_cols[k - 1] - 1]) < rownorm) {
		nextrow[L_lowest_cols[k - 1] - 1] = 0.0;
	    }
	}
	for (k = 2; k <= lenu ; ++k) {
	    if (fabs(nextrow[U_cols[k] - 1]) < rownorm) {
		nextrow[U_cols[k] - 1] = 0.0;
	    }
	}
	rownorm /= (double ) (lenu + lenl);

        /*  Main reduction loop:  calculate row i of L and U using only  */
        /*  previous rows corresponding to nonzero elements in the first */
        /*  part of the row. */

	heaplen = lenl;
	heap2len = 0;
L30:
	if (heaplen == 0) goto L60;

        /* row is smallest remaining column in row i */
	row = L_lowest_cols[0];
	multiplier = nextrow[row - 1] * val[row];
	nextrow[row - 1] = multiplier;
	if (fabs(multiplier) * unorm[row] < rownorm) {
	    pattern[row] = 0;
	    AZ_rm_heap_root(L_lowest_cols, &heaplen);
	    goto L30;
	}
	if (pattern[row] != 1) {
	    if (heap2len < fill) {
		AZ_put_in_dbl_heap(&row, nextrow, L_largest_vals, &heap2len);
	    } else if (fabs(multiplier) < fabs(nextrow[L_largest_vals[0]-1])) {
		pattern[row] = 0;
		AZ_rm_heap_root(L_lowest_cols, &heaplen);
		goto L30;
	    } else {
		pattern[L_largest_vals[0]] = 0;
		AZ_rm_dbl_heap_root(L_largest_vals, nextrow, &heap2len);
		AZ_put_in_dbl_heap(&row, nextrow, L_largest_vals, &heap2len);
	    }
	}

        /* Reduce current row */

	for (j = iu[row]; j <= bindx[row + 1] - 1; ++j) {
	    tmp = multiplier * val[j];
	    juj = bindx[j];
	    if (pattern[juj] != 0) {
		nextrow[juj - 1] -= tmp;
	    } else if (fabs(tmp) >= rownorm) {
		nextrow[juj - 1] = -tmp;
		pattern[juj] = 2;
		if (juj >= i) {
		    ++lenu;
		    U_cols[lenu] = juj;
		} else {
		    AZ_put_in_heap(L_lowest_cols, &juj, &heaplen);
		}
	    }
	}
	AZ_rm_heap_root(L_lowest_cols, &heaplen);
	goto L30;
L60:
	if (aptr + lenl + heap2len + orig_U_len + fill - 2 >= bindx[i + 1]) {
            AZ_printf_err("ERROR: not enough memory for ILUT. Decrease\n");
            AZ_printf_err("       fill-in or increase drop tolerance.\n");
            exit(-1);
        }
	k = 0;
	for (j = bindx_prev; j <= bindx[i + 1] - 1; ++j) {
	    if (bindx[j] < i) {
		++k;
		val[aptr] = nextrow[bindx[j] - 1];
		bindx[aptr] = bindx[j];
		pattern[bindx[j]] = 0;
		++aptr;
	    }
	}
	for (j = 1; j <= heap2len; ++j) {
	    val[aptr] = nextrow[L_largest_vals[j - 1] - 1];
	    bindx[aptr] = L_largest_vals[j - 1];
	    pattern[L_largest_vals[j - 1]] = 0;
	    ++aptr;
	}
	iu[i] = aptr;

        /*  Pick out the diagonal element, store its reciprocal. */
	if (nextrow[i - 1] == 0.) {
            AZ_printf_err("ilut: zero pivot encountered!\n");
	    nextrow[i - 1] = rownorm;
	}
	pattern[i] = 0;
	val[i] = 1. / nextrow[i - 1];
	unorm[i] = fabs(nextrow[i - 1]);
	count = 0;
	for (j = 2; j <= lenu; ++j) {
	    col = U_cols[j];
	    if (pattern[col] != 1) {
		if (count < fill) {
		    AZ_put_in_dbl_heap(&col, nextrow, U_largest_vals, &count);
		} else if (fabs(nextrow[col-1]) > 
                           fabs(nextrow[U_largest_vals[0] - 1])) {
		    AZ_rm_dbl_heap_root(U_largest_vals, nextrow, &count);
		    AZ_put_in_dbl_heap(&col, nextrow, U_largest_vals, &count);
		}
	    } else {
		val[aptr] = nextrow[col - 1];
		bindx[aptr] = col;
		unorm[i] += fabs(val[aptr]);
		pattern[col] = 0;
		++aptr;
	    }
	    pattern[col] = 0;
	}
	for (j = 1; j <= count; ++j) {
	    val[aptr] = nextrow[U_largest_vals[j - 1] - 1];
	    bindx[aptr] = U_largest_vals[j - 1];
	    unorm[i] += fabs(val[aptr]);
	    ++aptr;
	}
	bindx_prev = bindx[i + 1];
	bindx[i + 1] = aptr;
	unorm[i] /= (double ) (orig_U_len + count + 1);
    }

    iu[*n + 1] = aptr;

    AZ_free(U_largest_vals);
    AZ_free(L_largest_vals);

    if (*n == 0) *nz_used = 0;
    else *nz_used = aptr;

} /* az_ilut */


/***********************************************************/
/***********************************************************/
/***********************************************************/
void AZ_rm_heap_root(int heap[], int *length)
{
/*
 * Remove the first element from the heap.
 * 
 * Parameters
 * ==========
 *
 * heap :   On input, a heap (see Knuth Vol 2) where
 *          the root is stored in the first element
 *          and the children of heap[k] are heap[2*k]
 *          and heap[2*k+1].
 *          On output, the root element is removed and the
 *          heap is reorganized so that it remains a heap
 *
 * length:  Length of the heap array.
 *
 */

 int parent,flag, left,right,prev,child;

 heap--;   /* fortran c conversion */
 parent = 1;
 flag   = 1;

 while(flag) {
    left   = 2*parent;
    right  = left+1;
    prev   = parent;
    if (right <= *length) {
       if ( heap[left] < heap[right]) parent = left;
       else parent = right;
       heap[prev] = heap[parent];
    }
    else if (left == *length) {
       parent = left;
       heap[prev] = heap[parent];
    }  
    else flag = 0;
 }
   
 child = parent;
 if (parent == 1) { (*length)--; return;}
 parent = child/2;
 while( heap[parent] > heap[*length] ) {
      heap[child] = heap[parent];
      child = parent;
      parent = child/2;
 }
 heap[child] = heap[*length]; 
   
 (*length)--;
}
   
/***********************************************************/
/***********************************************************/
/***********************************************************/

void AZ_put_in_heap(int heap[], int *val,int *length)
{
/*
 * Put 'val' into the heap.
 * 
 * Parameters
 * ==========
 *
 * heap :   On input, a heap (see Knuth Vol 2) where
 *          the root is stored in the first element
 *          and the children of heap[k] are heap[2*k]
 *          and heap[2*k+1].
 *          On output, '*val' is placed in the heap.
 *          Note: the heap is reorganized so that it 
 *          remains a heap
 *
 * val :    Element to be added to the heap
 *
 * length:  Length of the heap array.
 *
 */
   int next,prev;

   heap--;    /* fortran/C conversion */

   prev  = *length+1;
   next  = prev/2;

   while( (next !=0) && (*val < heap[next]) ) {
       heap[prev] = heap[next];
       prev = next;
       next = prev/2;
   }
   heap[prev] = *val;
   (*length)++;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void AZ_rm_dbl_heap_root(int heap[], double vals[], int *length)
{
/*
 * Remove the first element from the heap.
 * NOTE: this routine differs from AZ_rm_heap_root() in that
 * the data is actually contained in vals[]. heap[] is
 * an index into vals[] which contains the real data.
 * 
 * Parameters
 * ==========
 *
 * heap :   On input, a POINTER INTO VALS representing a
 *          heap (see Knuth Vol 2) where the root is 
 *          vals[heap[0]] and the children of vals[heap[k]] 
 *          are vals[heap[2*k]] and vals[heap[2*k+1]].
 *          On output, the root element is removed and the
 *          heap is reorganized so that it remains a heap
 *
 * length:  Length of the heap array.
 *
 */
   int parent,flag, left,right,prev,child;

   vals--; heap--;   /* fortan/C conversion */
   parent = 1;
   flag   = 1;

   while(flag) {
      left   = 2*parent;
      right  = left+1;
      prev   = parent;
      if (right <= *length) {
         if ( fabs(vals[heap[left]]) < fabs(vals[heap[right]])) parent = left;
         else parent = right;
         heap[prev] = heap[parent];
      }
      else if (left == *length) {
         parent = left;
         heap[prev] = heap[parent];
      }  
      else flag = 0;
   }
   
   child = parent;
   if (parent == 1) { (*length)--;  return; }
   parent = child/2;
   while( fabs(vals[heap[parent]]) > fabs(vals[heap[*length]]) ) {
      heap[child] = heap[parent];
      child = parent;
      parent = child/2;
   }
   heap[child] = heap[*length]; 
   (*length)--;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/


void AZ_put_in_dbl_heap(int *row, double vals[], int heap[], 
	int *length)
{
/*
 * Put 'vals[row]' into the heap.
 * NOTE: this routine differs from AZ_put_in_heap() in that
 * the data is actually contained in vals[]. heap[] is
 * an index into vals[] which contains the real data.
 * 
 * Parameters
 * ==========
 *
 * heap :   On input, a POINTER INTO VALS representing a
 *          heap (see Knuth Vol 2) where the root is 
 *          vals[heap[0]] and the children of vals[heap[k]] 
 *          are vals[heap[2*k]] and vals[heap[2*k+1]].
 *          On output, vals[row] is incorporated into the
 *          heap.
 *
 * length:  Length of the heap array.
 *
 */
   int next,prev;

   heap--; vals--;   /* fortran/C conversion */

   prev  = *length+1;
   next  = prev/2;

   while( (next !=0) && (fabs(vals[*row]) < fabs(vals[heap[next]])) ) {
       heap[prev] = heap[next];
       prev = next;
       next = prev/2;
   }
   heap[prev] = *row;
   (*length)++;
}
