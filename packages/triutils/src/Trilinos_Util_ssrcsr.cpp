// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Trilinos_Util.h"

int Trilinos_Util_ssrcsr ( int job, int value2, int nrow, double *a,
                          int *ja, int *ia, int nzmax, 
                          double *ao, int *jao, int *iao, int *indu, 
                          int *iwk) {

    int ipos, i, j, k, klast, kosav, ko, kfirst;
    double tmp;
    int nnz;
/*
 -----------------------------------------------------------------------
 
     Symmetric Sparse Row to Compressed Sparse Row format 
 -----------------------------------------------------------------------
 
     This subroutine converts a given matrix in SSR format to regular 
     CSR format by computing Ao = A + A' - diag(A), where A' is A 
     transpose. 

     Typically this routine is used to expand the SSR matrix of 
     Harwell Boeing matrices, or to obtain a symmetrized graph of 
     unsymmetric matrices. 

     This routine is inplace, i.e., (Ao,jao,iao) may be same as 
     (a,ja,ia). 

     It is possible to input an arbitrary CSR matrix to this routine, 
     since there is no syntactical difference between CSR and SSR 
     format. It also removes duplicate entries and perform a partial 
     ordering. The output matrix has an order of lower half, main 
     diagonal and upper half after the partial ordering. 
 -----------------------------------------------------------------------
 
 on entry: 
 --------- 

 job   = options 
         0 -- duplicate entries are not removed. If the input matrix is 

             SSR (not an arbitary CSR) matrix, no duplicate entry should

              arise from this routine. 
         1 -- eliminate duplicate entries, zero entries. 
         2 -- eliminate duplicate entries and perform partial ordering. 

         3 -- eliminate duplicate entries, sort the entries in the 
              increasing order of clumn indices. 

 value2= will the values of A be copied? 
         0 -- only expand the graph (a, ao are not touched) 
         1 -- expand the matrix with the values. 

 nrow  = column dimension of inout matrix 
 a, 
 ia, 
 ja    = matrix in compressed sparse row format. 

 nzmax = size of arrays ao and jao. SSRCSR will abort if the storage 
          provided in ao, jao is not sufficient to store A. See ierr. 

 on return: 
 ---------- 
 ao, jao, iao 
       = output matrix in compressed sparse row format. The resulting 
         matrix is symmetric and is equal to A+A'-D. ao, jao, iao, 
         can be the same as a, ja, ia in the calling sequence. 

 indu  = integer array of length nrow. INDU will contain pointers 
         to the beginning of upper traigular part if job > 1. 
         Otherwise it is also used as a work array (size nrow). 

 iwk   = integer work space (size nrow+1). 

 return flag = integer. Serving as error message. If the length of the arrays 

         ao, jao exceeds nzmax, returns the minimum value 
         needed for nzmax. otherwise 0 is normal return. 

 -----------------------------------------------------------------------
*/ 
    for (i = 0; i <nrow; i++) {
	indu[i] = 0;
	iwk[i] = 0;
    }
    iwk[nrow] = 0;

     // .. compute number of elements in each row of (A'-D) 
     // put result in iwk(i+1)  for row i. 

    for (i = 0; i <nrow ; i++) {
	for (k = ia[i]; k <ia[i+1]; k++) {
	    j = ja[k];
	    if (j != i) {
		++iwk[j+1];
	    }
	}
    }

    // .. find addresses of first elements of ouput matrix. result in iwk 


    iwk[0] = 0;
    for (i = 0; i <nrow; i++) {
	indu[i] = iwk[i] + ia[i+1] - ia[i];
	iwk[i+1] += indu[i];
	--indu[i];
    }
 // .....Have we been given enough storage in ao, jao ? 
    nnz = iwk[nrow];
    if (nnz > nzmax) return(nnz);

     // .. copy the existing matrix (backwards). 

    kosav = iwk[nrow];
    for (i = nrow-1; i >= 0; i--) {
	klast = ia[i+1] - 1;
	kfirst = ia[i];
	iao[i+1] = kosav;
	kosav = iwk[i];
	ko = iwk[i] - kfirst;
	iwk[i] = ko + klast + 1;
	for (k = klast; k >= kfirst; k--) {
	    if (value2 != 0) {
		ao[k+ko] = a[k];
	    }
	    jao[k+ko] = ja[k];
	}
    }
    iao[0] = 0;

     // now copy (A'-D). Go through the structure of ao, jao, iao 
     // that has already been copied. iwk(i) is the address 
     // of the next free location in row i for ao, jao. 

    for (i = 0; i <nrow; i++) {
	for (k = iao[i]; k <= indu[i]; k++) {
	    j = jao[k];
	    if (j != i) {
		ipos = iwk[j];
		if (value2 != 0) {
		    ao[ipos] = ao[k];
		}
		jao[ipos] = i;
		iwk[j] = ipos + 1;
	    }
	}
    }
    if (job <= 0) {
	return(0);
    }

     // .. eliminate duplicate entries -- 
     // array INDU is used as marker for existing indices, it is also the 
     // location of the entry. 
     // IWK is used to stored the old IAO array. 
     // matrix is copied to squeeze out the space taken by the duplicated 
     // entries. 

    for (i = 0; i < nrow; i++) {
	indu[i] = 0;
	iwk[i] = iao[i];
    }

    iwk[nrow] = iao[nrow];
    k = 0;
    for (i = 0; i < nrow; i++) {
	iao[i] = k;
	ipos = iwk[i];
	klast = iwk[i+1];
	while (ipos < klast) {
	    j = jao[ipos];
	    if (indu[j] == 0) {
     // .. new entry .. 
		if (value2 != 0) {
		    if (ao[ipos] != 0.) {
			indu[j] = k;
			jao[k] = jao[ipos];
			ao[k] = ao[ipos];
			++k;
		    }
		} else {
		    indu[j] = k;
		    jao[k] = jao[ipos];
		    ++k;
		}
	    } else if (value2 != 0) {
     // .. duplicate entry .. 
		ao[indu[j]] += ao[ipos];
	    }
	    ++ipos;
	}
     // .. remove marks before working on the next row .. 
	for (ipos = iao[i]; ipos < k; ipos++) indu[jao[ipos]] = 0; 
    }
    iao[nrow] = k;
    if (job <= 1) {
	return 0;
    }

     // .. partial ordering .. 
     // split the matrix into strict upper/lower triangular 
     // parts, INDU points to the the beginning of the strict upper part. 


    for (i = 0; i < nrow; i++) {
	klast = iao[i+1] - 1;
	kfirst = iao[i];
	while (klast > kfirst) {
	    if (jao[klast] < i && jao[kfirst] >= i) {
     // .. swap klast with kfirst .. 
		j = jao[klast];
		jao[klast] = jao[kfirst];
		jao[kfirst] = j;
		if (value2 != 0) {
		    tmp = ao[klast];
		    ao[klast] = ao[kfirst];
		    ao[kfirst] = tmp;
		}
	    }
	    if (jao[klast] >= i) {
		--klast;
	    }
	    if (jao[kfirst] < i) {
		++kfirst;
	    }
	}

	if (jao[klast] < i) {
	    indu[i] = klast + 1;
	} else {
	    indu[i] = klast;
	}
    }
    if (job <= 2) {
	return 0;
    }

     // .. order the entries according to column indices 
     // bubble-sort is used 

    for (i = 0; i <nrow; i++) {
	for (ipos = iao[i]; ipos <indu[i]; ipos++) {
	    for (j = indu[i]-1; j >ipos; j--) {
		k = j - 1;
		if (jao[k] > jao[j]) {
		    ko = jao[k];
		    jao[k] = jao[j];
		    jao[j] = ko;
		    if (value2 != 0) {
			tmp = ao[k];
			ao[k] = ao[j];
			ao[j] = tmp;
		    }
		}
	    }
	}
	for (ipos = indu[i]; ipos <iao[i+1]; ipos++) {
	    for (j = iao[i+1]-1; j >ipos; j--) {
		k = j - 1;
		if (jao[k] > jao[j]) {
		    ko = jao[k];
		    jao[k] = jao[j];
		    jao[j] = ko;
		    if (value2 != 0) {
			tmp = ao[k];
			ao[k] = ao[j];
			ao[j] = tmp;
		    }
		}
	    }
	}
    }
    return 0;
}

