int Trilinos_Util_coocsr( int nrow, int nnz, double *a, int *ir, int *jc,
                          double *ao, int *jao, int *iao) {
    int i, j, k;
    double x;
    int k0, iad;
/*
 -----------------------------------------------------------------------
  Coordinate     to   Compressed Sparse Row 
 -----------------------------------------------------------------------
 
 converts a matrix that is stored in coordinate format 
  a, ir, jc into a row general sparse ao, jao, iao format. 

 on entry: 
 --------- 
 nrow	= dimension of the matrix 
 nnz	= number of nonzero elements in matrix 
 a, 
 ir, 
 jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz 

         nonzero elements of the matrix with a(k) = actual real value of
 
 	  the elements, ir(k) = its row number and jc(k) = its column 
 	  number. The order of the elements is arbitrary. 

 on return: 
 ----------- 
 ir 	is destroyed 

 ao, jao, iao = matrix in general sparse matrix format with ao 
 	continung the real values, jao containing the column indices, 
 	and iao being the pointer to the beginning of the row, 
 	in arrays ao, jao. 

 Notes: 
 ------ This routine is NOT in place.  See coicsr 

----------------------------------------------------------------------- -
*/
    for (k = 0; k <= nrow; k++) iao[k] = 0;

 // determine row-lengths. 

    for (k = 0; k <nnz ; k++) ++iao[ir[k]];

 // starting position of each row.. 

    k = 0;
    for (j = 0; j <= nrow; j++) {
	k0 = iao[j];
	iao[j] = k;
	k += k0;
    }
 
 // go through the structure  once more. Fill in output matrix. 

    for (k = 0; k <nnz; k++) {
	i = ir[k];
	j = jc[k];
	x = a[k];
	iad = iao[i];
	ao[iad] = x;
	jao[iad] = j;
	iao[i] = iad + 1;
    }
 
  // shift back iao 
    for (j = nrow-1; j >= 0; j--) {
	iao[j + 1] = iao[j];
    }
    iao[0] = 0;
    return 0;
}

