// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_ConfigDefs.h"
#ifdef HAVE_GALERI_TRIUTILS
#include "Galeri_ReadHB.h"
#include "Trilinos_Util_iohb.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include <vector>
#include <algorithm>

// Copied from Trilinos/packages/triutils/src, this file is here
// only to make it easier for users to create and read matrices.
// Future versions of Galeri might not have this file.
namespace Galeri {

int CSRCSC(int n, int n2, int job, int ipos, double * a, 
           int *ja, int *ia, double *ao, int *jao, int *iao)
{
     
    int next, i, j, k;

/*
 -----------------------------------------------------------------------
 Compressed Sparse Row     to      Compressed Sparse Column 
 (transposition operation)   Not in place. 
 -----------------------------------------------------------------------
 Rectangular version.  n is number of rows of CSR matrix, 
                       n2 (input) is number of columns of CSC matrix. 
 -----------------------------------------------------------------------
 -- not in place -- 
 this subroutine transposes a matrix stored in a, ja, ia format. 
 --------------- 
 on entry: 
 ---------- 
 n	= number of rows of CSR matrix. 
 n2    = number of columns of CSC matrix. 
 job	= integer to indicate whether to fill the values (job.eq.0) of the 

         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1) 

 ipos  = starting position in ao, jao of the transposed matrix. 
        the iao array takes this into account (thus iao(1) is set to ipo
s.)
        Note: this may be useful if one needs to append the data structu
re
         of the transpose to that of A. In this case use for example 
                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
 	  for any other normal usage, enter ipos=1. 
 a	= real array of length nnz (nnz=number of nonzero elements in input 

         matrix) containing the nonzero elements. 
 ja	= integer array of length nnz containing the column positions 
 	  of the corresponding elements in a. 
 ia	= integer of size n+1. ia(k) contains the position in a, ja of 
 	  the beginning of the k-th row. 

 on return: 
 ---------- 
 output arguments: 
 ao	= real array of size nzz containing the "a" part of the transpose 
 jao	= integer array of size nnz containing the column indices. 
 iao	= integer array of size n+1 containing the "ia" index array of 
 	  the transpose. 
 -----------------------------------------------------------------------
*/
 // ----------------- compute lengths of rows of transp(A) ----------------
    for (i = 0; i <= n2; i++) iao[i] = 0;

    for (i = 0; i < n; i++) 
	for (k = ia[i]; k <ia[i+1]; k++) 
	    ++iao[ja[k]+1];

// ---------- compute pointers from lengths ------------------------------
 
    iao[0] = ipos;
    for (i = 0; i < n2; i++)
	iao[i+1] = iao[i] + iao[i+1];

// --------------- now do the actual copying -----------------------------
 
    for (i = 0; i < n; i++) {
	for (k = ia[i]; k <ia[i+1]; k++) {
	    j = ja[k];
	    next = iao[j];
	    if (job == 0) {
		ao[next] = a[k];
	    }
	    jao[next] = i;
	    iao[j] = next + 1;
	}
    }
 
  // -------------------------- reshift iao and leave ----------------------
 
    for (i=n2-1; i >= 0; i--) iao[i+1] = iao[i];
    iao[0] = ipos;
 return(0);
}

int SSRCSR (int job, int value2, int nrow, double *a,
            int *ja, int *ia, int nzmax, 
            double *ao, int *jao, int *iao, int *indu, 
            int *iwk) 
{
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

void  SCSCMV(int isym, int m, int n, double *val, int *indx, int *pntr, 
             double *x, double *y)
{
    int i, j, jbgn, jend;


/*     Performs the matrix-vector operation

                      y = A*x 

       where x and y are vectors and A is a sparse matrix stored
       in (Harwell-Boeing) compress column format. */

/*     -------------------------- 
       First executable statement 
       -------------------------- */

/* .....initialize soln */

    for (i = 0; i < m; i++)
	y[i] = 0.0;

/* .....do a series of SPAXPYs (sparse saxpys) */

    for (j = 0; j <n ; j++) 
      {
	jbgn = pntr[j];
	jend = pntr[j + 1];

	for (i = jbgn; i < jend; i++)
	  {
	    y[indx[i]] += val[i] * x[j];
	    if (indx[i] != j && isym) y[j] += val[i]*x[indx[i]];
	  }
      }

    return;
} /* scscmv */

double SCSCRES(int isym, int m, int n, double *val, int *indx, int *pntr,
               double *x, double *b)
{
    int i, j, ibgn, iend;
    double norm_tmp = 0.0, norm_b = 0.0;
    double scaled_res_norm = 0.0, res_norm = 0.0, *tmp = NULL, max_norm = 0.0;


/*     Computes the residual

                      res = || b - A*x ||

       where x and b are vectors and A is a sparse matrix stored
       in MSR format. */

/*     -------------------------- 
       First executable statement 
       -------------------------- */

    /* Create tmp workspace */
    tmp = (double *) calloc(m,sizeof(double));

/* .....initialize soln */

    for (i = 0; i < m; i++)
	tmp[i] = b[i];

/* .....do a series of SPAXPYs (sparse saxpys) */

    for (j = 0; j < n ; j++) 
      {
	ibgn = pntr[j];
	iend = pntr[j + 1];
	
	for (i = ibgn; i < iend; i++)
	  {
	    tmp[indx[i]] -= val[i] * x[j];
 	    if (indx[i] != j && isym) tmp[j] -= val[i]*x[indx[i]];
	  }
     }
    for (i = 0; i < m; i++)
      {
        if (fabs(tmp[i]) > max_norm) max_norm = fabs(tmp[i]);
	norm_tmp += tmp[i]*tmp[i];
	norm_b += b[i]*b[i];
      }
   
    res_norm = sqrt(norm_tmp);
    printf("\n\nMax norm of residual        = %12.4g\n",max_norm);
    printf(    "Two norm of residual        = %12.4g\n",res_norm);
    if (norm_b > 1.0E-7) 
      {
	   scaled_res_norm = res_norm/sqrt(norm_b);
	   printf(    "Scaled two norm of residual = %12.4g\n",scaled_res_norm);
      }

    free((void *) tmp);

    return(scaled_res_norm);

} /* scscres */

void ReadHB(const char *data_file, const Epetra_Comm  &comm, 
            Epetra_Map *& map, Epetra_CrsMatrix *& A, 
            Epetra_Vector *& x, Epetra_Vector *& b,
            Epetra_Vector *&xexact) 
{
  // In modifying this function for 64 bit Epetra GIDS, we assume that no one
  // would read in a matrix from a file with more than INT_MAX rows.  Thus,
  // 64 bit GIDs are used only when 32 bit GIDs are disabled. 
  // -- Chetan

  FILE *in_file ;
  int numGlobalEquations=0, N_columns=0, n_entries=0, Nrhs=0;
  char Title[73], Key[9], Rhstype[4];
  char Type[4] = "XXX";
  char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
  int Ptrcrd, Indcrd, Valcrd, Rhscrd;
  
  int * bindx = NULL, * pntr = NULL, * indx1 = NULL, * pntr1 = NULL;
  double * val = NULL, * val1 = NULL, * hbx = NULL, * hbxexact, * hbb = NULL;

  hbb = 0; hbxexact = 0; hbb = 0;

  int ok = 1;

  if(comm.MyPID() == 0)  { 
    in_file = fopen( data_file, "r");
    if (in_file == NULL)
    {
      printf("Error: Cannot open file: %s\n",data_file);
      ok = 0;
    }
  }

  // we need to take care that all processors throw an exception here
  comm.Broadcast(&ok, 1, 0);

  if (ok == 0) throw(-1);

  if(comm.MyPID() == 0)  { 

      /* Get information about the array stored in the file specified in the  */
      /* argument list:                                                       */

      printf("Reading matrix info from %s...\n",data_file);
      
      in_file = fopen( data_file, "r");
      if (in_file == NULL)
	{
	  printf("Error: Cannot open file: %s\n",data_file);
	  throw(1);
	}
      
      readHB_header(in_file, Title, Key, Type, &numGlobalEquations, &N_columns, 
		    &n_entries, &Nrhs,
		    Ptrfmt, Indfmt, Valfmt, Rhsfmt,
		    &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
      fclose(in_file);

      if (Nrhs < 0 ) Nrhs = 0;

      printf("***************************************************************\n");
      printf("Matrix in file %s is %d x %d, \n",data_file, numGlobalEquations, N_columns);
      printf("with %d nonzeros with type %3s;\n", n_entries, Type);
      printf("***************************************************************\n");
      printf("Title: %72s\n",Title);
      printf("***************************************************************\n");
      /*Nrhs = 0; */
      printf("%d right-hand-side(s) available.\n",Nrhs);

      if (Type[0] != 'R') perror("Can only handle real valued matrices");
      int isym = 0;
      if (Type[1] == 'S') 
	{
	  printf("Converting symmetric matrix to nonsymmetric storage\n");
	  n_entries = 2*n_entries - N_columns;
	  isym = 1;
	}
      if (Type[2] != 'A') perror("Can only handle assembled matrices");
      if (N_columns != numGlobalEquations) perror("Matrix dimensions must be the same");
      
      /* Read the matrix information, generating the associated storage arrays  */
      printf("Reading the matrix from %s...\n",data_file);

      /* Allocate space.  Note that we add extra storage in case of zero
	 diagonals.  This is necessary for conversion to MSR format. */

      pntr   = (int    *) calloc(N_columns+1,sizeof(int)) ;
      bindx = (int    *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      val   = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      readHB_mat_double(data_file, pntr, bindx, val);

      /* Translate integer arrays to zero base */
      for (int i = 0; i <= numGlobalEquations; i++) pntr[i]--;
      {for (int i = 0; i <= n_entries; i++) bindx[i]--;}

      /* If a rhs is specified in the file, read one, 
	 generating the associate storage */
      if (Nrhs > 0 && Rhstype[2] =='X')
	{
	  printf("Reading right-hand-side vector(s) from %s...\n",data_file);
	  hbb = (double *) calloc(N_columns,sizeof(double));
	  readHB_aux_double(data_file, 'F', hbb);
	  printf("Reading exact solution  vector(s) from %s...\n",data_file);
	  hbxexact = (double *) calloc(N_columns,sizeof(double));
	      readHB_aux_double(data_file, 'X', hbxexact);

	}
      else
	{
	  
	  /* Set Xexact to a random vector */

	  printf("Setting  random exact solution  vector\n");
	  hbxexact = (double *) calloc(N_columns,sizeof(double));
	  
	  for (int i=0;i<numGlobalEquations;i++)	 hbxexact[i] = 
                                              ((double)
                                               rand())/((double) RAND_MAX);
	  
	  /* Compute b to match xexact */
	  
	 hbb = (double   *) calloc(N_columns,sizeof(double)) ;
	  if (hbb == NULL) perror("Error: Not enough space to create rhs");
 

          SCSCMV (isym, N_columns, N_columns, val, bindx, pntr, hbxexact, hbb);
	}

      /* Compute residual using CSC format */

      double res = SCSCRES(isym, numGlobalEquations, numGlobalEquations, val, bindx, pntr, 
		    hbxexact, hbb);
      printf(
	      "The residual using CSC format and exact solution is %12.4g\n",
	      res);

      
      /* Set initial guess to zero */
      
      hbx = (double   *) calloc(numGlobalEquations,sizeof(double)) ;
      
      if (hbx == NULL) 
	perror("Error: Not enough space to create guess");
      
      
      /* Set RHS to a random vector, initial guess to zero */
      {for (int i=0;i<numGlobalEquations;i++) hbx[i] = 0.0;}
      
      
      /* Allocate temporary space */
      
      pntr1 = (int   *) calloc(N_columns+1,sizeof(int)) ;
      indx1 = (int   *) calloc(n_entries+N_columns+1,sizeof(int)) ;
      val1 = (double *) calloc(n_entries+N_columns+1,sizeof(double)) ;
      
      
      /* Convert in the following way:
	 - CSC to CSR 
	 - CSR to MSR
      */
      CSRCSC(numGlobalEquations, numGlobalEquations, 0, 0, val, bindx, pntr, val1, indx1, pntr1);
      
      if (Type[1] == 'S') 
	{
	  int *indu, *iwk;
	  int ierr;
	  indu = new int[N_columns];
	  iwk = new int[N_columns+1];
	  ierr = SSRCSR(3, 1, N_columns, val1, indx1, pntr1, n_entries,
	  		  val1, indx1, pntr1, indu, iwk);
	  delete [] indu;
	  delete [] iwk;
	  if (ierr !=0 ) 
	    {
	    printf(" Error in converting from symmetric form\n  IERR = %d\n",ierr);
	    abort();
	    }
	}
  }
  comm.Broadcast(&numGlobalEquations, 1, 0);
  int nlocal = 0;
  if (comm.MyPID()==0) nlocal = numGlobalEquations;

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  // If 32 bit map constructor is available, this code will call it
  // else 64 bit constructor will be called (and it will be available
  // due to the macro checks above).
  map = new Epetra_Map(numGlobalEquations, nlocal, 0, comm); // Create map with all elements on PE 0
#else
  map = 0;
#endif

  A = new Epetra_CrsMatrix(Copy, *map, 0); // Construct matrix on PE 0
  if (comm.MyPID()==0) {
#if defined(EPETRA_NO_32BIT_GLOBAL_INDICES) && !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
    // Locally copy int indices to long long space.  Not the best way, but at least
    // the code will run.
    for (int i=0; i<numGlobalEquations; i++) {
      if(pntr1[i+1] > pntr1[i]) { // If something to insert
        std::vector<long long> local_LL_ids(pntr1[i+1]-pntr1[i]);
	    std::copy(indx1 + pntr1[i], indx1 + pntr1[i+1], local_LL_ids.begin());
        A->InsertGlobalValues(i, pntr1[i+1]-pntr1[i], val1+pntr1[i], &local_LL_ids[0]);
      }
	}
#elif !defined(EPETRA_NO_32BIT_GLOBAL_INDICES)
    for (int i=0; i<numGlobalEquations; i++)
      A->InsertGlobalValues(i, pntr1[i+1]-pntr1[i], val1+pntr1[i], indx1+pntr1[i]);
#endif
  }
  A->FillComplete();
  
  x = new Epetra_Vector(Copy, *map, hbx);
  b = new Epetra_Vector(Copy, *map, hbb);
  xexact = new Epetra_Vector(Copy, *map, hbxexact);

  Epetra_Vector bcomp(*map);

  A->Multiply(false, *xexact, bcomp);
  double residual;
  bcomp.Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of computed b = " << residual << endl;
  b->Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of given b    = " << residual << endl;
  bcomp.Update(-1.0, *b, 1.0);
  bcomp.Norm2(&residual);
  if (comm.MyPID()==0) cout << "Norm of difference between computed b and given b for xexact = " << residual << endl;
  
  /* Release unneeded space */
  if (comm.MyPID()==0) {
    if (hbb!=0) free((void *) hbb);
    if (hbx!=0) free((void *) hbx);
    if (hbxexact!=0) free((void *) hbxexact);
    free((void *) val);
    free((void *) bindx);
    free((void *) val1);
    free((void *) indx1);
    free((void *) pntr1);
    free((void *) pntr);
  }
  return;
}
#endif //ifdef HAVE_GALERI_TRIUTILS
} // namespace Galeri
