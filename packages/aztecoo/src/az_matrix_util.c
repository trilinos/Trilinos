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

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "az_aztec.h"
#include "az_blas_wrappers.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gmax_matrix_norm(double val[], int indx[],int bindx[], int rpntr[],
                           int cpntr[], int bpntr[], int proc_config[],
                           int data_org[])

/*******************************************************************************

  This routine returns maximum matrix norm for VBR matrix A: this is a parallel
  version for a distributed matrix.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     double, maximum matrix norm.
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see file
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  register int indx_ptr, irow, iblk_row, jcol, jblk_col, icol_blk, iblk = 0;
  int           num_blk_rows, num_col_blks, num_blk_cols;
  double        row_sum = 0.0, row_max = 0.0;
  int           k;
  int           j_last, bindx_row;

  /**************************** execution begins ******************************/

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {

    for (irow = 0; irow < data_org[AZ_N_internal] + data_org[AZ_N_border];
         irow++) {

      /* compute diagonal contribution */

      row_sum = fabs(val[irow]);

      /* nonzero off diagonal contibution */

      j_last  = bindx[irow+1] - bindx[irow];
      bindx_row = bindx[irow];

      for (jcol = 0; jcol < j_last; jcol++) {
        k        = bindx_row + jcol;
        row_sum += fabs(val[k]);
      }
      row_max = AZ_MAX(row_sum, row_max);
    }

    row_max = AZ_gmax_double(row_max, proc_config);

    return row_max;
  }

  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {

    /* loop over the block rows */

    for (iblk_row = 0;
         iblk_row < data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
         iblk_row++) {

      /* find out how many rows are in this block row */

      num_blk_rows = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* find out how many block columns are in this block row */

      num_col_blks = bpntr[iblk_row+1] - bpntr[iblk_row];

      /* loop over all the rows in this block row */

      for (irow = 0; irow < num_blk_rows; irow++) {

        /* loop over all the column blocks in this block row */

        for (jblk_col = 0; jblk_col < num_col_blks; jblk_col++) {

          /* find out which column block this is */

          icol_blk = bindx[iblk];
          indx_ptr = indx[iblk++];

          /* find out how many columns are in this block */

          num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

          /* loop over all the columns in this block */

          for (jcol = 0; jcol < num_blk_cols; jcol++)
            row_sum += fabs(val[indx_ptr + jcol*num_blk_rows + irow]);
        }

        iblk   -= num_col_blks;
        row_max = AZ_MAX(row_sum, row_max);
        row_sum = 0.0;
      }

      iblk += num_col_blks;
    }

    row_max = AZ_gmax_double(row_max, proc_config);

    return row_max;
  }

  else {
    (void) AZ_printf_err( "ERROR: invalid matrix type %d\n",
                   data_org[AZ_matrix_type]);
    exit(1);
  }
  return(0.0);

} /* AZ_gmax_matrix_norm */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gvector_norm(int n, int p, double *x, int proc_config[])

/*******************************************************************************

  Function which returns the lp norm of the vector x, i.e., if p = 2, the
  standard Euclidean norm is returned. NOTE, to get the l-infinity norm,
  set p = -1.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     double, requested norm value.
  ============

  Parameter list:
  ===============

  n:               Order of vector x.

  p:               Order of the norm to perform, i.e., ||x||p.

  x:               Vector of length n (on this processor).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  register int    i, j;
  register double sum, power;
  int             index, one = 1, *n_ptr;
  double          norm;

  /**************************** execution begins ******************************/

  /* error checking */

  if (p <= 0 && p != -1) return -1.0;

  /* initialize */

  n_ptr = &n;

  switch (p) {

  case -1:                      /* infinity norm */

    index = IDAMAX_F77(n_ptr, x, &one) - 1;
    if (index < 0) return -1.0;
    norm  = AZ_gmax_double(fabs(x[index]), proc_config);
    break;

  case 1:                       /* sum norm */
    sum  = DASUM_F77(n_ptr, x, &one);
    norm = AZ_gsum_double(sum, proc_config);
    break;

  case 2:                       /* Euclidean norm */
    norm = sqrt(AZ_gdot(n, x, x, proc_config));
    break;

  default:                      /* any other p-norm */

    sum = 0.0;
    for (i = 0; i < n; i++) {
      power = *x;
      for (j = 0; j < p; j++)
        power *= *x;
      sum += fabs(power);
      x++;
    }

    norm = pow(AZ_gsum_double(sum, proc_config), 1.0 / (double) p);
  }

  return norm;

} /* vector norm */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_vbr_matrix(int matrix_flag, int Proc, int itotal_nodes,
                      int ext_nodes, double val[], int indx[], int bindx[],
                      int rpntr[], int bpntr[])

/*******************************************************************************

  Prints out the VBR matrix and pointers.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  matrix_flag:     = 0 no matrix output, = 1 output matrix.

  Proc:            Current processor number.

  itotal_nodes:    Number of internal + border nodes on this processor.

  ext_nodes:       Number of external nodes.

  val:             Array containing the nonzero entries of the matrix (see file
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  int iblk, i, iblk_row, m1, n1, ipoint;
  int ival, jblk, j, ib1, ib2, jpoint;

  /**************************** execution begins ******************************/

  /* print out the VBR indexing information for the matrix */

  (void) AZ_printf_out("\n----- Proc: %d indx -----\n\n", Proc);

  for (iblk = 0; iblk < itotal_nodes; iblk++) {
    for (i = *(bpntr+iblk); i < *(bpntr+iblk+1); i++)
      (void) AZ_printf_out("%d ", *(indx+i));

    if (iblk == itotal_nodes - 1){
      (void) AZ_printf_out("%d\n", *(indx+i));
    }
    else
      (void) AZ_printf_out("\n");
  }

  (void) AZ_printf_out("\n----- Proc: %d bindx -----\n\n", Proc);

  for (iblk = 0; iblk < itotal_nodes; iblk++) {
    for (i = *(bpntr+iblk); i < *(bpntr+iblk+1); i++)
      (void) AZ_printf_out("%d ", *(bindx+i));
    (void) AZ_printf_out("\n");
  }

  (void) AZ_printf_out("\n----- Proc: %d rpntr -----\n\n", Proc);

  for (i = 0; i < itotal_nodes + ext_nodes + 1; i++)
    (void) AZ_printf_out("%d ", *(rpntr+i));
  (void) AZ_printf_out("\n");

  (void) AZ_printf_out("\n----- Proc: %d bpntr -----\n\n", Proc);

  for (i = 0; i < itotal_nodes + 1; i++)
    (void) AZ_printf_out("%d ", *(bpntr+i));
  (void) AZ_printf_out("\n");

  /* dump of matrix in a block output format */

  if (matrix_flag) {

    /* loop over block rows */

    for (iblk_row = 0; iblk_row < itotal_nodes; iblk_row++) {

      /* number of rows in the current row block */

      m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* starting index of current row block */

      ival = indx[bpntr[iblk_row]];

      /* loop over all the blocks in the current block-row */

      for (j = bpntr[iblk_row]; j < bpntr[iblk_row+1]; j++) {
        jblk = bindx[j];

        /* the starting point column index of the current block */

        ib1 = rpntr[jblk];

        /* ending point column index of the current block */

        ib2 = rpntr[jblk+1];

        /* number of columns in the current block */

        n1 = ib2 - ib1;

        (void) AZ_printf_out("\nProc: %d Block Row: %d Block Column: %d "
                      "Row Pointer: %d Column Pointer: %d\n", Proc, iblk_row,
                      jblk, rpntr[iblk_row], rpntr[jblk]);

        (void) AZ_printf_out("---------------------------------------"
                      "---------------------------------------\n");

        for (ipoint = 0; ipoint < m1; ipoint++) {
          for (jpoint = 0; jpoint < n1; jpoint++)
            (void) AZ_printf_out("val[%d]: %e ", ival+jpoint*m1+ipoint,
                          val[ival+jpoint*m1+ipoint]);
          (void) AZ_printf_out("\n");
        }

        ival += m1*n1;
      }
    }
  }

} /* print_vbr_matrix */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_dtrans(int *m, int *n, double *A)

/*******************************************************************************

  Perform an in-place transpose of a general m x n matrix stored in "A".

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m, n:            Number of rows, columns in "A", respectively.

  A:               Matrix to be transposed.

*******************************************************************************/

{

  /* local variables */

  register int  i, j, w_index = 0, A_index;
  int           itemp;
  double       *work;

  /**************************** execution begins ******************************/

  work = (double *) AZ_allocate(*m * *n * sizeof(double));

  /* do the transpose */

  for (i = 0; i < *n; i++)
    for (j = 0; j < *m; j++) {
      A_index = i + j * *n;
      *(work + w_index++) = *(A + A_index);
    }

  /* copy from "work" back to "A" */

  for (i = 0; i < *m * *n; i++)
    *(A + i) = *(work + i);


  AZ_free((void *) work);

  /* exchange m and n */

  itemp = *m;
  *m    = *n;
  *n    = itemp;

} /* dtrans */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_get_sym_indx(int iblk, int jblk, int *indx, int *bindx, int *bpntr)

/*******************************************************************************

  Given a block-row and block-column index, return the index to the symmetric
  starting point of the matrix.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     int, starting index into val of the symmetric portion of the
  ============          matrix.

  Parameter list:
  ===============

  iblk:            Current block row index.

  jblk:            Current block column index.

  indx,
  bindx,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

*******************************************************************************/

{

  register int i, icount = 0;
  int          itemp, pre_num_nz_blks = 0, num_nz_blks;
  int         *bindx_ptr;

  itemp     = bpntr[jblk];
  bindx_ptr = bindx + itemp;

  /*
   * Count up how many nonzero block precede the current iblk' in the current
   * block row.
   */

  num_nz_blks = bpntr[jblk+1] - itemp;

  for (i = 0; i < num_nz_blks; i++) {
    if (*(bindx_ptr + icount) == iblk) {
      pre_num_nz_blks = icount;
      break;
    }
    else
      icount++;
  }

  return indx[itemp + pre_num_nz_blks];

} /* AZ_get_sym_indx */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_out(int update_index[], int extern_index[], int update[], 
	int external[], double val[], int indx[], int bindx[], int rpntr[], 
	int cpntr[], int bpntr[], int proc_config[], int choice, int matrix, 
	int N_update, int N_external, int of)
{
/*******************************************************************************

  Print out the matrix in 1 of 3 formats described below.
  starting point of the matrix.

  Author:          Ray Tuminaro, SNL, 9222
  =======

  Return code:     none.
  ============ 

  Parameter list:
  ===============

  update_index,    On input, ordering of update and external locally on this
  extern_index:    processor. For example  'update_index[i]' gives the index
                   location of the block which has the global index 'update[i]'.
                   (AZ_global_mat only).
 
  update:          On input, blks updated on this node (AZ_global_mat only).

  external:        On input, list of external blocks (AZ_global_mat only).

  val,bindx        On input, matrix (MSR or VBR) arrays holding matrix values.
  indx, bnptr,     When using either AZ_input_form or AZ_explicit, these can
  rnptr, cnptr:    be either pre or post-AZ_transform() values depending on what
                   the user wants to see. When using AZ_global_form, these must
                   be post-AZ_transform() values.

  proc_config:     On input, processor information corresponding to:
                      proc_config[AZ_node] = name of this processor
                      proc_config[AZ_N_procs] = # of processors used
 
  choice:          On input, 'choice' determines the output to be printed
		   as described above.
 
  matrix:          On input, type of matrix (AZ_MSR_MATRIX or AZ_VBR_MATRIX).

  N_update:        On input, number of points/blks to be updated on this node.

  N_external:      On input, number of external points/blks needed by this node.

  of:              On input, an offset used with AZ_global_matrix and 
		   AZ_explicit. In particular, a(of,of) is printed for 
                   the matrix element stored as a(0,0).

*******************************************************************************/

   int type, neighbor, cflag;
   int ii,i = 1,j,tj;
   int iblk, jblk, m1, n1, ival, new_iblk, new_jblk;
   MPI_AZRequest request;  /* Message handle */

 
   type            = proc_config[AZ_MPI_Tag];
   proc_config[AZ_MPI_Tag] =(type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +AZ_MSG_TYPE;

   /* Synchronize things so that processor 0 prints first and then */
   /* sends a message to processor 1 so that he prints, etc.      */

   neighbor = proc_config[AZ_node] - 1;
   if ( proc_config[AZ_node] != 0) {
      mdwrap_iread((void *) &i, 0, &neighbor, &type, &request);
      mdwrap_wait((void *) &i, 0, &neighbor, &type, &cflag, &request);
   }
   AZ_printf_out("proc %d:\n",proc_config[AZ_node]);

   if (choice == AZ_input_form ) {
     if ( update != (int *) NULL) {
        AZ_printf_out("  N_update: %5d\n",N_update); AZ_printf_out("  update:");
        AZ_list_print(update, N_update, (double *) NULL , 0);
     }

     if (matrix == AZ_MSR_MATRIX) {
        AZ_printf_out("  bindx: ");
        AZ_list_print(bindx, bindx[N_update], (double *) NULL , 0);

        AZ_printf_out("  val:   ");
        AZ_list_print((int *) NULL , N_update, val , bindx[N_update]);
     }
     else if (matrix == AZ_VBR_MATRIX) {
        AZ_printf_out("  rpntr: ");
        AZ_list_print(rpntr, N_update+1, (double *) NULL , 0);
        if ( cpntr != (int *) NULL ) {
           AZ_printf_out("  cpntr: ");
           AZ_list_print(cpntr, N_update+1+ N_external, (double *) NULL , 0);
        }
        AZ_printf_out("  bpntr: ");
        AZ_list_print(bpntr, N_update+1, (double *) NULL , 0);
        AZ_printf_out("  bindx: ");
        AZ_list_print(bindx, bpntr[N_update], (double *) NULL , 0);
        AZ_printf_out("  indx:  ");
        AZ_list_print(indx, bpntr[N_update]+1, (double *) NULL , 0);
        AZ_printf_out("  val:   ");
        AZ_list_print((int *) NULL, indx[bpntr[N_update]], val, 0);
     }
   }
   else if (choice == AZ_global_mat ) {
     if ( matrix == AZ_MSR_MATRIX) {
        for (i = 0 ; i < N_update; i++ ) {
          ii = update_index[i];
          AZ_printf_out("   a(%d,%d) = %20.13e;\n",update[i]+of,update[i]+of,val[ii]);
          for (j = bindx[ii] ; j < bindx[ii+1] ; j++ ) {
            tj = AZ_find_simple(bindx[j], update_index, N_update, extern_index,
                              N_external,update,external);
            if (tj != -1) 
               AZ_printf_out("   a(%d,%d) = %20.13e;\n",update[i]+of,tj+of,val[j]);
            else (void) AZ_printf_err("col %d (= bindx[%d]) is undefined\n",
                                tj, j);
          }
        }
     }
     else if (matrix == AZ_VBR_MATRIX) {
        for (iblk= 0; iblk < N_update; iblk++) {
           new_iblk = update_index[iblk];

           m1 = rpntr[new_iblk+1] - rpntr[new_iblk];
 
           /* loop over blocks in the current block-row */
 
           for (ii = bpntr[new_iblk]; ii < bpntr[new_iblk+1]; ii++) {
              new_jblk = AZ_find_simple(bindx[ii], update_index, N_update, 
			 extern_index, N_external,update,external);
              if (new_jblk == -1) {
                 AZ_printf_out("local column %d not found\n",new_jblk);
                 exit(-1);
              }
              jblk = bindx[ii];
              ival =  indx[ii];

              n1 = cpntr[jblk+1] - cpntr[jblk];

              for (i = 0; i < m1; i++) {
                 for (j = 0; j < n1; j++)
                    (void) AZ_printf_out("   a(%d(%d),%d(%d)) = %20.13e;\n",update[iblk]+
				  of,i+of, new_jblk+of, j+of, val[ival+j*m1+i]);
              }
           }
        }
     }
   }
   else if (choice == AZ_explicit) {
     if ( matrix == AZ_MSR_MATRIX) {
        for (i = 0 ; i < N_update; i++ ) {
          if (update == NULL) tj = i+of;
          else tj = update[i] + of;
          AZ_printf_out("   a(%d,%d) = %20.13e;\n",tj,tj,val[i]);
          for (j = bindx[i] ; j < bindx[i+1] ; j++ ) {
               AZ_printf_out("   a(%d,%d) = %20.13e;\n",tj,bindx[j]+of,val[j]);
          }
        }
     }
     else if (matrix == AZ_VBR_MATRIX) {
        for (iblk = 0; iblk < N_update; iblk++) {
           if (update == NULL) tj = iblk+of;
           else tj = update[iblk] + of;

           m1 = rpntr[iblk+1] - rpntr[iblk];
 
           /* loop over blocks in the current block-row */
 
           for (ii = bpntr[iblk]; ii < bpntr[iblk+1]; ii++) {
              jblk = bindx[ii];
              ival =  indx[ii];
              n1 = (indx[ii+1]-ival)/m1;

              for (i = 0; i < m1; i++) {
                 for (j = 0; j < n1; j++)
                    (void) AZ_printf_out("   a(%d(%d),%d(%d)) = %20.13e;\n", tj, 
				  i+of, jblk+of, j+of, val[ival+j*m1+i]);
              }
 
           }
        }
     }
   }
   else (void) AZ_printf_err("AZ_matrix_out: output choice unknown\n");

   neighbor = proc_config[AZ_node] + 1;
   if ( proc_config[AZ_node] != proc_config[AZ_N_procs] - 1) 
      mdwrap_write((char *) &i, 0, neighbor, type, &cflag);

   i = AZ_gsum_int(i,proc_config);


}
int AZ_find_simple(int k, int *update_index, int N_update, int *extern_index,
   int N_external, int *update, int *external)
{
   int i;

   if (k < N_update) {
      for (i = 0 ; i < N_update ; i++ ) 
         if ( update_index[i] == k) return(update[i]);
   }
   else {
      for (i = 0 ; i < N_external ; i++ ) 
         if ( extern_index[i] == k) return(external[i]);
   }

   return(-1);
}
void AZ_list_print(int ivec[] , int length1, double dvec[],int length2)

{
   int i, count = 0;

   if (ivec != (int *) NULL) {
      for (i = 0 ; i < length1 ; i++ ) {
         AZ_printf_out("%8d ",ivec[i]); count += 8;
         if (count > 50) { count = 0; AZ_printf_out("\n         "); }
      }
   }
   else if (dvec != (double *) NULL) {
      for (i = 0 ; i < length1 ; i++ ) {
         AZ_printf_out("%8.1e ",dvec[i]); count += 8;
         if (count > 50) { count = 0; AZ_printf_out("\n         "); }
      }
      if (length2 != 0) { 
         AZ_printf_out("      -- "); count += 8;
         if (count > 50) { count = 0; AZ_printf_out("\n         "); }
      }

      for (i = length1+1 ; i < length2 ; i++ ) {
         AZ_printf_out("%8.1e ",dvec[i]); count += 8;
         if (count > 50) { count = 0; AZ_printf_out("\n         "); }
      }
   }
   AZ_printf_out("\n");
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
 
int AZ_compute_max_nz_per_row(AZ_MATRIX *Amat, int N, int Nb, int *largest_band)
{
/*******************************************************************************
 
  Compute (and return) the maximum number of nonzeros in 
  any matrix row. Also compute the largest band in the matrix 
  and return its value in '*largest_band'.
 
  Author:          Ray Tuminaro, SNL, 9222
  =======
 
  Return code:     int


  Parameter list:
  ===============
 
  matrix_type:     On input, indicates whether we have an MSR or VBR matrix.

  N:               On input, indicates order of matrix system.

  Nb:              On input, indicates number of blocks for a VBR matrix.

  bindx,bpntr      On input, matrix for which the maximum number of nonzeros
  cpntr:           in a row will be computed.

  largest_band:    On output, largest band in the matrix

*******************************************************************************/

   int i,j,largest = 0;
   int col, col_ptr;
   int kk,left_most,right_most;
   int matrix_type, *bindx, *bpntr, *cpntr;
 
   matrix_type = Amat->matrix_type;
   bindx       = Amat->bindx;
   *largest_band = -1;
   if (matrix_type == AZ_MSR_MATRIX) {
      for (i = 0 ; i < N; i++ ) {
         left_most  = i;
         right_most = i;
         j = bindx[i+1]-bindx[i];
         if (largest < j) largest = j;
         for (kk = bindx[i] ; kk < bindx[i+1] ; kk++ ) {
            if ( bindx[kk] < left_most  ) left_most  = bindx[kk];
            if ( bindx[kk] > right_most ) right_most = bindx[kk];
         }
         if (*largest_band <= right_most - left_most)
            *largest_band = right_most - left_most + 1;
      }
   }
   else if (matrix_type == AZ_VBR_MATRIX) {
      bpntr = Amat->bpntr;
      cpntr = Amat->cpntr;
      for  (i = 0; i < Nb ; i++) {
         j = 0;
         left_most  = cpntr[i];
         right_most = cpntr[i];
         for (col_ptr = bpntr[i]; col_ptr < bpntr[i+1]; col_ptr++ ) {
            col  = bindx[col_ptr];
            if (cpntr[col]   < left_most ) left_most  = cpntr[col];
            if (cpntr[col+1] > right_most) right_most = cpntr[col+1] -1;
            j   += (cpntr[col+1] - cpntr[col]);
         }
         if (*largest_band <= right_most - left_most) 
            *largest_band = right_most - left_most + 1;
         if (j > largest) largest = j;
      }
   }
   else {
      if ( (Amat->largest_band == -1) || (Amat->max_per_row == -1))
          AZ_matfree_Nnzs(Amat);
      *largest_band = Amat->largest_band;
      largest       = Amat->max_per_row - 1;
   }
   (*largest_band)++;
   largest++;
   return(largest);
}

/****************************************************************/
/****************************************************************/
/****************************************************************/
void AZ_check_block_sizes(int bindx[], int cpntr[], int Nrows, 
			  int *new_block)
{
/*
 * This routine is designed for mem_efficient_msr_2_vbr conversion.
 * It takes cpntr[] which is composed of a set of subsequences
 * (a subsequence is set of consecutive cpntr[] values that are 
 * the same) and splits these subsequences if they do not corrspond
 * to blocks in the MSR matrix (bindx,val). 
 *
 * Parameters
 * ==========
 *    bindx,          On input, an msr-LIKE matrix for which we wish
 *                    to discover the block structure. NONZEROS are
 *                    stored consecutively (row-wise) including the 
 *                    diagonal. The last column number in each row
 *                    is encoded as a negative (=  -1 - column).
 *
 *    cpntr[]         On input, cpntr[0] = 0. If the sparsity pattern 
 *                    of row i and row i-1 are identical, cpntr[i] = 
 *                    cpntr[i-1]. Otherwise, cpntr[i] = cpntr[i-1]+1.
 *                    On output, cpntr[0] = 0. If row i and row i-1 can
 *                    be grouped into a block, cpntr[i] = cpntr[i-1].
 *                    Otherwise, cpntr[i] = cpntr[i-1]+1.
 *
 *    Nrows           On input, number of rows in the matrix.
 * 
 *    *new_block      On input, number of unique values in cpntr[].
 *                    On output, number of unique values in the new cpntr[].
 */
 
   int start_blk_ptr, end_blk_ptr, start_blk_column, end_blk_column;
   int blk_name0, blk_name1, blk_name2, blk_name3;
   int prev_col, current_col, start_row, end_row;
   int row, tptr, max_col, last_nz_in_row, N_blks, current;

   max_col = Nrows-1; 

   if (Nrows == 0) return;

   start_blk_ptr = 0;
   row = 0;

   while (row < Nrows) {          /* for each row */
      start_row = start_blk_ptr;
      last_nz_in_row = 0;

      while(!last_nz_in_row) {    /* for each nonzero within a row */
         start_blk_column = bindx[start_blk_ptr];
         if (start_blk_column < 0) { 
            start_blk_column = -1 - start_blk_column; 
            last_nz_in_row = 1;
         }
         prev_col = start_blk_column;
         blk_name1   = cpntr[start_blk_column];

         end_blk_ptr = start_blk_ptr + 1;

         /* Find the end of the current block. If the column numbers */
         /* corresponding to this block are not all found (i.e.      */
         /* current_col != prev_col+1) cut the block off.            */

         while (!last_nz_in_row) {
            current_col = bindx[end_blk_ptr];
            if (current_col < 0) { 
                current_col = -1 - current_col; 
		last_nz_in_row = 1;
            }
            if (current_col != prev_col+1) break;
            if (cpntr[current_col] != blk_name1) break;
            end_blk_ptr++;
            prev_col = current_col;

         }
         end_blk_ptr--;

         /* check to see if the block name associated with the  */
         /* beginning and end of this block are different from  */
         /* the previous and next blocks. If they are the same  */
         /* create a new block name for this block.             */

         end_blk_column = bindx[end_blk_ptr];
         if (end_blk_column < 0) {
            end_blk_column = -1 - end_blk_column;
            last_nz_in_row = 1;
         }
         else last_nz_in_row = 0;

         blk_name2 = cpntr[end_blk_column];

         blk_name0 = -10; blk_name3 = -10;
         if (start_blk_column !=       0) blk_name0 = cpntr[start_blk_column-1];
         if (end_blk_column != max_col) blk_name3 = cpntr[end_blk_column+1];

         if (blk_name1 == blk_name0) {
            for ( tptr = start_blk_column; tptr <= end_blk_column; tptr++) 
               cpntr[tptr] = *new_block;
            (*new_block)++;
         }
         else if (blk_name2 == blk_name3) {
            for ( tptr = start_blk_column; tptr <= end_blk_column; tptr++) 
               cpntr[tptr] = *new_block;
            (*new_block)++;
         }
         start_blk_ptr = end_blk_ptr + 1;
      }
      end_row = end_blk_ptr;

      /* look for the next row which does not have the same */
      /* sparsity pattern as the previous row.              */

      row++;
      while ((row < Nrows) && (cpntr[row-1] == cpntr[row])) {
         row++;
         start_blk_ptr += (end_row - start_row + 1);
      }
   }

   /* Number each of the subsequences in cpntr (a subsequence  */
   /* is series of consecutive cpntr values that are the same) */
   /* uniquely starting from 0 and counting up.                */

   N_blks = 0;
   current = cpntr[0];
   cpntr[0] = N_blks;
   for (row = 1 ; row < Nrows ; row++ ) {
      if (cpntr[row] != current) {
         N_blks++;
         current = cpntr[row];
      }
      cpntr[row] = N_blks;
   }
   *new_block = N_blks;
}

#ifdef AZ_ENABLE_CAPTURE_MATRIX

void AZ_capture_matrix(AZ_MATRIX *Amat, int proc_config[],
                      int data_org[], double b[])

/*******************************************************************************

  Routine to capture matrix in a form acceptable to MATLAB.  Prints matrix out
  in i,j,a(i,j) format.  Currently implemented only for one processor.  Based
  on the routine AZ_output_matrix.
  Test to see if we should capture matrix, rhs and partitioning info
  in an ASCII data file.  If the file "AZ_write_matrix_now" exists in
  the current working directory, then the files 

  - AZ_capture_matrix.dat
  - AZ_capture_rhs.dat
  - AZ_capture_partition.dat (VBR only)
  
  will be appended with the current matrix in (i,j,val) format, the
  current RHS and the current partition information.  The existence
  of "AZ_write_matrix_now" is check each time.  Thus, capturing can
  be turned on and off at will during the run of a simulation.


  Author:          Michael A. Heroux, 9222, SNL.
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix.

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file params.txt).
  b:               Right hand side values.

*******************************************************************************/

{

  /* local variables */

  int  iblk_row, i, j, ib1, ib2, n1, jblk, m1, ipoint, jpoint;
  int  ival = 0;
  int  num_total_nonzeros;
  int  num_total_nodes, num_total_equations = 0;
  int  Proc, Num_Proc;
  double * x, * y;
  double * val;
  int * indx, * bindx, * rpntr,* cpntr, * bpntr;
  FILE *AZ_capture_flag;
  /********** execution begins **********/
  { 
    val = Amat->val;
    indx = Amat->indx;
    bindx =Amat->bindx;
    rpntr = Amat->rpntr;
    cpntr = Amat->cpntr;
    bpntr = Amat->bpntr;

  AZ_capture_flag = fopen("AZ_write_matrix_now","r");
  if(AZ_capture_flag)
    {

      Proc               = proc_config[AZ_node];
      Num_Proc           = proc_config[AZ_N_procs];
      if (Num_Proc != 1)
	AZ_printf_out("\nMatrix Capture Routine currently only works for one PE\n");
      else {

	AZ_print_sync_start(Proc, AZ_TRUE, proc_config);

	{ FILE *AZ_capt_matrix, *AZ_capture_rhs, *AZ_capture_partition;

	if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
	  num_total_nodes    = data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk];
	  num_total_equations = rpntr[num_total_nodes];
	  num_total_nonzeros = indx[bpntr[num_total_nodes]];

     /***** Print out the VBR partitioning information for the matrix *****/

	  AZ_capture_partition = fopen("AZ_capture_partition.dat","a");
	  
	  fprintf(AZ_capture_partition, "Start of partition\n");

	  for (i = 0; i < num_total_nodes + 1; i++)
	    fprintf(AZ_capture_partition,"%d\n", rpntr[i]);

	  fclose(AZ_capture_partition);

     /***** Print out the VBR i,j,a(i,j) information for the matrix *****/

	  AZ_capt_matrix = fopen("AZ_capture_matrix.dat","a");
	  fprintf(AZ_capt_matrix, "Start of VBR matrix\n");
	  fprintf(AZ_capt_matrix, "%d %d\n", 
		  num_total_equations, num_total_nonzeros);

	  /* loop over block rows */

	  for (iblk_row = 0; iblk_row < num_total_nodes; iblk_row++) {

	    /* the starting point row index of the current block */

	    ib1 = rpntr[iblk_row];      

	    /* number of rows in the current row block */

	    m1 = rpntr[iblk_row+1] - ib1;

	    /* starting index of current row block */

	    ival = indx[bpntr[iblk_row]];

	    /* loop over all the blocks in the current block-row */

	    for (j = bpntr[iblk_row]; j < bpntr[iblk_row+1]; j++) {
	      jblk = bindx[j];

	      /* the starting point column index of the current block */

	      ib2 = cpntr[jblk];

	      /* number of columns in the current block */

	      n1 = cpntr[jblk+1] - ib2;

	      for (jpoint = 0; jpoint < n1; jpoint++)
		for (ipoint = 0; ipoint < m1; ipoint++) {
		  fprintf(AZ_capt_matrix,"%d %d %22.16e\n", 
			  ib1+ipoint+1, 
			  ib2+jpoint+1, 
			  val[ival+jpoint*m1+ipoint]);
	  
		}

	      ival += m1*n1;
	    }
	  }
	  fclose(AZ_capt_matrix);
	}

	if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
	  num_total_equations = data_org[AZ_N_internal]+data_org[AZ_N_border];

	  num_total_nonzeros = bindx[num_total_equations]-1;

     /***** Print out the MSR i,j,a(i,j) information for the matrix *****/

	  AZ_capt_matrix = fopen("AZ_capture_matrix.dat","a");
#ifdef MB_MODIF
	  fprintf(AZ_capt_matrix, "* Start of MSR matrix\n");
#else
	  fprintf(AZ_capt_matrix, "Start of MSR matrix\n");
#endif
	  fprintf(AZ_capt_matrix, "%d %d\n", 
		  num_total_equations, num_total_nonzeros);
	  for (i = 0; i < num_total_equations; i++) {
	    fprintf(AZ_capt_matrix,"%d %d %22.16e\n", i+1, i+1, val[i]);
	    for (j = bindx[i]; j < bindx[i+1]; j++ )
	      fprintf(AZ_capt_matrix,"%d %d %22.16e\n", 
		      i+1, bindx[j]+1, val[j]);
	  }
	  fclose(AZ_capt_matrix);
	}

	else { /* Matrix-free mode so multiply by e_j vecs */
	  num_total_equations = data_org[AZ_N_internal]+data_org[AZ_N_border];

	  num_total_nonzeros = 0;

     /***** Print out the i,j,a(i,j) information for the matrix col-by-col *****/

	  AZ_capt_matrix = fopen("AZ_capture_matrix.dat","a");
#ifdef MB_MODIF
	  fprintf(AZ_capt_matrix, "* Start of Matrix-free matrix\n");
#else
	  fprintf(AZ_capt_matrix, "Start of Matrix-free matrix\nDim and NNZ info at bottom\n\n");
#endif
	  x = (double *) malloc(num_total_equations * sizeof(double));
	  y = (double *) malloc(num_total_equations * sizeof(double));
	  for (i = 0; i < num_total_equations; i++) x[i] = 0.0;
	  for (i = 0; i < num_total_equations; i++) y[i] = 0.0;
	  for (i = 0; i < num_total_equations; i++) {
	    /* For each column i, multiply by the vector e_i to extract the
	       coefficients from the ith column.  Then write these out. */
	    x[i] = 1.0;
	    Amat->matvec(x, y, Amat, proc_config);
	    for (j =0; j < num_total_equations; j++ )
	      if (y[j] !=0.0) {
		fprintf(AZ_capt_matrix,"%d %d %22.16e\n", 
			j+1, i+1, y[j]);
		num_total_nonzeros++;
	      }
	    x[i] = 0.0;
	  }
	  free ((void *) x);
	  free ((void *) y);
	  fprintf(AZ_capt_matrix, "%d %d\n", 
		  num_total_equations, num_total_nonzeros);
	  fclose(AZ_capt_matrix);
	}

     /***** Print out the RHS information for the matrix *****/

	AZ_capture_rhs = fopen("AZ_capture_rhs.dat","a");
	fprintf(AZ_capture_rhs, "Start of RHS\n");
	for (i = 0; i < num_total_equations; i++)
	  fprintf(AZ_capture_rhs,"%22.16e\n", b[i]);
	fclose(AZ_capture_rhs);

	AZ_print_sync_end(proc_config, AZ_TRUE);
	}
      }
      fclose(AZ_capture_flag);
    }
  }
} /* AZ_capture_matrix */
     
/*endif for 'ifdef AZ_ENABLE_CAPTURE_MATRIX' block */
#endif

   
struct submat_struct {
   int Nsub_rows, *sub_rows, Nsub_cols, *sub_cols;
};

int AZ_subMSR_getrow(int[], double[], int[], AZ_MATRIX *, int, int[], int);
void AZ_subMSR_matvec_mult (double *, double *, struct AZ_MATRIX_STRUCT *,
	int *);
int AZ_blockMSR_getrow(int[], double[], int[], AZ_MATRIX *, int, int[], 
	int);
void AZ_blockMSR_matvec_mult (double *, double *, struct AZ_MATRIX_STRUCT *,
	 int *);


AZ_MATRIX *AZ_submatrix_create(AZ_MATRIX *Amat, int Nsub_rows, int sub_rows[], 
				int Nsub_cols, int sub_cols[], int proc_config[])
{
	int i;
	
	AZ_MATRIX *sub_matrix;
	struct submat_struct *new_dataptr;

	
	sub_matrix = AZ_matrix_create(Nsub_rows);
	
	new_dataptr = (struct submat_struct *)malloc(sizeof(struct submat_struct));
	
	new_dataptr->Nsub_rows = Nsub_rows;
	new_dataptr->Nsub_cols = Nsub_cols;
	new_dataptr->sub_rows = (int *) malloc(Nsub_rows * sizeof(int));
	new_dataptr->sub_cols = (int *) malloc(Nsub_cols * sizeof(int));

	if ((new_dataptr->sub_rows == NULL) || (new_dataptr->sub_cols == NULL)) {
		AZ_printf_out("couldn't allocate memory for submatrix rows or cols\n");
		exit(-1);
	}

	for (i=0; i < Nsub_rows; i++)
		new_dataptr->sub_rows[i]=sub_rows[i];
	for (i = 0; i < Nsub_cols; i++)
		new_dataptr->sub_cols[i]=sub_cols[i];
	
	sub_matrix->bindx = Amat->bindx;
	sub_matrix->val = Amat->val;

	AZ_set_MATFREE(sub_matrix, new_dataptr, AZ_subMSR_matvec_mult);
	AZ_set_MATFREE_getrow(sub_matrix, new_dataptr, AZ_subMSR_getrow, NULL, 0, proc_config);

	return(sub_matrix);
}

void AZ_submatrix_destroy(AZ_MATRIX **submat)
{
	int *sub_rows, *sub_cols;
	struct submat_struct *data;

	data = (struct submat_struct *) AZ_get_matvec_data(*submat); 

	if (data != NULL) {
		sub_rows = data->sub_rows;
		sub_cols = data->sub_cols;
		
		free(sub_rows);
		free(sub_cols);
		free(data);
	}

	AZ_matrix_destroy(submat);
}
	
int  AZ_subMSR_getrow(int columns[], double values[], int row_lengths[],
											AZ_MATRIX *Amat, int N_requested_rows,
											int requested_rows[], int allocated_space)
{ 
	int    *bindx, i, j, count = 0, fullrow, count_row;
	int     Nsub_rows, Nsub_cols, *sub_rows, *sub_cols, subrow;
	int     fullcol, subcol;
	double *val;
	struct submat_struct *dataptr;
	
	bindx = Amat->bindx;
	val   = Amat->val;
	dataptr = (struct submat_struct *) AZ_get_matvec_data(Amat);
	sub_rows=dataptr->sub_rows;
	sub_cols=dataptr->sub_cols;
	Nsub_rows=dataptr->Nsub_rows;
	Nsub_cols=dataptr->Nsub_cols;
	
	for (i = 0; i < N_requested_rows; i++) {
		count_row = 0;
		subrow  = requested_rows[i];  /* row number in the local (sub) matrix */
		if (subrow < Nsub_rows)
			fullrow = sub_rows[subrow]; /* row number in the full matrix (of which
																		 this matrix is a part) */
		else {
			AZ_printf_out("getrow requested row %d of a submatrix with only %d rows\n",
						 subrow, Nsub_rows);
			exit(-1);
		}

		/* don't know the actual length of the returned row yet, but the length
			 of the row in the full matrix is a decent upper bound */
		row_lengths[i] = bindx[fullrow+1] - bindx[fullrow] + 1;
		if (count+row_lengths[i] > allocated_space) return(0);
		
		/* put in diagonal element if it exists */
		if (AZ_find_index(fullrow, sub_cols, Nsub_cols) >= 0) {
			columns[count + count_row  ] = subrow;
			values[count  + count_row++]  = val[fullrow];
		}


		/* off-diagonals */		
		for (j = bindx[fullrow] ; j < bindx[fullrow+1] ; j++) {
			fullcol = bindx[j];
			subcol = AZ_find_index(fullcol, sub_cols, Nsub_cols);
			if (subcol >= 0) {
				columns[count + count_row] = subcol;
				values[count  + count_row++] = val[j];
			}
		}
		/* now we know the real row length, so set it */
		row_lengths[i] = count_row;
		count += count_row;
	}
	return(1);
}

void AZ_subMSR_matvec_mult (double *b, double *c, struct AZ_MATRIX_STRUCT *Amat, 
														int proc_config[])
{
  double *val;
  int *data_org, *bindx, subrow, Nrows, Ncols, *rows, *cols;
  register int j, k, bindx_row;
  int          nzeros, fullrow, subcol;
	struct submat_struct *dataptr;


	dataptr = (struct submat_struct *) AZ_get_matvec_data(Amat);

	Nrows=dataptr->Nsub_rows;
	Ncols=dataptr->Nsub_cols;
	rows=dataptr->sub_rows;
	cols=dataptr->sub_cols;
  val = Amat->val;
  bindx = Amat->bindx;
  data_org = Amat->data_org;



  /* exchange boundary info */

  AZ_exchange_bdry(b, data_org, proc_config);

  for (subrow = 0; subrow < Nrows; subrow++) {
		fullrow = rows[subrow];  /* this row in the full matrix */

    /* compute diagonal contribution if there is one*/
		if (AZ_find_index(fullrow, cols, Ncols) >= 0)
			*c = val[fullrow] * b[subrow];
		else
			*c = 0.0;

    /* nonzero off diagonal contribution */

    bindx_row = bindx[fullrow];
    nzeros    = bindx[fullrow+1] - bindx_row;

		/* loop through the nonzeros in the full matrix and check to see if they're
			 in the submatrix.  If they are, add the appropriate product to the appropriate
			 entry */
    for (j = 0; j < nzeros; j++) {
      k   = bindx_row + j;
			subcol = AZ_find_index(bindx[k], cols, Ncols);
			if (subcol >= 0)				
				*c += val[k] * b[subcol];
    }
    c++;
  }
 
} /* AZ_subMSR_matvec_mult */



struct blockmat_struct {
	int Nblock_rows, Nblock_cols;
	int *Nsub_rows, **sub_rows, *Nsub_cols, **sub_cols, Nsub_mats;
	AZ_MATRIX **submat_list;
	int **submat_locs;
	int tot_Nrows;
};


AZ_MATRIX *AZ_blockmatrix_create(AZ_MATRIX **submat_list, int Nsub_mats, int **submat_locs, 
				int Nblock_rows, int Nblock_cols, int Nsub_rows[], int **sub_rows, int Nsub_cols[], 
				int **sub_cols, int proc_config[])
		 /* submat_locs is a Nsub_mats x 2 array of integers - the first column gives the block
				row of each submatrix and the second column gives the block column of each submatrix */
{
	int i, j, tot_Nrows;
	
	AZ_MATRIX *block_matrix;
	struct blockmat_struct *new_dataptr;

	
	tot_Nrows=0;
	for (i = 0; i < Nblock_rows; i++)
		tot_Nrows += Nsub_rows[i];
	
	block_matrix = AZ_matrix_create(tot_Nrows);
	
	new_dataptr = (struct blockmat_struct *)malloc(sizeof(struct blockmat_struct));
	
	new_dataptr->tot_Nrows   = tot_Nrows;
	new_dataptr->Nblock_rows = Nblock_rows;
	new_dataptr->Nblock_cols = Nblock_cols;
	new_dataptr->Nsub_mats   = Nsub_mats;

	/* allocate space for all the things we need to copy over */
	new_dataptr->submat_list = (AZ_MATRIX **) malloc(Nsub_mats * sizeof (AZ_MATRIX *));
	new_dataptr->submat_locs = (int **)  malloc(Nsub_mats   * sizeof(int *));
	new_dataptr->Nsub_rows   = (int *)  malloc(Nblock_rows * sizeof(int));
	new_dataptr->Nsub_cols   = (int *)  malloc(Nblock_cols * sizeof(int));
	new_dataptr->sub_rows    = (int **) malloc(Nblock_rows * sizeof(int *));
	new_dataptr->sub_cols    = (int **) malloc(Nblock_cols * sizeof(int *));
	
	if (new_dataptr->sub_cols == NULL) {
		AZ_printf_out("memory allocation error\n");
		exit(-1);
	}
	
	for (i = 0; i < Nsub_mats; i++) {
		new_dataptr->submat_list[i]=submat_list[i];
		new_dataptr->submat_locs[i]=(int *) malloc(2*sizeof(int));
		if (new_dataptr->submat_locs[i] == NULL) {
			AZ_printf_out("memory allocation error\n");
			exit(-1);
		}
		new_dataptr->submat_locs[i][0]=submat_locs[i][0];
		new_dataptr->submat_locs[i][1]=submat_locs[i][1];
	}

	for (i = 0; i < Nblock_rows; i++) {
		new_dataptr->Nsub_rows[i] = Nsub_rows[i];
		new_dataptr->sub_rows[i] = (int *) malloc(Nsub_rows[i] * sizeof(int));
		if (new_dataptr->sub_rows[i] == NULL) {
			AZ_printf_out("memory allocation error\n");
			exit(-1);
		}
		
		for (j = 0; j < Nsub_rows[i]; j++)
			new_dataptr->sub_rows[i][j]=sub_rows[i][j];
	}

	for (i = 0; i < Nblock_cols; i++) {
		new_dataptr->Nsub_cols[i] = Nsub_cols[i];
		new_dataptr->sub_cols[i] = (int *) malloc(Nsub_cols[i] * sizeof(int));
		if (new_dataptr->sub_cols[i] == NULL) {
			AZ_printf_out("memory allocation error\n");
			exit(-1);
		}
		
		for (j = 0; j < Nsub_cols[i]; j++)
			new_dataptr->sub_cols[i][j]=sub_cols[i][j];
	}
			
	
	AZ_set_MATFREE(block_matrix, new_dataptr, AZ_blockMSR_matvec_mult);
	AZ_set_MATFREE_getrow(block_matrix, new_dataptr, AZ_blockMSR_getrow, NULL, 0, proc_config);


	return(block_matrix);
}

void AZ_blockmatrix_destroy(AZ_MATRIX **blockmat)
{
	int **sub_rows, **sub_cols, *Nsub_rows, *Nsub_cols, **submat_locs;
	int i, Nblock_rows, Nblock_cols, Nsub_mats;
	struct blockmat_struct *data;
	AZ_MATRIX **submat_list;

	data = (struct blockmat_struct *) AZ_get_matvec_data(*blockmat);

	sub_rows = data->sub_rows;
	sub_cols = data->sub_cols;
	Nsub_rows = data->Nsub_rows;
	Nsub_cols = data->Nsub_cols;
	Nblock_rows = data->Nblock_rows;
	Nblock_cols = data->Nblock_cols;
	Nsub_mats = data->Nsub_mats;
	submat_locs = data->submat_locs;
	submat_list = data->submat_list;

	for (i=0; i < Nblock_rows; i++)
		free(sub_rows[i]);
	for (i=0; i < Nblock_cols; i++)
		free(sub_cols[i]);
	for (i=0; i < Nsub_mats; i++)
		free(submat_locs[i]);

	free(Nsub_rows);
	free(Nsub_cols);
	free(submat_list);
	free(submat_locs);
	free(sub_rows);
	free(sub_cols);
	free(data);
	AZ_matrix_destroy(blockmat);
}


int  AZ_blockMSR_getrow(int cols[], double vals[], int row_lengths[],
											AZ_MATRIX *Amat, int N_requested_rows,
											int requested_rows[], int allocated_space)
{ 
	int    i, j, k, ctr, count = 0, block_row, block_col, count_row;
	int     *Nsub_rows, **sub_rows, **sub_cols, full_req_row;
	int    local_req_row, tmp_row_len, *tmpcols, Nsub_mats, **submat_locs;
	int    tmp_alloc_space, max_nnz_per_row=500, tmp;
	double *tmpvals;
	struct blockmat_struct *dataptr;
	struct AZ_MATRIX_STRUCT *submat;

	dataptr = (struct blockmat_struct *) AZ_get_matvec_data(Amat);
	sub_rows=dataptr->sub_rows;
	sub_cols=dataptr->sub_cols;
	Nsub_rows=dataptr->Nsub_rows;
	Nsub_mats=dataptr->Nsub_mats;
	submat_locs=dataptr->submat_locs;

	tmpcols = (int *) malloc(max_nnz_per_row*sizeof(int));
	tmpvals = (double *) malloc(max_nnz_per_row*sizeof(double));
	if (tmpvals==NULL) {
		AZ_printf_out("memory allocation error\n");
		exit(-1);
	}
	tmp_alloc_space = max_nnz_per_row;


	for (i = 0; i < N_requested_rows; i++) {
		count_row = 0;
		full_req_row  = requested_rows[i];
		if (full_req_row > dataptr->tot_Nrows) {
			AZ_printf_out("Error: requested row %d of a matrix with %d rows\n", 
						 full_req_row+1, dataptr->tot_Nrows);
			exit(-1);
		}
			
		ctr=0;
		/* figure out which submatrix this row is in and which row it is 
			 within that submatrix */
		local_req_row=AZ_find_index(full_req_row, sub_rows[0],Nsub_rows[0]);
		while(local_req_row < 0) {
			ctr ++;
			local_req_row = AZ_find_index(full_req_row, sub_rows[ctr],
																		Nsub_rows[ctr]);			
		}

		block_row=ctr;  /* the block row where the requesed row is located */

		/* now loop through the submatrices and see which ones are in this block row: */
		for (j = 0; j < Nsub_mats; j++)
			if (submat_locs[j][0]==block_row) { /* this matrix contains part of the requested row */
				submat = dataptr->submat_list[j];
				/* figure out which block column this submatrix is in and get the relavent row */
				block_col=submat_locs[j][1];
				tmp = submat->getrow(tmpcols, tmpvals, &tmp_row_len, submat, 1, 
														 &local_req_row, tmp_alloc_space);
				/* in case we didn't allocate enough space for this row, keep checking the return
					 value of getrow until we get a positive result */
				while (tmp == 0) {
					free(tmpcols);
					free(tmpvals);
					max_nnz_per_row=max_nnz_per_row*2+1;
					tmp_alloc_space=max_nnz_per_row;
					tmpcols=(int *)malloc(tmp_alloc_space*sizeof(int));
					tmpvals=(double *)malloc(tmp_alloc_space*sizeof(double));
					tmp = submat->getrow(tmpcols, tmpvals, &tmp_row_len, submat, 1, 
															 &local_req_row, tmp_alloc_space);
				}

				for (k = 0; k < tmp_row_len; k++) {
					/* if we still have space in the output vector, go through and insert each
						 value from the submatrix row into the output vector */
					if (allocated_space > count+count_row) {
						cols[count+count_row]   = sub_cols[block_col][tmpcols[k]];
						vals[count+count_row++] = tmpvals[k];
					}
					else {
						free(tmpcols);
						free(tmpvals);
						return(0);
					}
				}
			}
		count += count_row;
		row_lengths[i]=count_row;
	}  

	free(tmpcols);
	free(tmpvals);
	return(1);
}

void AZ_blockMSR_matvec_mult (double *b, double *c, struct AZ_MATRIX_STRUCT *Amat, 
														int proc_config[])
{
	double *tmpb, *tmpc;
  int *data_org, Nrows;
  register int j;
  int          i;
	int *submat_loc, block_row, block_col, Nsub_mats, Nsub_rows, Nsub_cols;
	struct blockmat_struct *dataptr;
	AZ_MATRIX *submat;

  data_org = Amat->data_org;


  /* exchange boundary info */

  AZ_exchange_bdry(b, data_org, proc_config);


	dataptr = (struct blockmat_struct *) AZ_get_matvec_data(Amat);

	Nrows=dataptr->tot_Nrows;
	tmpb= (double *) malloc(Nrows * sizeof(double));
	tmpc= (double *) malloc(Nrows * sizeof(double));
	if (tmpc == NULL) {
		AZ_printf_out("memory allocation error\n");
		exit(-1);
	}

	for (i = 0; i < Nrows; i++)
		c[i] = 0.0;
	
	Nsub_mats= dataptr->Nsub_mats;

	/* loop through all the submatrices */
	for (i = 0; i < Nsub_mats; i++) {
		submat = dataptr->submat_list[i];
		submat_loc=dataptr->submat_locs[i];
		block_row = submat_loc[0];
		block_col = submat_loc[1];
		Nsub_rows = dataptr->Nsub_rows[block_row];
		Nsub_cols = dataptr->Nsub_cols[block_col];
		/* copy the rows of b that correspond to the columns of this
			 submatrix into a temp vector so we can run matvec on it */
		for (j = 0; j < Nsub_cols; j++)
			tmpb[j]=b[dataptr->sub_cols[block_col][j]];

		submat->matvec(tmpb, tmpc, submat, proc_config);

		/* now add the result to the appropriate rows of the output
			 vector */
		for (j = 0; j < Nsub_rows; j++)
			c[dataptr->sub_rows[block_row][j]] += tmpc[j];
	}
 
} /* AZ_blockMSR_matvec_mult */


void AZ_abs_matvec_mult(double *b, double *c,AZ_MATRIX *Amat,int proc_config[])


/******************************************************************************
  c = |A||b|:
  Sparse (square) overlapped matrix-vector multiply, using the  MSR
  data structure .

  Author:          Lydie Prevost, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============


  b:               Contains the vector b.

  c:               Contains the result vector c.



  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  Amat:            Structure used to represent the matrix (see file az_aztec.h
                   and Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

******************************************************************************/



{
  double *val, *x;
  int *data_org, *bindx;
  register int j, k, irow, bindx_row;
  int          N, nzeros, num_blks, bpoff, rpoff, iblk_row, ibpntr, 
               ibpntr_next = 0;
  int          m1, n1, irpntr, irpntr_next, ib1, ii, jj, iblk_size, jblk;
  double       *val_pntr, *c_pntr;
  int          *rpntr, *cpntr, *bpntr;

  val      = Amat->val; 
  bindx    = Amat->bindx;
  data_org = Amat->data_org;
  N = data_org[AZ_N_internal] + data_org[AZ_N_border];

  /* exchange boundary info */


  AZ_exchange_bdry(b, data_org, proc_config);

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
     for (irow = 0; irow < N; irow++) {
        *c = fabs(val[irow]) * fabs(b[irow]);
        bindx_row = bindx[irow];
        nzeros    = bindx[irow+1] - bindx_row;
        for (j = 0; j < nzeros; j++) {
           k   = bindx_row + j;
           *c += fabs(val[k]) * fabs(b[bindx[k]]);
        }
        c++;
     }
  }
  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
     cpntr = Amat->cpntr;
     rpntr = Amat->rpntr;
     bpntr = Amat->bpntr;
     num_blks = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

     bpoff    = *bpntr;
     rpoff    = *rpntr;
     val_pntr = val;
     for (j = 0; j < rpntr[num_blks] - rpoff; c[j++] = 0.0);
     irpntr_next = *rpntr++;
     bpntr++;
     c          -= rpoff;

     for (iblk_row = 0; iblk_row < num_blks; iblk_row++) {
        irpntr      = irpntr_next;
        irpntr_next = *rpntr++;
        ibpntr      = ibpntr_next;
        ibpntr_next = *bpntr++ - bpoff;
        c_pntr      = c + irpntr;
        m1          = irpntr_next - irpntr;

        /* loop over each block in the current row block */

        for (j = ibpntr; j < ibpntr_next; j++) {
           jblk = *(bindx+j);
           ib1  = *(cpntr+jblk);
           n1   = cpntr[jblk + 1] - ib1;
           iblk_size = m1*n1;
           x    = b + ib1;
           for (ii = 0; ii < m1; ii++) {
              for (jj = 0; jj < n1; jj++) {
                 c_pntr[ii] += fabs(val_pntr[n1*jj + ii])*fabs(x[jj]);
              }
           }
           val_pntr += iblk_size;
       }
     }


  }
  else {
     AZ_printf_out("Error: AZ_expected_value convergence options can only be done with MSR or VBR matrices\n");
     AZ_exit(1);
  }
 
} /* AZ_abs_matvec_mult */
