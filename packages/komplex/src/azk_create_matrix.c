/*
//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
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

#include <stdlib.h>
#include <stdio.h>
#include "az_aztec.h"
#include "azk_komplex.h"
/*#define DEBUG*/

/*! \file
\brief Creation routines for building Komplex matrices.

KOMPLEX is an add-on module to AZTEC that allows users to solve complex-valued
linear systems.

KOMPLEX solves a complex-valued linear system Ax = b by solving
an equivalent real-valued system of twice the dimension.  Specifically,
writing in terms of real and imaginary parts, we have

 \f[ (A_r + i*A_i)*(x_r + i*x_i) = (b_r + i*b_i) \f]

  or by separating into real and imaginary equations we have

\f[
  \left( \begin{array}{rr}
                                    A_r & -A_i\\
                                    A_i &  A_r
                             \end{array}
   \right)
   \left( \begin{array}{r}
                                    x_r\\
                                    x_i
                             \end{array}
   \right)
   =
   \left( \begin{array}{r}
                                    b_r\\
                                    b_i
                             \end{array}
   \right)
\f]
  which is a real-valued system of twice the size.  If we find xr and xi, we
  can form the solution to the original system as x = xr +i*xi.


KOMPLEX accept user linear systems in three forms with either global
or local index values.

1) The first form is true complex.  The user passes in an MSR or VBR
format matrix where the values are stored like Fortran complex
numbers.
Thus, the values array is of type double that is twice as long as the
number of complex values.  Each complex entry is stored with real part
followed by imaginary part (as in Fortran).

2) The second form stores real and imaginary parts separately, but the
pattern for each is identical.  Thus only the values of the imaginary
part are passed to the creation routines.

3) The third form accepts two real-valued matrices with no assumption
about the structure of the matrices.  Each matrix is multiplied by a
user-supplied complex constant.  This is the most general form.

Each of the above forms supports a global or local index set.  By this
we mean that the index values (stored in bindx) refer to the global
problem indices, or the local indices (for example after calling
AZ_transform).

*/


/*! \fn void AZK_create_matrix_c2k(int options[], double params[], 
                  int proc_config[], 
			   AZ_MATRIX *Amat_complex, AZ_MATRIX **Amat_komplex)

\brief Create Komplex matrix from Complex matrix.

Transforms a complex-valued matrix where double precision array hold 
the complex values of Amat_complex in Fortran complex format. 

\param options (In)
       Determines specific solution method and other parameters.
\param params (In)
       Drop tolerance and convergence tolerance info.
\param proc_config (In)
       Machine configuration.  proc_config[AZ_node] is the node
       number.  proc_config[AZ_N_procs] is the number of processors.
\param Amat_complex  (In)
       An AZ_MATRIX structure where Amat_complex->val contain the
       values of the complex matrix in Fortran complex format.

\param Amat_komplex (Out)
       Komplex version of matrix stored as an AZ_MATRIX structure.
*/


void AZK_create_matrix_c2k(int options[], double params[], int proc_config[], 
			   AZ_MATRIX *Amat_complex, AZ_MATRIX **Amat_komplex)
{
  AZ_KOMPLEX *pass_data;
  int i, j, ii, jj, kk;
  int nb_non_zeros = 0, n_non_zeros;
  int N, N_blk, N_equations, N_blk_equations = 0;
  int n_complex_nonzeros = 0, mat_type;
  int *bindx_complex, *bpntr_complex, *indx_complex, *rpntr_complex;
  double *val_complex;
  int *rpntr, *cpntr, *bpntr, *indx, *bindx;
  double *val;
  int *data_org_complex, *data_org;
  int *external, *update_index, *extern_index;
  int *update_complex, *global_bindx_complex;

#ifdef DEBUG
  int ione = 1;
#endif

  /* 

     NOTE:  This version works for Fortran complex storage of the values.
     It assumes the matrix index values are in global index space.
     This subroutine converts the real and imaginary parts of the complex
     matrix into a VBR matrix of twice the size where the block entries
     of of size 2 x 2.  The block entries are

     (1,1) =   real_part
     (1,2) = - imag_part
     (2,1) =   imag_part
     (2,2) =   real_part

  */

   /* Extract necessary data from data_org_complex */

  data_org_complex = Amat_complex->data_org;

  val_complex = Amat_complex->val;
  bindx_complex = Amat_complex->bindx;
  if (Amat_complex->has_global_indices)
    N = Amat_complex->N_local;
  else
    N = data_org_complex[AZ_N_internal] + data_org_complex[AZ_N_border];
  N_equations = 2 * N;

  if (Amat_complex->has_global_indices)
    mat_type = Amat_complex->matrix_type;
  else
    mat_type = data_org_complex[AZ_matrix_type];
  if (mat_type == AZ_VBR_MATRIX )
    {
      if (Amat_complex->has_global_indices)
	N_blk = Amat_complex->N_update;
      else
	N_blk = data_org_complex[AZ_N_int_blk] + data_org_complex[AZ_N_bord_blk];
      N_blk_equations = N_blk;
      nb_non_zeros = Amat_complex->bpntr[N_blk];
      n_complex_nonzeros = Amat_complex->indx[Amat_complex->bpntr[N_blk]];
    }
  else if (mat_type == AZ_MSR_MATRIX )
    {
      N_blk_equations = N;
      n_complex_nonzeros = Amat_complex->bindx[N]-1;
      nb_non_zeros = n_complex_nonzeros;
    }
  else
    {
      AZ_perror("Unsupported Matrix types");
    }

  if (!Amat_complex->has_global_indices)
    AZ_find_global_ordering (proc_config, Amat_complex, &global_bindx_complex,
			     &update_complex);
  else
    { 
      global_bindx_complex = bindx_complex;
      update_complex = Amat_complex->update;
    }

  n_non_zeros = n_complex_nonzeros*4;


  bindx = (int    *) AZ_allocate((nb_non_zeros+1)*sizeof(int));
  indx  = (int    *) AZ_allocate((nb_non_zeros+1)*sizeof(int));
  bpntr = (int    *) AZ_allocate((N_blk_equations+1)*sizeof(int));
  rpntr = (int    *) AZ_allocate((N_blk_equations+1)*sizeof(int));
  val   = (double *) AZ_allocate((n_non_zeros+1)*sizeof(double));
  if (val == NULL) AZ_perror("AZK_create_matrix_c2k: Out of memory.");

  /* Build VBR matrix with 2x-by-2x block entries */

  if (mat_type == AZ_VBR_MATRIX)
    {
      /* Extract pointers for real operator 
	 NOTE:  We are assuming that the real and imaginary parts have
	 same structure here!!!
      */

      bpntr_complex = Amat_complex->bpntr;
      indx_complex  = Amat_complex->indx;
      rpntr_complex  = Amat_complex->rpntr;


      /* global_bindx and bpntr are just copied */

      for (ii=0; ii < nb_non_zeros+1; ii++)
	bindx[ii] = global_bindx_complex[ii];
      for (ii=0; ii < N_blk_equations+1; ii++)
	bpntr[ii] = bpntr_complex[ii];

      /* Values in indx are 4x that of the real part */

      for (ii=0; ii < nb_non_zeros+1; ii++)
	indx[ii] = 4*indx_complex[ii];

      /* Values in rpntr are twice that of the real part */

      for (ii=0; ii<N_blk_equations+1; ii++) rpntr[ii] = 2*rpntr_complex[ii];

      /* Do for each row */
 
      for (ii=0; ii < N_blk_equations; ii++)
	{
	  
	  int istart, istop, nrow, ncol;
	  double *cur_complex, *cur_komplex;
	  istart = bpntr_complex[ii];
	  istop  = bpntr_complex[ii+1];
	  nrow = rpntr_complex[ii+1] - rpntr_complex[ii];

	  for (jj=istart; jj<istop ; jj++)
	    {
	      ncol = (indx_complex[jj+1]-indx_complex[jj])/nrow;
	      cur_complex = val_complex + 2*indx_complex[jj];
	      cur_komplex = val + indx[jj];
	      AZK_create_matrix_c2k_fill_entry (nrow, ncol, 
						cur_complex, cur_komplex);
	    }
	}
    }
  else
    {
      
      i = 0;
      j = 0;
      for (ii=0; ii < N_blk_equations; ii++)
	{
	  bpntr[ii] = i;
	  
	  /* Do lower triangle */
	  for (jj=global_bindx_complex[ii]; global_bindx_complex[jj]<ii && jj < 
		 global_bindx_complex[ii+1]; jj++)
	    {
	      bindx[i] = global_bindx_complex[jj];
	      indx[i] = j;
	      val[j  ] = val_complex[2*jj];
	      val[j+1] = val_complex[2*jj+1];
	      val[j+2] = - val[j+1];
	      val[j+3] = val[j];
	      i++;
	      j+=4;
	    }
	  /* Do main block diagonal */
	  bindx[i] = update_complex[ii];
	  indx[i] = j;
	  val[j  ] = val_complex[2*ii];
	  val[j+1] = val_complex[2*ii+1];
	  val[j+2] = - val[j+1];
	  val[j+3] = val[j];
	  i++;
	  j+=4;
	  
	  /* Do upper triangle */
	  kk = jj;
	  for (jj=kk; jj<global_bindx_complex[ii+1]; jj++)
	    {
	      bindx[i] = global_bindx_complex[jj];
	      indx[i] = j;
	      val[j  ] = val_complex[2*jj];
	      val[j+1] = val_complex[2*jj+1];
	      val[j+2] = - val[j+1];
	      val[j+3] = val[j];
	      i++;
	      j+=4;
	    }
	}
      bpntr[N_blk_equations] = n_complex_nonzeros;
      indx[i] = j;
      if (i != nb_non_zeros || i != n_complex_nonzeros)
	printf(" i = %d nb_non_zeros = %d j = %d n_complex_nonzeros = %d\n",
	       i, nb_non_zeros, j, n_complex_nonzeros);
      
      
      for (ii=0; ii<N_blk_equations+1; ii++) rpntr[ii] = ii*2;
    }

  /* Free temp space */

  if (!Amat_complex->has_global_indices)
    AZ_free((void *) global_bindx_complex);
 
  /* Transform komplex matrix from global to local indices */

  AZ_transform(proc_config, &external, bindx, val, update_complex,
	       &update_index, &extern_index, &data_org,
	       N_blk_equations, indx, bpntr, rpntr, &cpntr,
               AZ_VBR_MATRIX);
#ifdef DEBUG
  AZ_check_vbr(N_blk_equations, data_org[AZ_N_ext_blk], AZ_LOCAL, 
	       bindx, bpntr, cpntr, rpntr, proc_config);
#endif
  *Amat_komplex = AZ_matrix_create(N_equations);
  AZ_set_VBR((*Amat_komplex), rpntr, cpntr, bpntr, indx, bindx, val, data_org,
	     N_blk_equations, update_complex, AZ_LOCAL);
  pass_data = (AZ_KOMPLEX *) AZ_allocate(sizeof(AZ_KOMPLEX));
  if (pass_data == NULL) AZ_perror("AZK_create_matrix_c2k: Out of memory.");
  pass_data->external = external;
  pass_data->update_index = update_index;
  pass_data->extern_index = extern_index;
  pass_data->Form_of_Equations = AZK_True_Complex;

  if (Amat_complex->has_global_indices)
    pass_data->From_Global_Indices = 1;
  else
    pass_data->From_Global_Indices = 0;

  (*Amat_komplex)->aux_ptr = (void *) pass_data;  

/* end AZ_c2k
   */
}

void AZK_create_matrix_c2k_fill_entry(int nrow, int ncol, 
			   double *cur_complex, double *cur_komplex)
{
  int i, j;
  double *p11, *p12, *p21, *p22;
  p11 = cur_komplex;
  p21 = p11+1;
  p12 = cur_komplex+2*nrow;
  p22 = p12+1;
  
  for (j=0;j<2*ncol;j+=2)
    {
      for (i=0; i<2*nrow; i+=2)
	{
	  /*  printf("j = %d i = %d\n",j,i); */
	  p11[i] = *cur_complex;
	  p22[i] = *cur_complex++;
	  p12[i] = -*cur_complex;
	  p21[i] = *cur_complex++;
	}
      p11 += 4*nrow;
      p21 += 4*nrow;
      p12 += 4*nrow;
      p22 += 4*nrow;
    }
}

/*! \fn void AZK_create_matrix_g2k(int options[], double params[], int proc_config[], 
		  double c0r, double c0i, AZ_MATRIX *Amat_mat0, 
		  double c1r, double c1i, AZ_MATRIX *Amat_mat1, 
		  AZ_MATRIX **Amat_komplex)

\brief Create Komplex Matrix from General Matrix.

   Transforms a complex-valued Matrix

   (c0r+i*c0i)*A0 +(c1r+i*c1i)*A1)

   to a Komplex matrix.

\param options (In)
       Determines specific solution method and other parameters.
\param params (In)
       Drop tolerance and convergence tolerance info.
\param proc_config (In)
       Machine configuration.  proc_config[AZ_node] is the node
       number.  proc_config[AZ_N_procs] is the number of processors.
\param c0r (In)
       Real part of constant to be multiplied with first matrix.
\param c0i (In)
       Imaginary part of constant to be multiplied with first matrix.
\param Amat_mat0 (In)
       AZ_MATRIX object containing first real-valued matrix.
\param c1r (In)
       Real part of constant to be multiplied with second matrix.
\param c1i (In)
       Imaginary part of constant to be multiplied with second matrix.
\param Amat_mat1 (In)
       AZ_MATRIX object containing second real-valued matrix.

\param Amat_komplex (Out)
       Komplex version of matrix stored as an AZ_MATRIX structure.
*/


void AZK_create_matrix_g2k(int options[], double params[], int proc_config[], 
		  double c0r, double c0i, AZ_MATRIX *Amat_mat0, 
		  double c1r, double c1i, AZ_MATRIX *Amat_mat1, 
		  AZ_MATRIX **Amat_komplex)

{
  int ii, jj;
  int n_blk_nonzeros_mat0 = 0, n_blk_nonzeros_mat1 = 0;
  int N_equations = 0, N_blk_equations = 0;
  int is_VBR = 0;
  int n_nonzeros, n_blk_nonzeros;
  int N_update, N_mat0, N_mat1, N_blk_mat0, N_blk_mat1;
  int n_mat0_nonzeros, n_mat1_nonzeros;
  int mat_type_mat0, mat_type_mat1;

  int *data_org_mat0, *data_org_mat1;
  int *bindx_mat0, *bpntr_mat0, *indx_mat0, *rpntr_mat0;
  int *bindx_mat1, *bpntr_mat1, *indx_mat1, *rpntr_mat1;
  double *val_mat0, *val_mat1;
  int *global_bindx_mat0, *global_bindx_mat1;
  int *update_mat0, *update_mat1;

  int *cur_bindx, *cur_bindx_mat0, *cur_bindx_mat1;
  int *first_bindx_in_row, *last_bindx_in_row;
  int *map, *cur_map, *first_map_in_row;
  int nentries_mat0, nentries_mat1, nentries;
  int nrow, ncol, cur_global_blk_col;
  int *pntr_mat0, *pntr_mat1;

  double *mat0v, *mat1v, *komplex;
  int *new_bindx, *new_indx, *new_bpntr;

  int *bindx, *bpntr, *indx, *rpntr, *cpntr;
  int *external, *update_index, *extern_index, *data_org;
  double *val;
  AZ_KOMPLEX *pass_data;
#ifdef DEBUG
  int ione = 1;
  int n_neighbors;
  double fnrm, fnrm0, fnrm1;
#endif


  /* Extract necessary data from data_org */

  data_org_mat0 = Amat_mat0->data_org;
  data_org_mat1 = Amat_mat1->data_org;

  /* Extract pointers for mat0 and mat1 operators */
  
  val_mat0   = Amat_mat0->val;
  bindx_mat0 = Amat_mat0->bindx;
  bpntr_mat0 = Amat_mat0->bpntr;
  indx_mat0  = Amat_mat0->indx;
  rpntr_mat0 = Amat_mat0->rpntr;

  val_mat1   = Amat_mat1->val;
  bindx_mat1 = Amat_mat1->bindx;
  bpntr_mat1 = Amat_mat1->bpntr;
  indx_mat1  = Amat_mat1->indx;
  rpntr_mat1 = Amat_mat1->rpntr;

  if (Amat_mat0->has_global_indices)
    N_mat0 = Amat_mat0->N_local;
  else
    N_mat0 = data_org_mat0[AZ_N_internal] + data_org_mat0[AZ_N_border];
  if (Amat_mat1->has_global_indices)
    N_mat1 = Amat_mat1->N_local;
  else
    N_mat1 = data_org_mat1[AZ_N_internal] + data_org_mat1[AZ_N_border];
  if (N_mat0 != N_mat1)
    AZ_perror("Error: Dimensions of mat0 and mat1 part not equal");

  if (Amat_mat0->has_global_indices != Amat_mat1->has_global_indices)
    AZ_perror("Error: Indices of mat0 and mat1 must be of like type");

  if (Amat_mat0->has_global_indices)
    mat_type_mat0 = Amat_mat0->matrix_type;
  else
    mat_type_mat0 = data_org_mat0[AZ_matrix_type];

  if (Amat_mat1->has_global_indices)
    mat_type_mat1 = Amat_mat1->matrix_type;
  else
    mat_type_mat1 = data_org_mat1[AZ_matrix_type];

  if (mat_type_mat0 == AZ_VBR_MATRIX &&
      mat_type_mat1 == AZ_VBR_MATRIX )
    {
      is_VBR = 1;

      if (Amat_mat0->has_global_indices)
	N_blk_mat0 = Amat_mat0->N_update;
      else
	N_blk_mat0 = data_org_mat0[AZ_N_int_blk] + data_org_mat0[AZ_N_bord_blk];
      if (Amat_mat1->has_global_indices)
	N_blk_mat1 = Amat_mat1->N_update;
      else
	N_blk_mat1 = data_org_mat1[AZ_N_int_blk] + data_org_mat1[AZ_N_bord_blk];
      if (N_blk_mat0 != N_blk_mat1)
	AZ_perror("Error: Block dimensions of mat0 and mat1 part not equal\n");

      N_update = N_blk_mat0;
      N_equations = 2 * rpntr_mat0[N_update];
      N_blk_equations = N_blk_mat0;
      n_blk_nonzeros_mat0 = bpntr_mat0[N_blk_mat0];
      n_blk_nonzeros_mat1 = bpntr_mat1[N_blk_mat1];
      n_mat0_nonzeros = indx_mat0[bpntr_mat0[N_blk_mat0]];
      n_mat1_nonzeros = indx_mat1[bpntr_mat1[N_blk_mat1]];
    }
  else if (mat_type_mat0 == AZ_MSR_MATRIX &&
	   mat_type_mat1 == AZ_MSR_MATRIX)
    {
      is_VBR = 0;
      N_update = N_mat0;
      N_equations = 2 * N_update;
      N_blk_equations = N_update;
      n_mat0_nonzeros = bindx_mat0[N_update]-1;
      n_mat1_nonzeros = bindx_mat1[N_update]-1;
      n_blk_nonzeros_mat0 = n_mat0_nonzeros;
      n_blk_nonzeros_mat1 = n_mat1_nonzeros;
    }
  else
    {
      AZ_perror("Unsupported Matrix types");
    }
      
  /* After decoding structures we must have global index space in
     order to properly merge real and imaginary parts.  Must create
     an artificial one if matrix is in local index form.
  */

  if (!Amat_mat0->has_global_indices)
    AZ_find_global_ordering (proc_config, Amat_mat0, &global_bindx_mat0,
			     &update_mat0);
  else
    {
      global_bindx_mat0 = bindx_mat0;
      update_mat0 = Amat_mat0->update;
    }

  if (!Amat_mat1->has_global_indices)
    AZ_find_global_ordering (proc_config, Amat_mat1, &global_bindx_mat1,
			   &update_mat1); 
  else
    {
      global_bindx_mat1 = bindx_mat1;
      update_mat1 = Amat_mat1->update;
    }


#ifdef DEBUG
  if (proc_config[AZ_node] == 0)
    for (ii=0;ii<n_blk_nonzeros_mat1; ii++)
      printf("Processor %d of %d bindx_mat1[%d] = %d\n",
	     proc_config[AZ_node],proc_config[AZ_N_procs],ii,bindx_mat1[ii]);  
 
  {int n_mat0_nonzerosp1, n_mat1_nonzerosp1;
  if (mat_type_mat0 == AZ_VBR_MATRIX)
    {
      fnrm0 = ddot_(&n_mat0_nonzeros, (val_mat0), &ione, (val_mat0), &ione);
      fnrm1 = ddot_(&n_mat1_nonzeros, (val_mat1), &ione, (val_mat1), &ione);
    }
  else
    {
      n_mat0_nonzerosp1 = n_mat0_nonzeros + 1;
      n_mat1_nonzerosp1 = n_mat1_nonzeros + 1;
      fnrm0 = ddot_(&n_mat0_nonzerosp1, (val_mat0), &ione, (val_mat0), &ione) - 
	val_mat0[N_update]*val_mat0[N_update];
      fnrm1 = ddot_(&n_mat1_nonzerosp1, (val_mat1), &ione, (val_mat1), &ione) - 
	val_mat1[N_update]*val_mat1[N_update];
    }
  printf("Node %d:  Square of F-norm for mat0 matrix = %12.8e\n"
	 ,proc_config[AZ_node],fnrm0);
  printf("Node %d:  Number of nonzeros for mat0 matrix = %d\n"
	 ,proc_config[AZ_node],n_mat0_nonzeros);
  printf("Node %d:  Square of F-norm for mat1 matrix = %12.8e\n"
	 ,proc_config[AZ_node],fnrm1);
  printf("Node %d:  Number of nonzeros for mat1 matrix = %d\n"
	 ,proc_config[AZ_node],n_mat1_nonzeros);
  printf("Node %d:  Number of mat0 nonzeros =  %d\n",
	 proc_config[AZ_node], n_mat0_nonzeros);
  printf("Node %d:  Number of mat1 nonzeros =  %d\n",
	 proc_config[AZ_node], n_mat1_nonzeros);
  }
#endif

  /* Now allocate bindx to maximum possible size, will shrink later.
     Also, allocate a map vector for finding things in Amat_mat0 and Amat_mat1 */

  bindx = (int *)AZ_allocate((n_blk_nonzeros_mat0+n_blk_nonzeros_mat1)*sizeof(int));
  map   = (int *)AZ_allocate((n_blk_nonzeros_mat0+n_blk_nonzeros_mat1)*sizeof(int));
  if (map == NULL) AZ_perror("AZK_create_matrix_g2k: Out of memory.");


  /* 
     In the following loop set, bindx points to the merged set
     of column indices from the mat0 and mat1 parts.  Note that the
     mat0 indices were multiplied by 2 and the mat1 ones were  
     multiplied by 2 and shifted by one.  Thus, mat0 and mat1 parts
     will be interleaved as appropriate.  map records where the columns
     indices came from in global_bindx_mat0 and global_bindx_mat1
  */

  if (is_VBR)
    {
      pntr_mat0 = bpntr_mat0;
      pntr_mat1 = bpntr_mat1;
    }
  else
    {
      pntr_mat0 = global_bindx_mat0;
      pntr_mat1 = global_bindx_mat1;
    }
  
  n_nonzeros = 0;
  n_blk_nonzeros = 0;
  cur_bindx = bindx;
  cur_map = map;
  for (ii=0; ii < N_blk_equations; ii++)
    {
      cur_bindx_mat0 = global_bindx_mat0+pntr_mat0[ii];
      nentries_mat0  = pntr_mat0[ii+1] - pntr_mat0[ii];
      
      cur_bindx_mat1 = global_bindx_mat1+pntr_mat1[ii];
      nentries_mat1  = pntr_mat1[ii+1] - pntr_mat1[ii];
      
      if (is_VBR)
	{
	  nrow = rpntr_mat0[ii+1] - rpntr_mat0[ii];
	  if (nrow != rpntr_mat1[ii+1]-rpntr_mat1[ii])
	    AZ_perror("Error: Block partition of mat0 and mat1 part not equal\n");
	}
      else
	nrow = 1;
      
      nentries = nentries_mat0 + nentries_mat1;
      
      first_bindx_in_row = cur_bindx;
      
      first_map_in_row = cur_map;
      
      for (jj=0; jj < nentries_mat0 ; jj++)
	{
	  *cur_map++   = (int) (cur_bindx_mat0 - global_bindx_mat0);      
	  *cur_bindx++ = 2 * (*cur_bindx_mat0++);
	}
      
      for (jj=0; jj < nentries_mat1 ; jj++)
	{
	  *cur_map++   = (int) (cur_bindx_mat1 - global_bindx_mat1);
	  *cur_bindx++ = 1 + 2 * (*cur_bindx_mat1++);
	}
      
      if (!is_VBR) /* Add diagonal for MSR case */
	{
	  if (update_mat0[ii] != update_mat1[ii])
	    AZ_perror("Error: Update index  of mat0 and mat1 part not equal\n");
	  *cur_bindx++ = 2 * update_mat0[ii];
	  *cur_map++ = ii;
	  *cur_bindx++ = 1 + 2 * update_mat1[ii];
	  *cur_map++ = ii;
	  nentries += 2;
	}
      
      AZ_sort(first_bindx_in_row, nentries, first_map_in_row, NULL);
      /* 
	 Now we need to determine the number of komplex nonzero terms by
	 counting how many single mat0 parts, paired mat0/mat1 and single
	 mat1 parts there are.  We get both the block and point nonzero
	 count.
      */
      
      cur_bindx = first_bindx_in_row;
      cur_map = first_map_in_row;
      last_bindx_in_row = cur_bindx+nentries;
      
      while (cur_bindx<last_bindx_in_row)
	{
	  if (*cur_bindx%2==0)
	    {
	      if (is_VBR)
		n_nonzeros += 4*(indx_mat0[*cur_map+1]-indx_mat0[*cur_map]);
	      else
		n_nonzeros += 4;
	      
	      n_blk_nonzeros++;
	      cur_bindx++;
	      cur_map++;
	      if (cur_bindx<last_bindx_in_row && *cur_bindx%2==1) 
		{
		  cur_bindx++;
		  cur_map++;
		}
	    }
	  else
	    {
	      if (is_VBR)
		n_nonzeros += 4*(indx_mat1[*cur_map+1]-indx_mat1[*cur_map]);
	      else
		n_nonzeros += 4;
	      
	      n_blk_nonzeros++;
	      cur_bindx++;
	      cur_map++;
	    }
	} 
    }


  /* Now we can allocate the rest of the space needed for the komplex 
     matrix.
  */

  indx  = (int    *) AZ_allocate((n_blk_nonzeros+1)*sizeof(int));
  bpntr = (int    *) AZ_allocate((N_blk_equations+1)*sizeof(int));
  rpntr = (int    *) AZ_allocate((N_blk_equations+1)*sizeof(int));
  val   = (double *) AZ_allocate((n_nonzeros+1)*sizeof(double));
  
  if (val == NULL) AZ_perror("AZK_create_matrix_g2k: Out of memory.");
  
  /* Build VBR matrix with 2x-by-2x block entries */
  

  /* Values in rpntr are twice that of the mat0 part */
  
  if (is_VBR)
    for (ii=0; ii<N_blk_equations+1; ii++) rpntr[ii] = 2*rpntr_mat0[ii];
  else
    for (ii=0; ii<N_blk_equations+1; ii++) rpntr[ii] = ii*2;

  n_nonzeros = 0;
  n_blk_nonzeros = 0;
  cur_bindx = bindx;
  cur_map = map;
  new_bindx = bindx;
  new_indx = indx;
  new_bpntr = bpntr;
  
  *new_indx++  = 0;
  *new_bpntr++ = 0;
  
  for (ii=0; ii < N_blk_equations; ii++)
    {
      nentries_mat0  = pntr_mat0[ii+1] - pntr_mat0[ii];
      nentries_mat1  = pntr_mat1[ii+1] - pntr_mat1[ii];
      nentries = nentries_mat0 + nentries_mat1;

      if (!is_VBR) nentries +=2;     

      if (is_VBR)
	nrow = rpntr_mat0[ii+1] - rpntr_mat0[ii];
      else
	nrow = 1;
      
      first_bindx_in_row = cur_bindx;
      
      first_map_in_row = cur_map;
      
      cur_bindx = first_bindx_in_row;
      cur_map = first_map_in_row;
      last_bindx_in_row = cur_bindx+nentries;
      if (proc_config[AZ_node] == -2)
	printf("Processor %d of %d nentries[%d] = %d mat0 = %d mat1 = %d.\n",
	       proc_config[AZ_node],proc_config[AZ_N_procs],ii,nentries,
	       nentries_mat0,nentries_mat1);       
      while (cur_bindx<last_bindx_in_row)
	{
	  if (*cur_bindx%2==0) /* Indicates there is a mat0 part */
	    {
      if (proc_config[AZ_node] == -2)
	printf("Processor %d of %d cur_bindx[%d] = %d in mat0 part.\n",
	       proc_config[AZ_node],proc_config[AZ_N_procs],
	       (int) (cur_bindx-first_bindx_in_row),*cur_bindx);
	      if (is_VBR)
		{
		  ncol = (indx_mat0[*cur_map+1] - indx_mat0[*cur_map])/nrow;
		  cur_global_blk_col = global_bindx_mat0[*cur_map];
		  mat0v = val_mat0 + indx_mat0[*cur_map++];
		}
	      else
		{
		  ncol = 1;
		  if (*cur_map == ii)
		    cur_global_blk_col = update_mat0[ii];
		  else
		    cur_global_blk_col = global_bindx_mat0[*cur_map];

		  mat0v = val_mat0 + *cur_map++;
		}
	      cur_bindx++;
	      
	      if (cur_bindx<last_bindx_in_row && *cur_bindx%2==1) 
		{
      if (proc_config[AZ_node] == -2)
	printf("Processor %d of %d cur_bindx[%d] = %d in mat1 of mat0 part.\n",
	       proc_config[AZ_node],proc_config[AZ_N_procs],
	       (int) (cur_bindx-first_bindx_in_row),*cur_bindx);
		  cur_bindx++;
		  if (is_VBR)
		      mat1v = val_mat1 + indx_mat1[*cur_map++];
		  else
		      mat1v = val_mat1 + *cur_map++;
		}
	      else
		mat1v = NULL;
	  
	    }

	  else           /* Otherwise there is only a mat1 part */
	    {
      if (proc_config[AZ_node] == -2)
	printf("Processor %d of %d cur_bindx[%d] = %d in mat1 part.\n",
	       proc_config[AZ_node],proc_config[AZ_N_procs],
	       (int) (cur_bindx-first_bindx_in_row),*cur_bindx);
	      mat0v = NULL;
	      
	      if (is_VBR)
		{
		  ncol = (indx_mat1[*cur_map+1] - indx_mat1[*cur_map])/nrow;
		  cur_global_blk_col = global_bindx_mat1[*cur_map];
		  mat1v = val_mat1 + indx_mat1[*cur_map++];
		}
	      else
		{
		  ncol = 1;
		  if (*cur_map == ii)
		    cur_global_blk_col = update_mat1[ii];
		  else
		    cur_global_blk_col = global_bindx_mat1[*cur_map];

		  mat1v = val_mat1 + *cur_map++;
		}
	      *cur_bindx++;
	    }
	  komplex = val + n_nonzeros;
	  n_nonzeros += 4*nrow*ncol;
	  *new_indx++  = n_nonzeros;
	  *new_bindx++ = cur_global_blk_col;
	  n_blk_nonzeros++;
	  AZK_create_matrix_g2k_fill_entry (nrow, ncol, c0r, c0i, mat0v, 
				 c1r, c1i, mat1v, komplex);
	  
	}
      *new_bpntr++ = n_blk_nonzeros;
      /*   printf("Processor %d of %d bpntr[%d] = %d.\n",
	   proc_config[AZ_node],proc_config[AZ_N_procs],ii,n_blk_nonzeros) ; */
    }

  /* Resize and destroy space */
  bindx = (int *) AZ_realloc((void *) bindx, (n_blk_nonzeros+1)*sizeof(int));
  AZ_free ((void *) map);

#ifdef DEBUG
  fnrm = ddot_(&n_nonzeros, val, &ione, val, &ione);
  printf("Node %d:  Square of F-norm for local matrix = %12.8e\n"
	 ,proc_config[AZ_node],fnrm);
  printf("Node %d:  Number of nonzeros for local matrix = %d\n"
	 ,proc_config[AZ_node],n_nonzeros);
  printf("Node %d:  Building RHS and Initial guess\n",proc_config[AZ_node]);
#endif

  /* Free temp space */

  if (!Amat_mat0->has_global_indices)
    {
      AZ_free((void *) global_bindx_mat0);
      AZ_free((void *) global_bindx_mat1);
      AZ_free((void *) update_mat1);
    }
  
  /* Transform komplex matrix from global to local indices */

  AZ_transform(proc_config, &external, bindx, val, update_mat0,
            &update_index, &extern_index, &data_org,
            N_blk_equations, indx, bpntr, rpntr,  &cpntr,
               AZ_VBR_MATRIX);
#ifdef DEBUG
  n_neighbors = data_org[AZ_N_neigh];

  printf("Processor %d of %d has %d neighbors.\n",
	 proc_config[AZ_node],proc_config[AZ_N_procs],n_neighbors) ;

  for (ii=0; ii<n_neighbors; ii++)
    {
      printf("Processor %d of %d receiving %d elements from Processor %d.\n",
	     proc_config[AZ_node],proc_config[AZ_N_procs],
	     data_org[AZ_rec_length+ii],data_org[AZ_neighbors+ii]) ;
      printf("Processor %d of %d sending %d elements to Processor %d.\n",
	     proc_config[AZ_node],proc_config[AZ_N_procs],
	     data_org[AZ_send_length+ii],data_org[AZ_neighbors+ii]) ;
    }

  /* Note:  This assumes mat0,mat1 interleaving.  Must change if we do more
     sophisticated reorderings. */
  printf("Processor %d of %d sending a total of %d elements.\n",
	 proc_config[AZ_node],proc_config[AZ_N_procs],
	 data_org[AZ_total_send]);

  for (ii=0; ii<data_org[AZ_total_send]/2; ii++)
    {
      printf("Processor %d of %d sending elements %d and %d.\n",
	     proc_config[AZ_node],proc_config[AZ_N_procs],
	     data_org[AZ_send_list+2*ii  ],data_org[AZ_send_list+2*ii+1]) ;
    }
  AZ_check_vbr(N_blk_equations, data_org[AZ_N_ext_blk], AZ_LOCAL, 
	       bindx, bpntr, cpntr, rpntr, proc_config);
#endif

  *Amat_komplex = AZ_matrix_create(N_equations);
  AZ_set_VBR((*Amat_komplex), rpntr, cpntr, bpntr, indx, bindx, val, data_org,
	     N_blk_equations, update_mat0, AZ_LOCAL);
  pass_data = (AZ_KOMPLEX *) AZ_allocate(sizeof(AZ_KOMPLEX));
  if (pass_data == NULL) AZ_perror("AZK_create_matrix_g2k: Out of memory.");
  pass_data->external = external;
  pass_data->update_index = update_index;
  pass_data->extern_index = extern_index;
  pass_data->Form_of_Equations = AZK_Komplex_General;

  if (Amat_mat0->has_global_indices)
    pass_data->From_Global_Indices = 1;
  else
    pass_data->From_Global_Indices = 0;

  (*Amat_komplex)->aux_ptr = (void *) pass_data;


  /* end AZ_complex_to_komplex
   */
}

void AZK_create_matrix_g2k_fill_entry(int nrow, int ncol, 
			   double c0r, double c0i, double *mat0v, 
			   double c1r, double c1i, double *mat1v, 
			   double *komplex)
{
  int i, j;
  double *p11, *p12, *p21, *p22;
  double m0, m1;
  p11 = komplex;
  p21 = p11+1;
  p12 = komplex+2*nrow;
  p22 = p12+1;


  if (mat0v != NULL && mat1v != NULL)
    {
      for (j=0;j<2*ncol;j+=2)
	{
	  for (i=0; i<2*nrow; i+=2)
	    {
	      /*  printf("j = %d i = %d\n",j,i); */
	      m0 = *mat0v++;
	      m1 = *mat1v++;
	      p11[i] = c0r*m0 + c1r*m1;
	      p22[i] = p11[i];
	      p21[i] = c1i*m1 + c0i*m0;
	      p12[i] = -p21[i];
	    }
	  p11 += 4*nrow;
	  p21 += 4*nrow;
	  p12 += 4*nrow;
	  p22 += 4*nrow;
	}
    }
  else if (mat0v != NULL)
    {
      for (j=0;j<2*ncol;j+=2)
	{
	  for (i=0; i<2*nrow; i+=2)
	    {
	      /*  printf("j = %d i = %d\n",j,i); */
	      m0 = *mat0v++;
	      m1 = 0.0;
	      p11[i] = c0r*m0 + c1r*m1;
	      p22[i] = p11[i];
	      p21[i] = c1i*m1 + c0i*m0;
	      p12[i] = -p21[i];
	    }
	  p11 += 4*nrow;
	  p21 += 4*nrow;
	  p12 += 4*nrow;
	  p22 += 4*nrow;
	}
    }
  else /* if (mat1v != NULL) */
    {
      for (j=0;j<2*ncol;j+=2)
	{
	  for (i=0; i<2*nrow; i+=2)
	    {
	      /*  printf("j = %d i = %d\n",j,i); */
	      m0 = 0.0;
	      m1 = *mat1v++;
	      p11[i] = c0r*m0 + c1r*m1;
	      p22[i] = p11[i];
	      p21[i] = c1i*m1 + c0i*m0;
	      p12[i] = -p21[i];
	    }
	  p11 += 4*nrow;
	  p21 += 4*nrow;
	  p12 += 4*nrow;
	  p22 += 4*nrow;
	}
    }
}
 
/*! \fn void AZK_create_matrix_ri2k(int options[], double params[], int proc_config[], 
			    AZ_MATRIX *Amat_real, double *val_imag, 
			    AZ_MATRIX **Amat_komplex)

\brief Create Komplex Matrix from Real and Imaginary Parts.

   Transforms a complex-valued matrix

         (Ar +i*Ai) 

   where double precision arrays hold the real and imaginary parts
separately.  The pattern of the imaginary part matches the real part. Thus no
structure for the imaginary part is passed in.

\param options (In)
       Determines specific solution method and other parameters.
\param params (In)
       Drop tolerance and convergence tolerance info.
\param proc_config (In)
       Machine configuration.  proc_config[AZ_node] is the node
       number.  proc_config[AZ_N_procs] is the number of processors.
\param Amat_real (In)
       AZ_MATRIX object containing real matrix.
\param val_imag (In)
       Double arrya containing the values ONLY for imaginary matrix.

\param Amat_komplex (Out)
       Komplex version of matrix stored as an AZ_MATRIX structure.

*/ 
void AZK_create_matrix_ri2k(int options[], double params[], int proc_config[], 
			    AZ_MATRIX *Amat_real, double *val_imag, 
			    AZ_MATRIX **Amat_komplex)
{
  AZ_KOMPLEX *pass_data;
  int i, j, ii, jj, kk;
  int nb_non_zeros = 0, n_non_zeros;
  int N, N_blk, N_equations, N_blk_equations = 0;
  int n_real_nonzeros = 0, mat_type;
  int *bindx_real, *bpntr_real, *indx_real, *rpntr_real;
  double *val_real;
  int *rpntr, *cpntr, *bpntr, *indx, *bindx;
  double *val;
  int *data_org_real, *data_org;
  int *external, *update_index, *extern_index;
  int *update_real, *global_bindx_real;

#ifdef DEBUG
  int ione = 1;
#endif
  /* 

     NOTE:  This version works for Fortran real storage of the values.
     It assumes the matrix index values are in global index space.
     This subroutine converts the real and imaginary parts of the real
     matrix into a VBR matrix of twice the size where the block entries
     of of size 2 x 2.  The block entries are

     (1,1) =   real_part
     (1,2) = - imag_part
     (2,1) =   imag_part
     (2,2) =   real_part

  */

   /* Extract necessary data from data_org_real */

  data_org_real = Amat_real->data_org;

  val_real = Amat_real->val;
  bindx_real = Amat_real->bindx;
  if (Amat_real->has_global_indices)
    N = Amat_real->N_local;
  else
    N = data_org_real[AZ_N_internal] + data_org_real[AZ_N_border];
  N_equations = 2 * N;

  if (Amat_real->has_global_indices)
    mat_type = Amat_real->matrix_type;
  else
    mat_type = data_org_real[AZ_matrix_type];
  if (mat_type == AZ_VBR_MATRIX )
    {
      if (Amat_real->has_global_indices)
	N_blk = Amat_real->N_update;
      else
	N_blk = data_org_real[AZ_N_int_blk] + data_org_real[AZ_N_bord_blk];
      N_blk_equations = N_blk;
      nb_non_zeros = Amat_real->bpntr[N_blk];
      n_real_nonzeros = Amat_real->indx[Amat_real->bpntr[N_blk]];
    }
  else if (mat_type == AZ_MSR_MATRIX )
    {
      N_blk_equations = N;
      n_real_nonzeros = Amat_real->bindx[N]-1;
      nb_non_zeros = n_real_nonzeros;
    }
  else
    {
      AZ_perror("Unsupported Matrix types");
    }

  if (!Amat_real->has_global_indices)
    AZ_find_global_ordering (proc_config, Amat_real, &global_bindx_real,
			     &update_real);
  else
    { 
      global_bindx_real = bindx_real;
      update_real = Amat_real->update;
    }

  n_non_zeros = n_real_nonzeros*4;


  bindx = (int    *) AZ_allocate((nb_non_zeros+1)*sizeof(int));
  indx  = (int    *) AZ_allocate((nb_non_zeros+1)*sizeof(int));
  bpntr = (int    *) AZ_allocate((N_blk_equations+1)*sizeof(int));
  rpntr = (int    *) AZ_allocate((N_blk_equations+1)*sizeof(int));
  val   = (double *) AZ_allocate((n_non_zeros+1)*sizeof(double));
  if (val == NULL) AZ_perror("AZK_create_matrix_ri2k: Out of memory.");

  /* Build VBR matrix with 2x-by-2x block entries */

  if (mat_type == AZ_VBR_MATRIX)
    {
      /* Extract pointers for real operator 
	 NOTE:  We are assuming that the real and imaginary parts have
	 same structure here!!!
      */

      bpntr_real = Amat_real->bpntr;
      indx_real  = Amat_real->indx;
      rpntr_real  = Amat_real->rpntr;


      /* global_bindx and bpntr are just copied */

      for (ii=0; ii < nb_non_zeros+1; ii++)
	bindx[ii] = global_bindx_real[ii];
      for (ii=0; ii < N_blk_equations+1; ii++)
	bpntr[ii] = bpntr_real[ii];

      /* Values in indx are 4x that of the real part */

      for (ii=0; ii < nb_non_zeros+1; ii++)
	indx[ii] = 4*indx_real[ii];

      /* Values in rpntr are twice that of the real part */

      for (ii=0; ii<N_blk_equations+1; ii++) rpntr[ii] = 2*rpntr_real[ii];
 
      /* Do for each row */
 
      for (ii=0; ii < N_blk_equations; ii++)
	{
	  
	  int istart, istop, nrow, ncol;
	  double *realv, *imagv, *komplex;
	  istart = bpntr_real[ii];
	  istop  = bpntr_real[ii+1];
	  nrow = rpntr_real[ii+1] - rpntr_real[ii];
	    
	  for (jj=istart; jj<istop ; jj++)
	    {
	      ncol = (indx_real[jj+1]-indx_real[jj])/nrow;
	      realv = val_real + indx_real[jj];
	      imagv = val_imag + indx_real[jj];
	      komplex = val + indx[jj];
	      AZK_create_matrix_ri2k_fill_entry (nrow, ncol, realv, imagv, komplex);
	    }
	}
    }
  else
    {
      
      i = 0;
      j = 0;
      for (ii=0; ii < N_blk_equations; ii++)
	{
	  bpntr[ii] = i;
	  
	  /* Do lower triangle */
	  for (jj=global_bindx_real[ii]; global_bindx_real[jj]<ii && jj < 
		 global_bindx_real[ii+1]; jj++)
	    {
	      bindx[i] = global_bindx_real[jj];
	      indx[i] = j;
	      val[j  ] = val_real[jj];
	      val[j+1] = val_imag[jj];
	      val[j+2] = - val[j+1];
	      val[j+3] = val[j];
	      i++;
	      j+=4;
	    }
	  /* Do main block diagonal */
	  bindx[i] = update_real[ii];
	  indx[i] = j;
	  val[j  ] = val_real[ii];
	  val[j+1] = val_imag[ii];
	  val[j+2] = - val[j+1];
	  val[j+3] = val[j];
	  i++;
	  j+=4;
	  
	  /* Do upper triangle */
	  kk = jj;
	  for (jj=kk; jj<global_bindx_real[ii+1]; jj++)
	    {
	      bindx[i] = global_bindx_real[jj];
	      indx[i] = j;
	      val[j  ] = val_real[jj];
	      val[j+1] = val_imag[jj];
	      val[j+2] = - val[j+1];
	      val[j+3] = val[j];
	      i++;
	      j+=4;
	    }
	}
      bpntr[N_blk_equations] = n_real_nonzeros;
      indx[i] = j;
      if (i != nb_non_zeros || i != n_real_nonzeros)
	printf(" i = %d nb_non_zeros = %d j = %d n_real_nonzeros = %d\n",
	       i, nb_non_zeros, j, n_real_nonzeros);
      
      
      for (ii=0; ii<N_blk_equations+1; ii++) rpntr[ii] = ii*2;
    }

  /* Free temp space */

  if (!Amat_real->has_global_indices)
    AZ_free((void *) global_bindx_real);

  /* Transform komplex matrix from global to local indices */

  AZ_transform(proc_config, &external, bindx, val, update_real,
	       &update_index, &extern_index, &data_org,
	       N_blk_equations, indx, bpntr, rpntr, &cpntr,
               AZ_VBR_MATRIX);
#ifdef DEBUG
  AZ_check_vbr(N_blk_equations, data_org[AZ_N_ext_blk], AZ_LOCAL, 
	       bindx, bpntr, cpntr, rpntr, proc_config);
#endif
  *Amat_komplex = AZ_matrix_create(N_equations);
  AZ_set_VBR((*Amat_komplex), rpntr, cpntr, bpntr, indx, bindx, val, data_org,
	     N_blk_equations, update_real, AZ_LOCAL);
  pass_data = (AZ_KOMPLEX *) AZ_allocate(sizeof(AZ_KOMPLEX));
  if (pass_data == NULL) AZ_perror("AZK_create_matrix_ri2k: Out of memory.");
  pass_data->external = external;
  pass_data->update_index = update_index;
  pass_data->extern_index = extern_index;
  pass_data->Form_of_Equations = AZK_Komplex_Same_Structure;

  if (Amat_real->has_global_indices)
    pass_data->From_Global_Indices = 1;
  else
    pass_data->From_Global_Indices = 0;

  (*Amat_komplex)->aux_ptr = (void *) pass_data;  

/* end AZ_ri2k
   */
}

void AZK_create_matrix_ri2k_fill_entry(int nrow, int ncol, 
			   double *realv, double *imagv, double *komplex)
{
  int i, j;
  double *p11, *p12, *p21, *p22;
  p11 = komplex;
  p21 = p11+1;
  p12 = komplex+2*nrow;
  p22 = p12+1;
  
  for (j=0;j<2*ncol;j+=2)
    {
      for (i=0; i<2*nrow; i+=2)
	{
	  /*  printf("j = %d i = %d\n",j,i); */
	  p11[i] = *realv;
	  p22[i] = *realv++;
	  p12[i] = -*imagv;
	  p21[i] = *imagv++;
	}
      p11 += 4*nrow;
      p21 += 4*nrow;
      p12 += 4*nrow;
      p22 += 4*nrow;
    }
}
 
