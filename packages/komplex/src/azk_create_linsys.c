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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "azk_komplex.h"
/*! \file 
\brief Creation routines for building Komplex systems.

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
format matrix where the values are stored like Fortran complex numbers.
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


/*! \fn void AZK_create_linsys_c2k(double *xc, double *bc,
                      int *options, double *params, int *proc_config,
                      AZ_MATRIX *Amat_complex,
                      double **x, double **b, AZ_MATRIX **Amat_komplex)

\brief Create Komplex System from Complex System.

Transforms a complex-valued system

         Amat_complex * xc = bc

   where double precision arrays hold the complex values of Amat_complex, xc
   and bc in Fortran complex format, i.e., if dimension of complex system is N
   then xc is of length 2*N and the first complex value is stored with the
   real part in xc[0] and the imaginary part in xc[1] and so on.

\param xc (In) 
       Contains the complex initial guess/solution vector with the
       real/imag parts interleaved as in Fortran complex format.
\param bc (In) 
       RHS in Fortran complex format.
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

\param x (Out)
       Komplex version of initial guess and solution.
\param b (Out)
       Komplex version of RHS.
\param Amat_komplex (Out)
       Komplex version of matrix stored as an AZ_MATRIX structure.

*/
void AZK_create_linsys_c2k(double *xc, double *bc,
                      int *options, double *params, int *proc_config,
                      AZ_MATRIX *Amat_complex,
                      double **x, double **b, AZ_MATRIX **Amat_komplex)
{

  /* First executable statement */

  /* Build Komplex matrix */

  AZK_create_matrix_c2k(options, params, proc_config, 
				Amat_complex, Amat_komplex);

  /* Build RHS and initial guess */

  AZK_create_vector_c2k( options, params, proc_config, 
				 (*Amat_komplex), xc, x);
  AZK_create_vector_c2k( options, params, proc_config, 
				 (*Amat_komplex), bc, b);

  /* Permute K system for better numerical stability */

  AZK_permute_ri( options, params, proc_config, (*b), (*Amat_komplex));
  return;
}

/*! \fn void AZK_create_linsys_g2k(double *xr, double *xi, double *br, 
                      double *bi, 
			       int *options, double *params, int *proc_config,
			       double c0r, double c0i, AZ_MATRIX *Amat_mat0, 
			       double c1r, double c1i, AZ_MATRIX *Amat_mat1,
			       double **x, double **b, AZ_MATRIX **Amat_komplex)

\brief Create Komplex System from General System.

   Transforms a complex-valued system 

   (c0r+i*c0i)*A0 +(c1r+i*c1i)*A1) * (xr+i*xi) = (br+i*bi)

   to a Komplex system.

\param xr (In)
       Real part of initial guess.
\param xi (In)
       Imaginary part of initial guess.
\param br (In)
       Real part of right hand side of linear system.
\param bi (In)
       Imaginary part of right hand side of linear system.
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
\param c1r (In)
       Real part of constant to be multiplied with second matrix.
\param c1i (In)
       Imaginary part of constant to be multiplied with second matrix.
\param Amat_mat0 (In)
       AZ_MATRIX object containing first real-valued matrix.
\param Amat_mat1 (In)
       AZ_MATRIX object containing second real-valued matrix.

\param x (Out)
       Komplex version of initial guess and solution.
\param b (Out)
       Komplex version of RHS.
\param Amat_komplex (Out)
       Komplex version of matrix stored as an AZ_MATRIX structure.
*/

void AZK_create_linsys_g2k(double *xr, double *xi, double *br, double *bi, 
			       int *options, double *params, int *proc_config,
			       double c0r, double c0i, AZ_MATRIX *Amat_mat0, 
			       double c1r, double c1i, AZ_MATRIX *Amat_mat1,
			       double **x, double **b, AZ_MATRIX **Amat_komplex)

{
  /* First executable statement */

  /* Create Komplex matrix */

  AZK_create_matrix_g2k( options, params, proc_config, c0r, c0i, Amat_mat0, 
			    c1r, c1i, Amat_mat1, Amat_komplex);

  /* Build RHS and initial guess */

  AZK_create_vector_g2k( options, params, proc_config, 
				  (*Amat_komplex), xr, xi, x);

  AZK_create_vector_g2k( options, params, proc_config, 
				  (*Amat_komplex), br, bi, b);

  /* Permute K system for better numerical stability */

  AZK_permute_ri( options, params, proc_config, (*b), (*Amat_komplex));

  return;
}
/*! 
\fn void AZK_create_linsys_ri2k(double *xr, double *xi, double *br, double *bi, 
			      int *options, double *params, int *proc_config,
			      AZ_MATRIX *Amat_real, double *val_imag,
			      double **x, double **b, AZ_MATRIX **Amat_komplex)

\brief Create Komplex System from Real and Imaginary Parts.

   Transforms a complex-valued system 

         (Ar +i*Ai) * (xr + i*xi) = (br + i*bi)

   where double precision arrays hold the real and imaginary parts separately.
   The pattern of the imaginary part matches the real part. Thus no structure 
   for the imaginary part is passed in.

\param xr (In)
       Real part of initial guess.
\param xi (In)
       Imaginary part of initial guess.
\param br (In)
       Real part of right hand side of linear system.
\param bi (In)
       Imaginary part of right hand side of linear system.
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

\param x (Out)
       Komplex version of initial guess and solution.
\param b (Out)
       Komplex version of RHS.
\param Amat_komplex (Out)
       Komplex version of matrix stored as an AZ_MATRIX structure.

*/

void AZK_create_linsys_ri2k(double *xr, double *xi, double *br, double *bi, 
			      int *options, double *params, int *proc_config,
			      AZ_MATRIX *Amat_real, double *val_imag,
			      double **x, double **b, AZ_MATRIX **Amat_komplex)
{
  /* First executable statement */

  /* Create Komplex matrix */

  AZK_create_matrix_ri2k( options, params, proc_config, Amat_real, 
			  val_imag, Amat_komplex);

  /* Build RHS and initial guess */

  AZK_create_vector_ri2k( options, params, proc_config, 
				  (*Amat_komplex), xr, xi, x);
  AZK_create_vector_ri2k( options, params, proc_config, 
				  (*Amat_komplex), br, bi, b);

  /* Permute K system for better numerical stability */

  AZK_permute_ri( options, params, proc_config, (*b),(*Amat_komplex));

  return;
}
#ifndef DOXYGEN_SHOULD_SKIP_THIS
void AZK_create_linsys_no_copy(double *xr, double *xi, double *br, double *bi, 
			      int *options, double *params, int *proc_config,
			      AZ_MATRIX *Amat_real, AZ_MATRIX *Amat_imag,
			      double **x, double **b, AZ_MATRIX **Amat)
{
/* 
   Transforms a complex-valued system 

         (Ar +i*Ai) * (xr + i*xi) = (br + i*bi)

   where double precision arrays hold the real and imaginary parts separately.

 Input arguments:
 ================
     
 xr,xi:         On input, contains the initial guess, real part in xr and
                imaginary part in xi. On output contains the solution to 
                the linear system.

 br,bi:         Right hand side of linear system.

 Output arguments:
 =================

 x:             Komplex version of initial guess and solution.
 b:             Komplex version of RHS.
 Amat:          Komplex version of matrix stored as an AZ_MATRIX structure.

*/
  AZ_KOMPLEX *linsys_pass_data;
  int N_equations, N_blk_equations, N_real, N_external;
  int *data_org_real, *data_org_imag;
  int *komplex_to_real, *komplex_to_imag;
  int i;


  if (Amat_real->has_global_indices || Amat_imag->has_global_indices)
    AZ_perror("AZK_create_linsys_no_copy requires local indices");

  linsys_pass_data = (AZ_KOMPLEX *) AZ_allocate(sizeof(AZ_KOMPLEX));
  if (linsys_pass_data == NULL)
    AZ_perror("AZK_create_linsys_no_copy: Out of memory.");
  data_org_real = Amat_real->data_org;
  data_org_imag = Amat_imag->data_org;
  N_real = data_org_real[AZ_N_internal] + data_org_real[AZ_N_border];
  
  N_equations = 2 * N_real;
  N_blk_equations = N_equations;

  N_external = AZ_MAX(data_org_real[AZ_N_external], data_org_imag[AZ_N_external]);

  if (Amat_real->matrix_type == AZ_MSR_MATRIX) {
      Amat_real->data_org[AZ_N_int_blk] = Amat_real->data_org[AZ_N_internal];
      Amat_real->data_org[AZ_N_bord_blk] = Amat_real->data_org[AZ_N_border];
      N_blk_equations = N_equations;
  }

  else if (Amat_real->matrix_type == AZ_VBR_MATRIX) 
    {
      N_blk_equations = data_org_real[AZ_N_int_blk] + data_org_real[AZ_N_bord_blk];
    }
  else if (Amat_real->matrix_type == AZ_USER_MATRIX){
      Amat_real->data_org[AZ_N_int_blk] = Amat_real->data_org[AZ_N_internal];
      Amat_real->data_org[AZ_N_bord_blk] = Amat_real->data_org[AZ_N_border];
      N_blk_equations = N_equations;
  }
  else
     AZ_perror("AZK_create_linsys_no_copy: Unknown matrix type.");
 
  if (Amat_imag->matrix_type == AZ_MSR_MATRIX) {
      Amat_imag->data_org[AZ_N_int_blk] = Amat_imag->data_org[AZ_N_internal];
      Amat_imag->data_org[AZ_N_bord_blk] = Amat_imag->data_org[AZ_N_border];
  }

  else if (Amat_imag->matrix_type == AZ_USER_MATRIX){
      Amat_imag->data_org[AZ_N_int_blk] = Amat_imag->data_org[AZ_N_internal];
      Amat_imag->data_org[AZ_N_bord_blk] = Amat_imag->data_org[AZ_N_border];
  }
 
  (*Amat) = AZ_create_matrix(N_equations, 0,
                AZ_USER_MATRIX, N_blk_equations,
                AZ_NOT_USING_AZTEC_MATVEC);

  /* Merge real and imaginary parts into K matrix order */
  komplex_to_real = (int *) AZ_allocate (N_real*sizeof(int));
  komplex_to_imag = (int *) AZ_allocate (N_real*sizeof(int));
  (*x) = (double *) AZ_allocate((N_equations+N_external)*sizeof(double));
  (*b) = (double *) AZ_allocate((N_equations+N_external)*sizeof(double));
  if ((*b) == NULL)
    AZ_perror("AZK_create_linsys_no_copy: Out of memory.");
  for (i=0; i <N_real; i++)
    {
      komplex_to_real[i] = 2*i;
      komplex_to_imag[i] = 2*i+1;
      (*x)[komplex_to_real[i]] = xr[i];
      (*x)[komplex_to_imag[i]] = xi[i];
      (*b)[komplex_to_real[i]] = br[i];
      (*b)[komplex_to_imag[i]] = bi[i];
    }
  linsys_pass_data->Amat_real = Amat_real;
  linsys_pass_data->Amat_imag = Amat_imag;
  linsys_pass_data->komplex_to_real = komplex_to_real;
  linsys_pass_data->komplex_to_imag = komplex_to_imag;

  linsys_pass_data->c11 = 1.0;
  linsys_pass_data->c12 = 0.0;
  linsys_pass_data->c21 = 0.0;
  linsys_pass_data->c22 = 1.0;
  linsys_pass_data->Form_of_Equations = AZK_Komplex_No_Copy;
  linsys_pass_data->From_Global_Indices = 0;

  (*Amat)->matvec  = AZK_matvec_no_copy;
  (*Amat)->aux_ptr = (void *)  linsys_pass_data;

  return;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
