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

/*! \file
\brief Include file for Aztec Komplex Library.

The file, along with az_aztec.h must be included in every file that
uses the Komplex library.

KOMPLEX is an add-on module to AZTEC that allows users to solve
complex-valued
linear systems.

KOMPLEX solves a complex-valued linear system Ax = b by solving
an equivalent real-valued system of twice the dimension.
Specifically,
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
  which is a real-valued system of twice the size.  If we find xr and
xi, we
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
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

#ifndef TRILINOS_NO_CONFIG_H

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 * KL 11/25/02
 */
#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include <Komplex_config.h>
#ifdef HAVE_MPI

#ifndef AZTEC_MPI
#define AZTEC_MPI
#endif

#endif
#endif

#include "az_aztec.h"

/* Define some constants for the user */
#define AZK_True_Complex 	     1
#define AZK_Komplex_Same_Structure   2
#define AZK_Komplex_General          3
#define AZK_Komplex_No_Copy 	     4
struct AZ_KOMPLEX_STRUCT {                  
/* Data passing structure. This user */
/* defined data structure is used to pass information through   */
/* Aztec and back into the user's matrix-vector product. */
  int Form_of_Equations;
  int From_Global_Indices;
  AZ_MATRIX *Amat_real;
  AZ_MATRIX *Amat_imag;
  AZ_MATRIX *Amat_complex;
  int *komplex_to_real;
  int *komplex_to_imag;
  int *komplex_to_complex;
  int AZK_precond;
  int *external, *update_index, *extern_index;

  double c11, c12, c21, c22;
};

typedef struct AZ_KOMPLEX_STRUCT AZ_KOMPLEX;

void AZK_create_linsys_c2k(double *xc, double *bc,
                      int *options, double *params, int *proc_config,
                      AZ_MATRIX *Amat_complex,
                      double **x, double **b, AZ_MATRIX **Amat_komplex);

void AZK_create_linsys_g2k(double *xr, double *xi, double *br, double *bi,
                      int *options, double *params, int *proc_config,
                      double c0r, double c0i, AZ_MATRIX *Amat_mat0,
                      double c1r, double c1i, AZ_MATRIX *Amat_mat1,
                      double **x, double **b, AZ_MATRIX **Amat_komplex);

void AZK_create_linsys_ri2k(double *xr, double *xi, double *br, double *bi,
                     int *options, double *params, int *proc_config,
                     AZ_MATRIX *Amat_real, double *val_imag,
                     double **x, double **b, AZ_MATRIX **Amat_komplex);

void AZK_create_linsys_no_copy(double *xr, double *xi, double *br, double *bi,
                     int *options, double *params, int *proc_config,
                     AZ_MATRIX *Amat_real, AZ_MATRIX *Amat_imag,
                     double **x, double **b, AZ_MATRIX **Amat);

void AZK_create_matrix_g2k(int options[], double params[], int proc_config[], 
		  double c0r, double c0i, AZ_MATRIX *Amat_mat0, 
		  double c1r, double c1i, AZ_MATRIX *Amat_mat1, 
		  AZ_MATRIX **Amat_komplex);

void AZK_create_matrix_g2k_fill_entry(int nrow, int ncol,
                  double c11, double c12, double *mat0v,
                  double c21, double c22, double *mat1v,
                  double *komplex);

void AZK_create_matrix_c2k(int options[], double params[], int proc_config[],
                  AZ_MATRIX *Amat_complex, AZ_MATRIX **Amat_komplex);

void AZK_create_matrix_c2k_fill_entry(int nrow, int ncol,
                  double *cur_complex, double *cur_komplex);

void AZK_create_matrix_ri2k(int options[], double params[], int proc_config[],
                   AZ_MATRIX *Amat_real, double *val_imag,
                   AZ_MATRIX **Amat_komplex);

void AZK_create_matrix_ri2k_fill_entry(int nrow, int ncol,
                  double *realv, double *imagv, double *komplex);

void AZK_create_precon(int *options, double *params,
                 int *proc_config,double *x, double *b,
                 AZ_MATRIX *Amat, AZ_PRECOND **Prec);

void AZK_create_vector_c2k(int *options, double *params,
                 int *proc_config, AZ_MATRIX *Amat_komplex,
                       double *vc, double **vk);

void AZK_create_vector_g2k(int *options, double *params,
                 int *proc_config, AZ_MATRIX *Amat_komplex,
                       double *vr, double *vi, double **vk);

void AZK_create_vector_ri2k(int *options, double *params,
                 int *proc_config, AZ_MATRIX *Amat_komplex,
                       double *vr, double *vi, double **vk);

void AZK_destroy_linsys( int *options, double *params,
                 int *proc_config, 
                 double **x, double **b,
                 AZ_MATRIX **Amat_komplex);

void AZK_destroy_matrix(int options[], double params[], int proc_config[],
            AZ_MATRIX **Amat_komplex);

void AZK_destroy_precon(int *options, double *params, int *proc_config,
                 AZ_MATRIX *Amat, AZ_PRECOND **Prec);

void AZK_destroy_vector(int *options, double *params,
                 int *proc_config, AZ_MATRIX *Amat_komplex,
                       double **vk);

void AZK_extract_solution_k2c(int *options, double *params, int *proc_config,
                    AZ_MATRIX *Amat_komplex, AZ_PRECOND *Prec, double *vk,
                    double *vc);

void AZK_extract_solution_k2g(int *options, double *params, int *proc_config,
                    AZ_MATRIX *Amat_komplex, AZ_PRECOND *Prec, double *vk,
                    double *vr, double *vi );

void AZK_extract_solution_k2ri(int *options, double *params, int *proc_config,
                    AZ_MATRIX *Amat_komplex, AZ_PRECOND *Prec, double *vk,
                    double *vr, double *vi );

void AZK_iterate(double *xr, double *xi, double *br, double *bi,
               int *options, double *params,
               double *status, int *proc_config,
               AZ_MATRIX *Amat_real, AZ_MATRIX *Amat_imag);

void AZK_matvec_no_copy(double *x, double *y, AZ_MATRIX *Amat, 
               int proc_config[]);

void AZK_permute_ri(int *options, double *params,
                   int *proc_config, double *b,
                   AZ_MATRIX *Amat_komplex);

void AZK_precon(double x[], int options[],  
     int proc_config[], double params[], AZ_MATRIX *Amat,
        AZ_PRECOND *Prec);

double AZK_residual_norm_no_copy(double *xr, double *xi, double *br, double *bi, 
                    int *options, double *params,
                    int *proc_config,
                    AZ_MATRIX *Amat_real, AZ_MATRIX *Amat_imag);

double AZK_residual_norm(double *xk, double *bk,
                int *options, double *params,
                int *proc_config,
                AZ_MATRIX *Amat_komplex);


