/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "azk_komplex.h"

/*! \file
\brief Permutation routine that checks real and imaginary parts and
swaps if needed for better numerical stability.

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

/*! \fn void AZK_permute_ri(int *options, double *params,
			    int *proc_config, double *b,
			    AZ_MATRIX *Amat_komplex)

\brief Permute a Komplex system for better numerical stability.

An alternative to the standard Komplex formulation is to permute the
block rows so that the imaginary part is on the main diagonal.  For
example:

\f[
  \left( \begin{array}{rr}
                                    A_i &  A_r\\
                                    A_r & -A_i
                             \end{array}
   \right)
   \left( \begin{array}{r}
                                    x_r\\
                                    x_i
                             \end{array}
   \right)
   =
   \left( \begin{array}{r}
                                    b_i\\
                                    b_r
                             \end{array}
   \right)
\f]

This
action may be desirable, or necessary in situations where the real
part has small or zero diagonal entries.  This routine looks at each
real/imaginary pair and, based on a heuristic may swap the real and
imaginary parts.  This action does not affect the sparsity pattern,
but only the mapping from the complex (or real/imaginary) mapping to
the komplex mapping, and back.

\param options (In)
       Determines specific solution method and other parameters.
\param params (In)
       Drop tolerance and convergence tolerance info.
\param proc_config (In)
       Machine configuration.  proc_config[AZ_node] is the node
       number.  proc_config[AZ_N_procs] is the number of processors.

\param b (Out)
       Komplex version of RHS, possibly permuted.
\param Amat_komplex (Out)
       Komplex version of matrix stored as an AZ_MATRIX structure,
       possibly permuted.

*/

void AZK_permute_ri(int *options, double *params,
			    int *proc_config, double *b,
			    AZ_MATRIX *Amat_komplex)
{
  /*
  int N_equations, N_real;
  int *data_org;
  AZ_KOMPLEX *linsys_pass_data;
  */

  /* First executable statement */

  /*
  data_org = Amat_komplex->data_org;

  N_equations = data_org[AZ_N_internal] + data_org[AZ_N_border];
  N_real = N_equations/2;

  linsys_pass_data = (AZ_KOMPLEX *) Amat_komplex->aux_ptr;
  */
  return;
}


