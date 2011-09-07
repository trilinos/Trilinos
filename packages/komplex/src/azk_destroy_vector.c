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
\brief Destruction routine for deleting Komplex vectors.

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


/*! \fn void AZK_destroy_vector(int *options, double *params,
		       int *proc_config, AZ_MATRIX *Amat_komplex,
                       double **vk)
\brief Destroy a Komplex vector.

Destroys a komplex vector created by any of the AZK_create_vector
functions.  Deletes any memory allocated by creation routine.

\param options (In)
       Determines specific solution method and other parameters.
\param params (In)
       Drop tolerance and convergence tolerance info.
\param proc_config (In)
       Machine configuration.  proc_config[AZ_node] is the node
       number.  proc_config[AZ_N_procs] is the number of processors.
\param Amat_komplex (In)
       Komplex version of matrix stored as an AZ_MATRIX structure.

\param vk (Out)
       Deleted komplex version of a vector.  Remember to call
       AZK_extract_solution_[k2c,g2k,ri2k] before calling this
       routine.

*/


void AZK_destroy_vector(int *options, double *params,
		       int *proc_config, AZ_MATRIX *Amat_komplex,
                       double **vk)
{

  /* First executable statement */

  AZ_free((void *) (*vk));
  *vk = NULL;
  
  return;
}
