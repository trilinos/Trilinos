/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include "ifp_Matrix.h"
#include "ifp_Precon.h"
#include "ifp_ifpack.h"

double ifp_Precon::condest()
{
  // This routine computes a  bound of the infinity-norm condition number 
  // of the preconditioner.
  // It is equal to infinity-norm of M^-1 e, where e is the vector of all
  // ones.

  int i;

  int m = dimrow(); // Object is a matrix, so it has row/col dimensions.
  int n = dimcol();
  double *u = new double[n];
  double *v = new double[m];

    for (i=0; i<n; i++) u[i] = 1.0;

    apply(n, 1, u, n, v, n);

    double cond_number = 0.0;
    for (i=0; i<m; i++)
        if (ABS(v[i]) > cond_number)
           cond_number = ABS(v[i]);

    delete [] u;
    delete [] v;

    return(cond_number);
}
