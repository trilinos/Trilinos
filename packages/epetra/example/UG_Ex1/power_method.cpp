//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
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
// ************************************************************************
//@HEADER

#include <stdio.h>
#include <stdlib.h>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
// Simple Power method algorithm
double power_method(const Epetra_CrsMatrix& A) {  
  // variable needed for iteration
  double lambda = 0.0;
  int niters = A.RowMap().NumGlobalElements()*10;
  double tolerance = 1.0e-10;
  // Create vectors
  Epetra_Vector q(A.RowMap());
  Epetra_Vector z(A.RowMap());
  Epetra_Vector resid(A.RowMap());
  // Fill z with random Numbers
  z.Random();
  // variable needed for iteration
  double normz;
  double residual = 0;
  int iter = 0;
  while (iter==0 || (iter < niters && residual > tolerance)) {
    z.Norm2(&normz); // Compute 2-norm of z
    q.Scale(1.0/normz, z);
    A.Multiply(false, q, z); // Compute z = A*q
    q.Dot(z, &lambda); // Approximate maximum eigenvalue
    if (iter%10==0 || iter+1==niters) {
      // Compute A*q - lambda*q every 10 iterations
      resid.Update(1.0, z, -lambda, q, 0.0);
      resid.Norm2(&residual);
      if (q.Map().Comm().MyPID()==0)
	cout << "Iter = " << iter << "  Lambda = " << lambda 
	     << "  Two-norm of A*q - lambda*q = " 
	     << residual << endl;
    } 
    iter++;
  }
  return(lambda);
}
