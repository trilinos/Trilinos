//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
