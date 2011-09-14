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
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"

// *************************************************************
// Function prototypes
// *************************************************************

int invIteration(Epetra_CrsMatrix& A, double &lambda, bool verbose); 
int applyInverseSetup(Epetra_CrsMatrix &A, Ifpack_CrsRiluk * & M);
int applyInverse(Epetra_CrsMatrix &A, Epetra_Vector & x, Epetra_Vector & b, Ifpack_CrsRiluk * M, 
		 bool verbose);
int applyInverseDestroy(Ifpack_CrsRiluk * M);
void BiCGSTAB(Epetra_CrsMatrix &A, 
	      Epetra_Vector &x, 
	      Epetra_Vector &b, 
	      Ifpack_CrsRiluk *M, 
	      int Maxiter, 
	      double Tolerance, bool verbose);

// *************************************************************
// main program - This benchmark code reads a Harwell-Boeing data
//                set and finds the minimal eigenvalue of the matrix
//                using inverse iteration.
// *************************************************************
int main(int argc, char *argv[]) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  cout << Comm << endl;

  int MyPID = Comm.MyPID();

  bool verbose = false; 
  if (MyPID==0) verbose = true; // Print out detailed results (turn off for best performance)

  if(argc != 2) {
    if (verbose) cerr << "Usage: " << argv[0] << " HB_data_file" << endl;
    exit(1); // Error
  }

  // Define pointers that will be set by HB read function

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
   
  // Call function to read in HB problem
  Trilinos_Util_ReadHb2Epetra(argv[1], Comm, readMap, readA, readx, readb, readxexact);

  // Not interested in x, b or xexact for an eigenvalue problem
  delete readx;
  delete readb;
  delete readxexact;

#ifdef EPETRA_MPI // If running in parallel, we need to distribute matrix across all PEs.

  // Create uniform distributed map
  Epetra_Map map(readMap->NumGlobalElements(), 0, Comm);

  // Create Exporter to distribute read-in matrix and vectors

  Epetra_Export exporter(*readMap, map);
  Epetra_CrsMatrix A(Copy, map, 0);

  A.Export(*readA, exporter, Add);
  assert(A.FillComplete()==0);

  delete readA;
  delete readMap;

#else // If not running in parallel, we do not need to distribute the matrix
  Epetra_CrsMatrix & A = *readA;
#endif

  // Create flop counter to collect all FLOPS
  Epetra_Flops counter;
  A.SetFlopCounter(counter);

  double lambda = 0; // Minimal eigenvalue returned here
  // Call inverse iteration solver
  Epetra_Time timer(Comm);
  invIteration(A, lambda, verbose);
  double elapsedTime = timer.ElapsedTime();
  double totalFlops = counter.Flops();
  double MFLOPS = totalFlops/elapsedTime/1000000.0;


  cout << endl
       << "*************************************************" << endl
       << " Approximate smallest eigenvalue = " << lambda << endl
       << "    Total Time    = " << elapsedTime << endl
       << "    Total FLOPS   = " << totalFlops << endl
       << "    Total MFLOPS  = " << MFLOPS << endl
       << "*************************************************" << endl;

  // All done
#ifdef EPETRA_MPI
  MPI_Finalize();
#else
  delete readA;
  delete readMap;  
#endif

return (0);
}

// The following functions are used to solve the problem:


// *************************************************************
// Computes the smallest eigenvalue of the given matrix A.
// Returns result in lambda
// *************************************************************

int invIteration(Epetra_CrsMatrix& A, double &lambda, bool verbose) {  

  Ifpack_CrsRiluk * M;
  applyInverseSetup(A, M);
  Epetra_Vector q(A.RowMap());
  Epetra_Vector z(A.RowMap());
  Epetra_Vector resid(A.RowMap());

  Epetra_Flops * counter = A.GetFlopCounter();
  if (counter!=0) {
    q.SetFlopCounter(A);
    z.SetFlopCounter(A);
    resid.SetFlopCounter(A);
  }

  // Fill z with random Numbers
  z.Random();

  // variable needed for iteration
  double normz, residual;

  int niters = 100;
  double tolerance = 1.0E-6;
  int ierr = 1;

  for (int iter = 0; iter < niters; iter++)
    {
      if (verbose) 
	cout << endl
	     << " ***** Performing step " << iter << " of inverse iteration ***** " << endl;

      z.Norm2(&normz); // Compute 2-norm of z
      q.Scale(1.0/normz, z);
      applyInverse(A, z, q, M, verbose); // Compute z such that Az = q
      q.Dot(z, &lambda); // Approximate maximum eigenvalue
      if (iter%10==0 || iter+1==niters)
	{
	  resid.Update(1.0, z, -lambda, q, 0.0); // Compute A(inv)*q - lambda*q
	  resid.Norm2(&residual);
	  cout << endl
	       << "***** Inverse Iteration Step " << iter+1 << endl
	       << "  Lambda = " << 1.0/lambda << endl
	       << "  Residual of A(inv)*q - lambda*q = " 
	       << residual << endl;
	} 
      if (residual < tolerance) {
	ierr = 0;
	break;
      }
    }
  // lambda is the largest eigenvalue of A(inv).  1/lambda is smallest eigenvalue of A.
  lambda = 1.0/lambda;

  // Compute A*q - lambda*q explicitly
  A.Multiply(false, q, z);
  resid.Update(1.0, z, -lambda, q, 0.0); // Compute A*q - lambda*q
  resid.Norm2(&residual);
  cout << "  Explicitly computed residual of A*q - lambda*q = " << residual << endl;
  applyInverseDestroy(M);
  return(ierr);
}

// Performs any setup that can be done once to reduce the cost of the applyInverse function
int applyInverseSetup(Epetra_CrsMatrix &A, Ifpack_CrsRiluk * & M) {
  int LevelFill = 4;
  int Overlap = 0;
  Ifpack_IlukGraph * IlukGraph = new Ifpack_IlukGraph(A.Graph(), LevelFill, Overlap);
    assert(IlukGraph->ConstructFilledGraph()==0);
    M = new Ifpack_CrsRiluk(A, *IlukGraph);
    M->SetFlopCounter(A);
    assert(M->InitValues()==0);
    assert(M->Factor()==0);
  return(0);
}

// *************************************************************
// Solves a problem Ax = b, for a given A and b.  
// M is a preconditioner computed in applyInverseSetup.
// *************************************************************

int applyInverse(Epetra_CrsMatrix &A, 
		 Epetra_Vector & x, Epetra_Vector & b, Ifpack_CrsRiluk * M, bool verbose) {

  int Maxiter = 1000;
  double Tolerance = 1.0E-16;
  BiCGSTAB(A, x, b, M, Maxiter, Tolerance, verbose);

  return(0);
}


// *************************************************************
// Releases any memory associated with the preconditioner M
// *************************************************************

int applyInverseDestroy(Ifpack_CrsRiluk * M) {
  Ifpack_IlukGraph & IlukGraph = (Ifpack_IlukGraph &) M->Graph();
  delete M;
  delete &IlukGraph;
  return(0);
}

// *************************************************************
// Uses the Bi-CGSTAB iterative method to solve Ax = b
// *************************************************************

void BiCGSTAB(Epetra_CrsMatrix &A, 
	      Epetra_Vector &x, 
	      Epetra_Vector &b, 
	      Ifpack_CrsRiluk *M, 
	      int Maxiter, 
	      double Tolerance, bool verbose) {

  // Allocate vectors needed for iterations
  Epetra_Vector phat(x.Map()); phat.SetFlopCounter(x);
  Epetra_Vector p(x.Map()); p.SetFlopCounter(x);
  Epetra_Vector shat(x.Map()); shat.SetFlopCounter(x);
  Epetra_Vector s(x.Map()); s.SetFlopCounter(x);
  Epetra_Vector r(x.Map()); r.SetFlopCounter(x);
  Epetra_Vector rtilde(x.Map()); rtilde.Random(); rtilde.SetFlopCounter(x);
  Epetra_Vector v(x.Map()); v.SetFlopCounter(x);
  

  A.Multiply(false, x, r); // r = A*x

  r.Update(1.0, b, -1.0); // r = b - A*x

  double r_norm, b_norm, scaled_r_norm, rhon, rhonm1 = 1.0;
  double alpha = 1.0, omega = 1.0, sigma;
  double omega_num, omega_den;
  r.Norm2(&r_norm);
  b.Norm2(&b_norm);
  scaled_r_norm = r_norm/b_norm;
  r.Dot(rtilde,&rhon);

  if (verbose) cout << "Initial residual = " << r_norm
		    << " Scaled residual = " << scaled_r_norm << endl;


  for (int i=0; i<Maxiter; i++) { // Main iteration loop   

    double beta = (rhon/rhonm1) * (alpha/omega);
    rhonm1 = rhon;

    /* p    = r + beta*(p - omega*v)       */
    /* phat = M^-1 p                       */
    /* v    = A phat                       */

    double dtemp = - beta*omega;

    p.Update(1.0, r, dtemp, v, beta);
    if (M==0) 
      phat.Scale(1.0, p);
    else
      M->Solve(false, p, phat);
    A.Multiply(false, phat, v);

    
    rtilde.Dot(v,&sigma);
    alpha = rhon/sigma;    

    /* s = r - alpha*v                     */
    /* shat = M^-1 s                       */
    /* r = A shat (r is a tmp here for t ) */

    s.Update(-alpha, v, 1.0, r, 0.0);
    if (M==0) 
      shat.Scale(1.0, s);
    else
      M->Solve(false, s, shat);
    A.Multiply(false, shat, r);

    r.Dot(s, &omega_num);
    r.Dot(r, &omega_den);
    omega = omega_num / omega_den;

    /* x = x + alpha*phat + omega*shat */
    /* r = s - omega*r */

    x.Update(alpha, phat, omega, shat, 1.0);
    r.Update(1.0, s, -omega); 
    
    r.Norm2(&r_norm);
    scaled_r_norm = r_norm/b_norm;
    r.Dot(rtilde,&rhon);

    if (verbose) cout << "Iter "<< i << " residual = " << r_norm
		      << " Scaled residual = " << scaled_r_norm << endl;

    if (scaled_r_norm < Tolerance) break;
  }
  return;
}
