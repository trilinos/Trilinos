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


#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_JadMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Flops.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "Epetra_Version.h"
#include <vector>
#include <algorithm>
#include <string>

// prototypes

int checkValues( double x, double y, string message = "", bool verbose = false) { 
  if (fabs((x-y)/x) > 0.01) {
    return(1); 
    if (verbose) cout << "********** " << message << " check failed.********** " << endl;
  }
  else {
    if (verbose) cout << message << " check OK." << endl;    
    return(0);
  }
}

int checkMultiVectors( Epetra_MultiVector & X, Epetra_MultiVector & Y, string message = "", bool verbose = false) {
  int numVectors = X.NumVectors();
  int length = Y.MyLength();
  int badvalue = 0;
  int globalbadvalue = 0;
  for (int j=0; j<numVectors; j++)
    for (int i=0; i< length; i++) 
      if (checkValues(X[j][i], Y[j][i])==1) badvalue = 1;
  X.Map().Comm().MaxAll(&badvalue, &globalbadvalue, 1);

  if (verbose) {
    if (globalbadvalue==0) cout << message << " check OK." << endl;
    else cout << "********* " << message << " check failed.********** " << endl;
  }
  return(globalbadvalue);
}

int check(Epetra_RowMatrix & A, Epetra_RowMatrix & B, bool verbose);

int powerMethodTests(Epetra_RowMatrix & A, Epetra_RowMatrix & JadA, Epetra_Map & Map, 
		     Epetra_Vector & q, Epetra_Vector & z, Epetra_Vector & resid, bool verbose);

int power_method(bool TransA, Epetra_RowMatrix& A, 
		 Epetra_Vector& q,
		 Epetra_Vector& z0, 
		 Epetra_Vector& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose);

int main(int argc, char *argv[])
{
  int ierr = 0, i, forierr = 0;
#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;


  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  // Redefine verbose to only print on PE 0
  if(verbose && rank!=0) 
		verbose = false;

  int NumMyEquations = 10000;
  long long NumGlobalEquations = (NumMyEquations * NumProc) + EPETRA_MIN(NumProc,3);
  if(MyPID < 3) 
    NumMyEquations++;

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_Map Map(NumGlobalEquations, NumMyEquations, 0LL, Comm);
  
  // Get update list and number of local equations from newly created Map
  vector<long long> MyGlobalElements(Map.NumMyElements());
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  vector<int> NumNz(NumMyEquations);

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for(i = 0; i < NumMyEquations; i++)
    if((MyGlobalElements[i] == 0) || (MyGlobalElements[i] == NumGlobalEquations - 1))
      NumNz[i] = 1;
    else
      NumNz[i] = 2;

  // Create a Epetra_Matrix

  Epetra_CrsMatrix A(Copy, Map, &NumNz[0]);
  EPETRA_TEST_ERR(A.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(A.IndicesAreLocal(),ierr);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  vector<double> Values(2);
  Values[0] = -1.0; 
	Values[1] = -1.0;
	vector<long long> Indices(2);
  double two = 2.0;
  int NumEntries;

  forierr = 0;
  for(i = 0; i < NumMyEquations; i++) {
    if(MyGlobalElements[i] == 0) {
			Indices[0] = 1;
			NumEntries = 1;
		}
    else if (MyGlobalElements[i] == NumGlobalEquations-1) {
			Indices[0] = NumGlobalEquations-2;
			NumEntries = 1;
		}
    else {
			Indices[0] = MyGlobalElements[i]-1;
			Indices[1] = MyGlobalElements[i]+1;
			NumEntries = 2;
		}
		forierr += !(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0])==0);
		forierr += !(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, &MyGlobalElements[i])>0); // Put in the diagonal entry
  }
  EPETRA_TEST_ERR(forierr,ierr);

  // Finish up
  A.FillComplete();
  A.OptimizeStorage();

  Epetra_JadMatrix JadA(A);
  Epetra_JadMatrix JadA1(A);
  Epetra_JadMatrix JadA2(A);

  // Create vectors for Power method

  Epetra_Vector q(Map);
  Epetra_Vector z(Map); z.Random();
  Epetra_Vector resid(Map);

  Epetra_Flops flopcounter;
  A.SetFlopCounter(flopcounter);
  q.SetFlopCounter(A);
  z.SetFlopCounter(A);
  resid.SetFlopCounter(A);
  JadA.SetFlopCounter(A);
  JadA1.SetFlopCounter(A);
  JadA2.SetFlopCounter(A);
  

  if (verbose) cout << "=======================================" << endl
		    << "Testing Jad using CrsMatrix as input..." << endl
		    << "=======================================" << endl;

  A.ResetFlops();
  powerMethodTests(A, JadA, Map, q, z, resid, verbose);

  // Increase diagonal dominance

  if (verbose) cout << "\n\nIncreasing the magnitude of first diagonal term and solving again\n\n"
		    << endl;

  
  if (A.MyGlobalRow(0)) {
    int numvals = A.NumGlobalEntries(0);
    vector<double> Rowvals(numvals);
    vector<long long> Rowinds(numvals);
    A.ExtractGlobalRowCopy(0, numvals, numvals, &Rowvals[0], &Rowinds[0]); // Get A[0,0]

    for (i=0; i<numvals; i++) if (Rowinds[i] == 0) Rowvals[i] *= 10.0;
    
    A.ReplaceGlobalValues(0, numvals, &Rowvals[0], &Rowinds[0]);
  }
  JadA.UpdateValues(A);
  A.ResetFlops();
  powerMethodTests(A, JadA, Map, q, z, resid, verbose);

  if (verbose) cout << "================================================================" << endl
		          << "Testing Jad using Jad matrix as input matrix for construction..." << endl
		          << "================================================================" << endl;
  JadA1.ResetFlops();
  powerMethodTests(JadA1, JadA2, Map, q, z, resid, verbose);

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return ierr ;
}

int powerMethodTests(Epetra_RowMatrix & A, Epetra_RowMatrix & JadA, Epetra_Map & Map, 
		     Epetra_Vector & q, Epetra_Vector & z, Epetra_Vector & resid, bool verbose) {

  // variable needed for iteration
  double lambda = 0.0;
  // int niters = 10000;
  int niters = 300;
  double tolerance = 1.0e-2;
  int ierr = 0;

  /////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Iterate

  Epetra_Time timer(Map.Comm());
	
  double startTime = timer.ElapsedTime();
  EPETRA_TEST_ERR(power_method(false, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  double elapsed_time = timer.ElapsedTime() - startTime;
  double total_flops = q.Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;
  double lambdaref = lambda;
  double flopsref = total_flops;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for reference first solve = " << MFLOPs << endl
		  <<     "Total FLOPS                            = " <<total_flops <<endl<<endl;

  lambda = 0.0;
  startTime = timer.ElapsedTime();
  EPETRA_TEST_ERR(power_method(false, JadA, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = timer.ElapsedTime() - startTime;
  total_flops = q.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for candidate first solve = " << MFLOPs << endl
		  <<     "Total FLOPS                            = " <<total_flops <<endl<<endl;

  EPETRA_TEST_ERR(checkValues(lambda,lambdaref," No-transpose Power Method result", verbose),ierr);
  EPETRA_TEST_ERR(checkValues(total_flops,flopsref," No-transpose Power Method flop count", verbose),ierr);

  /////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << endl;
  // Iterate
  lambda = 0.0;
  startTime = timer.ElapsedTime();
  EPETRA_TEST_ERR(power_method(true, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = timer.ElapsedTime() - startTime;
  total_flops = q.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  lambdaref = lambda;
  flopsref = total_flops;

  if (verbose) 
	 cout << "\n\nTotal MFLOPs for reference transpose solve = " << MFLOPs << endl
		 <<     "Total FLOPS                                = " <<total_flops <<endl<<endl;

  lambda = 0.0;
  startTime = timer.ElapsedTime();
  EPETRA_TEST_ERR(power_method(true, JadA, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = timer.ElapsedTime() - startTime;
  total_flops = q.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) 
	  cout << "\n\nTotal MFLOPs for candidate transpose solve = " << MFLOPs << endl
		  <<     "Total FLOPS                                = " <<total_flops <<endl<<endl;

  EPETRA_TEST_ERR(checkValues(lambda,lambdaref,"Transpose Power Method result", verbose),ierr);
  EPETRA_TEST_ERR(checkValues(total_flops,flopsref,"Transpose Power Method flop count", verbose),ierr);

  EPETRA_TEST_ERR(check(A, JadA, verbose),ierr);

  return(0);
}
int power_method(bool TransA, Epetra_RowMatrix& A, Epetra_Vector& q, Epetra_Vector& z0, 
		 Epetra_Vector& resid, double* lambda, int niters, double tolerance, bool verbose) 
{  
	
  // Fill z with random Numbers
  Epetra_Vector z(z0);
	
  // variable needed for iteration
  double normz, residual;

  int ierr = 1;
	
  for(int iter = 0; iter < niters; iter++) {
		z.Norm2(&normz); // Compute 2-norm of z
		q.Scale(1.0/normz, z);
		A.Multiply(TransA, q, z); // Compute z = A*q // SEGFAULT HAPPENS HERE
		q.Dot(z, lambda); // Approximate maximum eigenvaluE
		if(iter%100==0 || iter+1==niters) {
			resid.Update(1.0, z, -(*lambda), q, 0.0); // Compute A*q - lambda*q
			resid.Norm2(&residual);
			if(verbose) cout << "Iter = " << iter << "  Lambda = " << *lambda 
											 << "  Residual of A*q - lambda*q = " << residual << endl;
		}
		if(residual < tolerance) {
			ierr = 0;
			break;
		}
	}
  return(ierr);
}

int check(Epetra_RowMatrix& A, Epetra_RowMatrix & B, bool verbose)  {  

  int ierr = 0;
  EPETRA_TEST_ERR(!A.Comm().NumProc()==B.Comm().NumProc(),ierr);
  EPETRA_TEST_ERR(!A.Comm().MyPID()==B.Comm().MyPID(),ierr);
  EPETRA_TEST_ERR(!A.Filled()==B.Filled(),ierr);
  EPETRA_TEST_ERR(!A.HasNormInf()==B.HasNormInf(),ierr);
  EPETRA_TEST_ERR(!A.LowerTriangular()==B.LowerTriangular(),ierr);
  EPETRA_TEST_ERR(!A.Map().SameAs(B.Map()),ierr);
  EPETRA_TEST_ERR(!A.MaxNumEntries()==B.MaxNumEntries(),ierr);
  EPETRA_TEST_ERR(!A.NumGlobalCols64()==B.NumGlobalCols64(),ierr);
  EPETRA_TEST_ERR(!A.NumGlobalDiagonals64()==B.NumGlobalDiagonals64(),ierr);
  EPETRA_TEST_ERR(!A.NumGlobalNonzeros64()==B.NumGlobalNonzeros64(),ierr);
  EPETRA_TEST_ERR(!A.NumGlobalRows64()==B.NumGlobalRows64(),ierr);
  EPETRA_TEST_ERR(!A.NumMyCols()==B.NumMyCols(),ierr);
  EPETRA_TEST_ERR(!A.NumMyDiagonals()==B.NumMyDiagonals(),ierr);
  EPETRA_TEST_ERR(!A.NumMyNonzeros()==B.NumMyNonzeros(),ierr);
  for (int i=0; i<A.NumMyRows(); i++) {
    int nA, nB;
    A.NumMyRowEntries(i,nA); B.NumMyRowEntries(i,nB);
    EPETRA_TEST_ERR(!nA==nB,ierr);
  }
  EPETRA_TEST_ERR(!A.NumMyRows()==B.NumMyRows(),ierr);
  EPETRA_TEST_ERR(!A.OperatorDomainMap().SameAs(B.OperatorDomainMap()),ierr);
  EPETRA_TEST_ERR(!A.OperatorRangeMap().SameAs(B.OperatorRangeMap()),ierr);
  EPETRA_TEST_ERR(!A.RowMatrixColMap().SameAs(B.RowMatrixColMap()),ierr);
  EPETRA_TEST_ERR(!A.RowMatrixRowMap().SameAs(B.RowMatrixRowMap()),ierr);
  EPETRA_TEST_ERR(!A.UpperTriangular()==B.UpperTriangular(),ierr);
  EPETRA_TEST_ERR(!A.UseTranspose()==B.UseTranspose(),ierr);

  int NumVectors = 5;
  { // No transpose case
    Epetra_MultiVector X(A.OperatorDomainMap(), NumVectors);
    Epetra_MultiVector YA1(A.OperatorRangeMap(), NumVectors);
    Epetra_MultiVector YA2(YA1);
    Epetra_MultiVector YB1(YA1);
    Epetra_MultiVector YB2(YA1);
    X.Random();
    
    bool transA = false;
    A.SetUseTranspose(transA);
    B.SetUseTranspose(transA);
    A.Apply(X,YA1);
    A.Multiply(transA, X, YA2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YA2,"A Multiply and A Apply", verbose),ierr);
    B.Apply(X,YB1);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB1,"A Multiply and B Multiply", verbose),ierr);
    B.Multiply(transA, X, YB2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB2,"A Multiply and B Apply", verbose), ierr);
    
  }
  {// transpose case
    Epetra_MultiVector X(A.OperatorRangeMap(), NumVectors);
    Epetra_MultiVector YA1(A.OperatorDomainMap(), NumVectors);
    Epetra_MultiVector YA2(YA1);
    Epetra_MultiVector YB1(YA1);
    Epetra_MultiVector YB2(YA1);
    X.Random();
    
    bool transA = true;
    A.SetUseTranspose(transA);
    B.SetUseTranspose(transA);
    A.Apply(X,YA1);
    A.Multiply(transA, X, YA2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YA2, "A Multiply and A Apply (transpose)", verbose),ierr);
    B.Apply(X,YB1);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB1, "A Multiply and B Multiply (transpose)", verbose),ierr);
    B.Multiply(transA, X,YB2);
    EPETRA_TEST_ERR(checkMultiVectors(YA1,YB2, "A Multiply and B Apply (transpose)", verbose),ierr);
    
  }

  Epetra_Vector diagA(A.RowMatrixRowMap());
  EPETRA_TEST_ERR(A.ExtractDiagonalCopy(diagA),ierr);
  Epetra_Vector diagB(B.RowMatrixRowMap());
  EPETRA_TEST_ERR(B.ExtractDiagonalCopy(diagB),ierr);
  EPETRA_TEST_ERR(checkMultiVectors(diagA,diagB, "ExtractDiagonalCopy", verbose),ierr);

  Epetra_Vector rowA(A.RowMatrixRowMap());
  EPETRA_TEST_ERR(A.InvRowSums(rowA),ierr);
  Epetra_Vector rowB(B.RowMatrixRowMap());
  EPETRA_TEST_ERR(B.InvRowSums(rowB),ierr)
  EPETRA_TEST_ERR(checkMultiVectors(rowA,rowB, "InvRowSums", verbose),ierr);

  Epetra_Vector colA(A.RowMatrixColMap());
  EPETRA_TEST_ERR(A.InvColSums(colA),ierr);
  Epetra_Vector colB(B.RowMatrixColMap());
  EPETRA_TEST_ERR(B.InvColSums(colB),ierr);
  EPETRA_TEST_ERR(checkMultiVectors(colA,colB, "InvColSums", verbose),ierr);

  EPETRA_TEST_ERR(checkValues(A.NormInf(), B.NormInf(), "NormInf before scaling", verbose), ierr);
  EPETRA_TEST_ERR(checkValues(A.NormOne(), B.NormOne(), "NormOne before scaling", verbose),ierr);

  EPETRA_TEST_ERR(A.RightScale(colA),ierr);  
  EPETRA_TEST_ERR(B.RightScale(colB),ierr);


  EPETRA_TEST_ERR(A.LeftScale(rowA),ierr);
  EPETRA_TEST_ERR(B.LeftScale(rowB),ierr);


  EPETRA_TEST_ERR(checkValues(A.NormInf(), B.NormInf(), "NormInf after scaling", verbose), ierr);
  EPETRA_TEST_ERR(checkValues(A.NormOne(), B.NormOne(), "NormOne after scaling", verbose),ierr);

  vector<double> valuesA(A.MaxNumEntries());
  vector<int> indicesA(A.MaxNumEntries());  
  vector<double> valuesB(B.MaxNumEntries());
  vector<int> indicesB(B.MaxNumEntries());
  return(0);
  for (int i=0; i<A.NumMyRows(); i++) {
    int nA, nB;
    EPETRA_TEST_ERR(A.ExtractMyRowCopy(i, A.MaxNumEntries(), nA, &valuesA[0], &indicesA[0]),ierr); 
    EPETRA_TEST_ERR(B.ExtractMyRowCopy(i, B.MaxNumEntries(), nB, &valuesB[0], &indicesB[0]),ierr);
    EPETRA_TEST_ERR(!nA==nB,ierr);
    for (int j=0; j<nA; j++) {
      double curVal = valuesA[j];
      int curIndex = indicesA[j];
      bool notfound = true;
      int jj = 0;
      while (notfound && jj< nB) {
	if (!checkValues(curVal, valuesB[jj])) notfound = false;
	jj++;
      }
      EPETRA_TEST_ERR(notfound, ierr);
      vector<int>::iterator p = find(indicesB.begin(),indicesB.end(),curIndex);  // find curIndex in indicesB
      EPETRA_TEST_ERR(p==indicesB.end(), ierr);
    }

  }
  if (verbose) cout << "RowMatrix Methods check OK" << endl;

  return (ierr);
}
