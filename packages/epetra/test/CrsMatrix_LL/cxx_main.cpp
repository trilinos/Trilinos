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

// prototypes

int check(Epetra_CrsMatrix& A, int NumMyRows1, long long NumGlobalRows1, int NumMyNonzeros1,
	  long long NumGlobalNonzeros1, long long * MyGlobalElements, bool verbose);

int power_method(bool TransA, Epetra_CrsMatrix& A, 
		 Epetra_Vector& q,
		 Epetra_Vector& z, 
		 Epetra_Vector& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose);

int check_graph_sharing(Epetra_Comm& Comm);

int main(int argc, char *argv[])
{
  int ierr = 0, forierr = 0;
  bool debug = false;

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
  //  if (rank==0) cout << "Press any key to continue..."<< std::endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << std::endl << std::endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  bool verbose1 = verbose;

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
  long long* MyGlobalElements = new long long[Map.NumMyElements()];
  Map.MyGlobalElements(MyGlobalElements);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  int* NumNz = new int[NumMyEquations];

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (int i = 0; i < NumMyEquations; i++)
    if((MyGlobalElements[i] == 0) || (MyGlobalElements[i] == NumGlobalEquations - 1))
      NumNz[i] = 1;
    else
      NumNz[i] = 2;

  // Create a Epetra_Matrix

  Epetra_CrsMatrix A(Copy, Map, NumNz);
  EPETRA_TEST_ERR(A.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(A.IndicesAreLocal(),ierr);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double* Values = new double[2];
  Values[0] = -1.0; 
  Values[1] = -1.0;
  long long* Indices = new long long[2];
  double two = 2.0;
  int NumEntries;

  forierr = 0;
  for (int i = 0; i < NumMyEquations; i++) {
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
		forierr += !(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
		forierr += !(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i)>0); // Put in the diagonal entry
  }
  EPETRA_TEST_ERR(forierr,ierr);

  int * indexOffsetTmp;
  int * indicesTmp;
  double * valuesTmp;
  // Finish up
  EPETRA_TEST_ERR(!(A.IndicesAreGlobal()),ierr);
  EPETRA_TEST_ERR(!(A.ExtractCrsDataPointers(indexOffsetTmp, indicesTmp, valuesTmp)==-1),ierr);  // Should fail
  EPETRA_TEST_ERR(!(A.FillComplete(false)==0),ierr);
  EPETRA_TEST_ERR(!(A.ExtractCrsDataPointers(indexOffsetTmp, indicesTmp, valuesTmp)==-1),ierr);  // Should fail
  EPETRA_TEST_ERR(!(A.IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR(A.StorageOptimized(),ierr);
  A.OptimizeStorage();
  EPETRA_TEST_ERR(!(A.StorageOptimized()),ierr);
  EPETRA_TEST_ERR(!(A.ExtractCrsDataPointers(indexOffsetTmp, indicesTmp, valuesTmp)==0),ierr);  // Should succeed
  const Epetra_CrsGraph & GofA = A.Graph();
  EPETRA_TEST_ERR((indicesTmp!=GofA[0] || valuesTmp!=A[0]),ierr); // Extra check to see if operator[] is consistent
  EPETRA_TEST_ERR(A.UpperTriangular(),ierr);
  EPETRA_TEST_ERR(A.LowerTriangular(),ierr);
	
  int NumMyNonzeros = 3 * NumMyEquations;
  if(A.LRID(0) >= 0) 
		NumMyNonzeros--; // If I own first global row, then there is one less nonzero
  if(A.LRID(NumGlobalEquations-1) >= 0) 
		NumMyNonzeros--; // If I own last global row, then there is one less nonzero
  EPETRA_TEST_ERR(check(A, NumMyEquations, NumGlobalEquations, NumMyNonzeros, 3*NumGlobalEquations-2, 
	       MyGlobalElements, verbose),ierr);
  forierr = 0;
  for (int i = 0; i < NumMyEquations; i++) 
		forierr += !(A.NumGlobalEntries(MyGlobalElements[i])==NumNz[i]+1);
  EPETRA_TEST_ERR(forierr,ierr);
  forierr = 0;
  for (int i = 0; i < NumMyEquations; i++) 
		forierr += !(A.NumMyEntries(i)==NumNz[i]+1);
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "\n\nNumEntries function check OK" << std::endl<< std::endl;

  EPETRA_TEST_ERR(check_graph_sharing(Comm),ierr);

  // Create vectors for Power method

  Epetra_Vector q(Map);
  Epetra_Vector z(Map);
  Epetra_Vector resid(Map);

  // variable needed for iteration
  double lambda = 0.0;
  // int niters = 10000;
  int niters = 200;
  double tolerance = 1.0e-1;

  /////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Iterate

  Epetra_Flops flopcounter;
  A.SetFlopCounter(flopcounter);
  q.SetFlopCounter(A);
  z.SetFlopCounter(A);
  resid.SetFlopCounter(A);
	

  Epetra_Time timer(Comm);
  EPETRA_TEST_ERR(power_method(false, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  double elapsed_time = timer.ElapsedTime();
  double total_flops = A.Flops() + q.Flops() + z.Flops() + resid.Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for first solve = " << MFLOPs << std::endl<< std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////
	
  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << std::endl;
  // Iterate
  lambda = 0.0;
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  EPETRA_TEST_ERR(power_method(true, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops() + q.Flops() + z.Flops() + resid.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for transpose solve = " << MFLOPs << std::endl<< endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Increase diagonal dominance

  if (verbose) cout << "\n\nIncreasing the magnitude of first diagonal term and solving again\n\n"
		    << endl;

  
  if (A.MyGlobalRow(0)) {
    int numvals = A.NumGlobalEntries(0);
    double * Rowvals = new double [numvals];
    long long * Rowinds = new long long[numvals];
    A.ExtractGlobalRowCopy(0, numvals, numvals, Rowvals, Rowinds); // Get A[0,0]

    for (int i=0; i<numvals; i++) if (Rowinds[i] == 0) Rowvals[i] *= 10.0;
    
    A.ReplaceGlobalValues(0, numvals, Rowvals, Rowinds);
    delete [] Rowvals;
    delete [] Rowinds;
  }
  // Iterate (again)
  lambda = 0.0;
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  EPETRA_TEST_ERR(power_method(false, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops() + q.Flops() + z.Flops() + resid.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for second solve = " << MFLOPs << endl<< endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << endl;

  // Iterate (again)
  lambda = 0.0;
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  EPETRA_TEST_ERR(power_method(true, A, q, z, resid, &lambda, niters, tolerance, verbose),ierr);
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops() + q.Flops() + z.Flops() + resid.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;


  if (verbose) cout << "\n\nTotal MFLOPs for tranpose of second solve = " << MFLOPs << endl<< endl;

  if (verbose) cout << "\n\n*****Testing constant entry constructor" << endl<< endl;

  Epetra_CrsMatrix AA(Copy, Map, 5);
  
  if (debug) Comm.Barrier();

  double dble_one = 1.0;
  for (int i=0; i< NumMyEquations; i++) AA.InsertGlobalValues(MyGlobalElements[i], 1, &dble_one, MyGlobalElements+i);

  // Note:  All processors will call the following Insert routines, but only the processor
  //        that owns it will actually do anything

  long long One = 1;
  if (AA.MyGlobalRow(0)) {
    EPETRA_TEST_ERR(!(AA.InsertGlobalValues(0, 0, &dble_one, &One)==0),ierr);
  }
  else EPETRA_TEST_ERR(!(AA.InsertGlobalValues(0, 1, &dble_one, &One)==-1),ierr);
  EPETRA_TEST_ERR(!(AA.FillComplete(false)==0),ierr);
  EPETRA_TEST_ERR(AA.StorageOptimized(),ierr);
  EPETRA_TEST_ERR(!(AA.UpperTriangular()),ierr);
  EPETRA_TEST_ERR(!(AA.LowerTriangular()),ierr);
  
  if (debug) Comm.Barrier();
  EPETRA_TEST_ERR(check(AA, NumMyEquations, NumGlobalEquations, NumMyEquations, NumGlobalEquations, 
	       MyGlobalElements, verbose),ierr);

  if (debug) Comm.Barrier();

  forierr = 0;
  for (int i=0; i<NumMyEquations; i++) forierr += !(AA.NumGlobalEntries(MyGlobalElements[i])==1);
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "\n\nNumEntries function check OK" << endl<< endl;

  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing copy constructor" << endl<< endl;

  Epetra_CrsMatrix B(AA);
  EPETRA_TEST_ERR(check(B, NumMyEquations, NumGlobalEquations, NumMyEquations, NumGlobalEquations, 
	       MyGlobalElements, verbose),ierr);

  forierr = 0;
  for (int i=0; i<NumMyEquations; i++) forierr += !(B.NumGlobalEntries(MyGlobalElements[i])==1);
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "\n\nNumEntries function check OK" << endl<< endl;

  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing local view constructor" << endl<< endl;

  Epetra_CrsMatrix BV(View, AA.RowMap(), AA.ColMap(), 0);

  forierr = 0;
  int* Inds;
  double* Vals;
  for (int i = 0; i < NumMyEquations; i++) {
    forierr += !(AA.ExtractMyRowView(i, NumEntries, Vals, Inds)==0);
    forierr += !(BV.InsertMyValues(i, NumEntries, Vals, Inds)==0);
  }
  BV.FillComplete(false);
  EPETRA_TEST_ERR(check(BV, NumMyEquations, NumGlobalEquations, NumMyEquations, NumGlobalEquations, 
												MyGlobalElements, verbose),ierr);

  forierr = 0;
  for (int i=0; i<NumMyEquations; i++) forierr += !(BV.NumGlobalEntries(MyGlobalElements[i])==1);
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "\n\nNumEntries function check OK" << endl<< endl;

  if (debug) Comm.Barrier();
  if (verbose) cout << "\n\n*****Testing post construction modifications" << endl<< endl;

  EPETRA_TEST_ERR(!(B.InsertGlobalValues(0, 1, &dble_one, &One)==-2),ierr);


  // Release all objects
  delete [] NumNz;
  delete [] Values;
  delete [] Indices;
  delete [] MyGlobalElements;
			

  if (verbose1) {
    // Test ostream << operator (if verbose1)
    // Construct a Map that puts 2 equations on each PE
    
    int NumMyElements1 = 2;
    int NumMyEquations1 = NumMyElements1;
    int NumGlobalEquations1 = NumMyEquations1*NumProc;

    Epetra_Map Map1((long long)-1, NumMyElements1, 0LL, Comm);
    
    // Get update list and number of local equations from newly created Map
    long long* MyGlobalElements1 = new long long[Map1.NumMyElements()];
    Map1.MyGlobalElements(MyGlobalElements1);
    
    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor
    
    int * NumNz1 = new int[NumMyEquations1];
    
    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    
    for (int i=0; i<NumMyEquations1; i++)
      if (MyGlobalElements1[i]==0 || MyGlobalElements1[i] == NumGlobalEquations1-1)
	NumNz1[i] = 1;
      else
	NumNz1[i] = 2;
    
    // Create a Epetra_Matrix
    
    Epetra_CrsMatrix A1(Copy, Map1, NumNz1);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    
    double *Values1 = new double[2];
    Values1[0] = -1.0; Values1[1] = -1.0;
    long long *Indices1 = new long long[2];
    double two1 = 2.0;
    int NumEntries1;

    forierr = 0;
    for (int i=0; i<NumMyEquations1; i++)
      {
	if (MyGlobalElements1[i]==0)
	  {
	    Indices1[0] = 1;
	    NumEntries1 = 1;
	  }
	else if (MyGlobalElements1[i] == NumGlobalEquations1-1)
	  {
	    Indices1[0] = NumGlobalEquations1-2;
	    NumEntries1 = 1;
	  }
	else
	  {
	    Indices1[0] = MyGlobalElements1[i]-1;
	    Indices1[1] = MyGlobalElements1[i]+1;
	    NumEntries1 = 2;
	  }
	forierr += !(A1.InsertGlobalValues(MyGlobalElements1[i], NumEntries1, Values1, Indices1)==0);
	forierr += !(A1.InsertGlobalValues(MyGlobalElements1[i], 1, &two1, MyGlobalElements1+i)>0); // Put in the diagonal entry
      }
    EPETRA_TEST_ERR(forierr,ierr);
    delete [] Indices1;
    delete [] Values1;
    
    // Finish up
    EPETRA_TEST_ERR(!(A1.FillComplete(false)==0),ierr);
    
    // Test diagonal extraction function

    Epetra_Vector checkDiag(Map1);
    EPETRA_TEST_ERR(!(A1.ExtractDiagonalCopy(checkDiag)==0),ierr);

    forierr = 0;
    for (int i=0; i<NumMyEquations1; i++) forierr += !(checkDiag[i]==two1);
    EPETRA_TEST_ERR(forierr,ierr);

    // Test diagonal replacement method

    forierr = 0;
    for (int i=0; i<NumMyEquations1; i++) checkDiag[i]=two1*two1;
    EPETRA_TEST_ERR(forierr,ierr);

    EPETRA_TEST_ERR(!(A1.ReplaceDiagonalValues(checkDiag)==0),ierr);

    Epetra_Vector checkDiag1(Map1);
    EPETRA_TEST_ERR(!(A1.ExtractDiagonalCopy(checkDiag1)==0),ierr);

    forierr = 0;
    for (int i=0; i<NumMyEquations1; i++) forierr += !(checkDiag[i]==checkDiag1[i]);
    EPETRA_TEST_ERR(forierr,ierr);

    if (verbose) cout << "\n\nDiagonal extraction and replacement OK.\n\n" << endl;

    double orignorm = A1.NormOne();
    EPETRA_TEST_ERR(!(A1.Scale(4.0)==0),ierr);
    EPETRA_TEST_ERR(!(A1.NormOne()!=orignorm),ierr);
    
    if (verbose) cout << "\n\nMatrix scale OK.\n\n" << endl;

    if (verbose) cout << "\n\nPrint out tridiagonal matrix, each part on each processor.\n\n" << endl;
    cout << A1 << endl;
   

  // Release all objects
  delete [] NumNz1;
  delete [] MyGlobalElements1;

  }

  if (verbose) cout << "\n\n*****Testing LeftScale and RightScale" << endl << endl;

  int NumMyElements2 = 7;
  int NumMyRows2 = 1;//This value should not be changed without editing the
		// code below.
  Epetra_Map RowMap((long long)-1,NumMyRows2,0LL,Comm);
  Epetra_Map ColMap((long long)NumMyElements2,NumMyElements2,0LL,Comm);
  // The DomainMap needs to be different from the ColMap for the test to 
  // be meaningful.
  Epetra_Map DomainMap((long long)NumMyElements2,0LL,Comm);
  int NumMyRangeElements2 = 0;
  // We need to distribute the elements differently for the range map also.
  if (MyPID % 2 == 0)
    NumMyRangeElements2 = NumMyRows2*2; //put elements on even number procs 
  if (NumProc % 2 == 1 && MyPID == NumProc-1)
    NumMyRangeElements2 = NumMyRows2; //If number of procs is odd, put
			// the last NumMyElements2 elements on the last proc
  Epetra_Map RangeMap((long long)-1,NumMyRangeElements2,0LL,Comm);
  Epetra_CrsMatrix A2(Copy,RowMap,ColMap,NumMyElements2);
  double * Values2 = new double[NumMyElements2];
  int * Indices2 = new int[NumMyElements2]; 

  for (int i=0; i<NumMyElements2; i++) {
    Values2[i] = i+MyPID;
    Indices2[i]=i;
  }

  A2.InsertMyValues(0,NumMyElements2,Values2,Indices2);
  A2.FillComplete(DomainMap,RangeMap,false);
  Epetra_CrsMatrix A2copy(A2);

  double * RowLeftScaleValues = new double[NumMyRows2];
  double * ColRightScaleValues = new double[NumMyElements2];
  long long RowLoopLength = RowMap.MaxMyGID64()-RowMap.MinMyGID64()+1;
  for (long long i=0; i<RowLoopLength; i++)
    RowLeftScaleValues[i] = (double) ((i + RowMap.MinMyGID64() ) % 2 + 1);
  // For the column map, all procs own all elements
  for (int  i=0; i<NumMyElements2;i++)
    ColRightScaleValues[i] = i % 2 + 1;

  long long RangeLoopLength = RangeMap.MaxMyGID64()-RangeMap.MinMyGID64()+1;
  double * RangeLeftScaleValues = new double[(std::size_t) RangeLoopLength];
  long long DomainLoopLength = DomainMap.MaxMyGID64()-DomainMap.MinMyGID64()+1;
   double * DomainRightScaleValues = new double[(std::size_t) DomainLoopLength];
  for (long long i=0; i<RangeLoopLength; i++)
    RangeLeftScaleValues[i] = 1.0/((i + RangeMap.MinMyGID64() ) % 2 + 1);
  for (long long  i=0; i<DomainLoopLength;i++)
    DomainRightScaleValues[i] = 1.0/((i + DomainMap.MinMyGID64() ) % 2 + 1);
                                                                                
  Epetra_Vector xRow(View,RowMap,RowLeftScaleValues);
  Epetra_Vector xCol(View,ColMap,ColRightScaleValues);
  Epetra_Vector xRange(View,RangeMap,RangeLeftScaleValues);
  Epetra_Vector xDomain(View,DomainMap,DomainRightScaleValues);

  double A2infNorm = A2.NormInf();
  double A2oneNorm = A2.NormOne();

  if (verbose1) cout << A2;
  EPETRA_TEST_ERR(A2.LeftScale(xRow),ierr);
  double A2infNorm1 = A2.NormInf();
  double A2oneNorm1 = A2.NormOne();
  bool ScalingBroke = false;
  if (A2infNorm1>2*A2infNorm||A2infNorm1<A2infNorm) {
    EPETRA_TEST_ERR(-31,ierr);
    ScalingBroke = true;
  }
  if (A2oneNorm1>2*A2oneNorm||A2oneNorm1<A2oneNorm) {

    EPETRA_TEST_ERR(-32,ierr);
    ScalingBroke = true;
  }
  if (verbose1) cout << A2;
  EPETRA_TEST_ERR(A2.RightScale(xCol),ierr);
  double A2infNorm2 = A2.NormInf();
  double A2oneNorm2 = A2.NormOne();
  if (A2infNorm2>=2*A2infNorm1||A2infNorm2<=A2infNorm1) {
    EPETRA_TEST_ERR(-33,ierr);
    ScalingBroke = true;
  }
  if (A2oneNorm2>2*A2oneNorm1||A2oneNorm2<=A2oneNorm1) {
    EPETRA_TEST_ERR(-34,ierr);
    ScalingBroke = true;
  }
  if (verbose1) cout << A2;
  EPETRA_TEST_ERR(A2.RightScale(xDomain),ierr);
  double A2infNorm3 = A2.NormInf();
  double A2oneNorm3 = A2.NormOne();
  // The last two scaling ops cancel each other out
  if (A2infNorm3!=A2infNorm1) {
    EPETRA_TEST_ERR(-35,ierr)
    ScalingBroke = true;
  }
  if (A2oneNorm3!=A2oneNorm1) {
    EPETRA_TEST_ERR(-36,ierr)
    ScalingBroke = true;
  }
  if (verbose1) cout << A2;
  EPETRA_TEST_ERR(A2.LeftScale(xRange),ierr);
  double A2infNorm4 = A2.NormInf();
  double A2oneNorm4 = A2.NormOne();
  // The 4 scaling ops all cancel out
  if (A2infNorm4!=A2infNorm) {
    EPETRA_TEST_ERR(-37,ierr)
    ScalingBroke = true;
  }
  if (A2oneNorm4!=A2oneNorm) {
    EPETRA_TEST_ERR(-38,ierr)
    ScalingBroke = true;
  }

  //
  //  Now try changing the values underneath and make sure that 
  //  telling one process about the change causes NormInf() and 
  //  NormOne() to recompute the norm on all processes.
  //
  
  double *values; 
  int num_my_rows = A2.NumMyRows() ; 
  int num_entries;

  for ( int  i=0 ; i< num_my_rows; i++ ) {
    EPETRA_TEST_ERR( A2.ExtractMyRowView( i, num_entries, values ), ierr );
    for ( int j = 0 ; j <num_entries; j++ ) {
      values[j] *= 2.0; 
    }
  }


  if ( MyPID == 0 )
    A2.SumIntoGlobalValues((long long) 0, 0, 0, 0 ) ; 

  double A2infNorm5 = A2.NormInf();
  double A2oneNorm5 = A2.NormOne();

  if (A2infNorm5!=2.0 * A2infNorm4) {
    EPETRA_TEST_ERR(-39,ierr)
    ScalingBroke = true;
  }
  if (A2oneNorm5!= 2.0 * A2oneNorm4) {
    EPETRA_TEST_ERR(-40,ierr)
    ScalingBroke = true;
  }

  //
  //  Restore the values underneath
  //
  for ( int  i=0 ; i< num_my_rows; i++ ) {
    EPETRA_TEST_ERR( A2.ExtractMyRowView( i, num_entries, values ), ierr );
    for ( int j = 0 ; j <num_entries; j++ ) {
      values[j] /= 2.0; 
    }
  }

  if (verbose1) cout << A2;

  if (ScalingBroke) {
    if (verbose) cout << endl << "LeftScale and RightScale tests FAILED" << endl << endl;
  }
  else {
    if (verbose) cout << endl << "LeftScale and RightScale tests PASSED" << endl << endl;
  }

  Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing InvRowMaxs and InvColMaxs" << endl << endl;

  if (verbose1) cout << A2 << endl;
  EPETRA_TEST_ERR(A2.InvRowMaxs(xRow),ierr);
  EPETRA_TEST_ERR(A2.InvRowMaxs(xRange),ierr);
  if (verbose1) cout << xRow << endl << xRange << endl;

  if (verbose) cout << "\n\n*****Testing InvRowSums and InvColSums" << endl << endl;
  bool InvSumsBroke = false;
// Works!
  EPETRA_TEST_ERR(A2.InvRowSums(xRow),ierr);
  if (verbose1) cout << xRow;
  EPETRA_TEST_ERR(A2.LeftScale(xRow),ierr);
  float A2infNormFloat = A2.NormInf();
  if (verbose1) cout << A2 << endl;
  if (fabs(1.0-A2infNormFloat) > 1.e-5) {
    EPETRA_TEST_ERR(-41,ierr);
    InvSumsBroke = true;
  }

  // Works
  int expectedcode = 1;
  if (Comm.NumProc()>1) expectedcode = 0;
  EPETRA_TEST_ERR(!(A2.InvColSums(xDomain)==expectedcode),ierr); // This matrix has a single row, the first column has a zero, so a warning is issued.
  if (verbose1) cout << xDomain << endl;
  EPETRA_TEST_ERR(A2.RightScale(xDomain),ierr);
  float A2oneNormFloat2 = A2.NormOne();
  if (verbose1) cout << A2;
  if (fabs(1.0-A2oneNormFloat2)>1.e-5) {
    EPETRA_TEST_ERR(-42,ierr)
    InvSumsBroke = true;
  }

// Works!
  EPETRA_TEST_ERR(A2.InvRowSums(xRange),ierr);

  if (verbose1) cout << xRange;
  EPETRA_TEST_ERR(A2.LeftScale(xRange),ierr);
  float A2infNormFloat2 = A2.NormInf(); // We use a float so that rounding error
	// will not prevent the sum from being 1.0.
  if (verbose1) cout << A2;
  if (fabs(1.0-A2infNormFloat2)>1.e-5) {
    cout << "InfNorm should be = 1, but InfNorm = " << A2infNormFloat2 << endl;
    EPETRA_TEST_ERR(-43,ierr);
    InvSumsBroke = true;
  }

  // Doesn't work - may not need this test because column ownership is not unique
  /*  EPETRA_TEST_ERR(A2.InvColSums(xCol),ierr);
cout << xCol;
  EPETRA_TEST_ERR(A2.RightScale(xCol),ierr);
  float A2oneNormFloat = A2.NormOne();
cout << A2;
  if (fabs(1.0-A2oneNormFloat)>1.e-5) {
    EPETRA_TEST_ERR(-44,ierr);
    InvSumsBroke = true;
  }
  */
  delete [] ColRightScaleValues;
  delete [] DomainRightScaleValues;
  if (verbose) cout << "Begin partial sum testing." << endl;
  // Test with a matrix that has partial sums for a subset of the rows 
  // on multiple processors. (Except for the serial case, of course.)
  int NumMyRows3 = 2; // Changing this requires further changes below
  long long * myGlobalElements = new long long[NumMyRows3];
  for (int i=0; i<NumMyRows3; i++) myGlobalElements[i] = MyPID+i;
  Epetra_Map RowMap3((long long)NumProc*2, NumMyRows3, myGlobalElements, 0LL, Comm);
  int NumMyElements3 = 5;
  Epetra_CrsMatrix A3(Copy, RowMap3, NumMyElements3);
  double * Values3 = new double[NumMyElements3];
  long long * Indices3 = new long long[NumMyElements3];
  for (int i=0; i < NumMyElements3; i++) {
    Values3[i] = (int) (MyPID + (i+1));
    Indices3[i]=i;
  }
  for (int i=0; i<NumMyRows3; i++) {
    A3.InsertGlobalValues(myGlobalElements[i],NumMyElements3,Values3,Indices3);
  }
  Epetra_Map RangeMap3((long long)NumProc+1, 0LL, Comm);
  Epetra_Map DomainMap3((long long)NumMyElements3, 0LL, Comm);
  EPETRA_TEST_ERR(A3.FillComplete(DomainMap3, RangeMap3,false),ierr);
  if (verbose1) cout << A3;
  Epetra_Vector xRange3(RangeMap3,false);
  Epetra_Vector xDomain3(DomainMap3,false);

  EPETRA_TEST_ERR(A3.InvRowSums(xRange3),ierr);

  if (verbose1) cout << xRange3;
  EPETRA_TEST_ERR(A3.LeftScale(xRange3),ierr);
  float A3infNormFloat = A3.NormInf();
  if (verbose1) cout << A3;
  if (1.0!=A3infNormFloat) {
    cout << "InfNorm should be = 1, but InfNorm = " << A3infNormFloat <<endl;
    EPETRA_TEST_ERR(-61,ierr);
    InvSumsBroke = true;
  }
  // we want to take the transpose of our matrix and fill in different values.
  int NumMyColumns3 = NumMyRows3;
  Epetra_Map ColMap3cm(RowMap3); 
  Epetra_Map RowMap3cm(A3.ColMap());

  Epetra_CrsMatrix A3cm(Copy,RowMap3cm,ColMap3cm,NumProc+1);
  double *Values3cm = new double[NumMyColumns3];
  long long * Indices3cm = new long long[NumMyColumns3];
  for (int i=0; i<NumMyColumns3; i++) {
    Values3cm[i] = MyPID + i + 1;
    Indices3cm[i]= i + MyPID;
  }
  for (int ii=0; ii<NumMyElements3; ii++) {
    A3cm.InsertGlobalValues(ii, NumMyColumns3, Values3cm, Indices3cm);
  }

  // The DomainMap and the RangeMap from the last test will work fine for 
  // the RangeMap and DomainMap, respectively, but I will make copies to
  // avaoid confusion when passing what looks like a DomainMap where we
  // need a RangeMap and vice vera.
  Epetra_Map RangeMap3cm(DomainMap3);
  Epetra_Map DomainMap3cm(RangeMap3);
  EPETRA_TEST_ERR(A3cm.FillComplete(DomainMap3cm,RangeMap3cm),ierr);
  if (verbose1) cout << A3cm << endl;

  // Again, we can copy objects from the last example.
  //Epetra_Vector xRange3cm(xDomain3); //Don't use at this time
  Epetra_Vector xDomain3cm(DomainMap3cm,false);

  EPETRA_TEST_ERR(A3cm.InvColSums(xDomain3cm),ierr);

  if (verbose1) cout << xDomain3cm << endl;

  EPETRA_TEST_ERR(A3cm.RightScale(xDomain3cm),ierr);
  float A3cmOneNormFloat  = A3cm.NormOne();
  if (verbose1) cout << A3cm << endl;
  if (1.0!=A3cmOneNormFloat) {
    cout << "OneNorm should be = 1, but OneNorm = " << A3cmOneNormFloat << endl;
    EPETRA_TEST_ERR(-62,ierr);
    InvSumsBroke = true;
  }
  
  if (verbose) cout << "End partial sum testing" << endl;
  if (verbose) cout << "Begin replicated testing" << endl;

  // We will now view the shared row as a repliated row, rather than one 
  // that has partial sums of its entries on mulitple processors.
  // We will reuse much of the data used for the partial sum tesitng.
  Epetra_Vector xRow3(RowMap3,false); 
  Epetra_CrsMatrix A4(Copy, RowMap3, NumMyElements3);
  for (int ii=0; ii < NumMyElements3; ii++) {
    Values3[ii] = (int)((ii*.6)+1.0);
  }
  for (int ii=0; ii<NumMyRows3; ii++) {
    A4.InsertGlobalValues(myGlobalElements[ii],NumMyElements3,Values3,Indices3);
  }
  EPETRA_TEST_ERR(A4.FillComplete(DomainMap3, RangeMap3,false),ierr);
  if (verbose1) cout << A4 << endl;
  // The next two lines should be expanded into a verifiable test.
  EPETRA_TEST_ERR(A4.InvRowMaxs(xRow3),ierr);
  EPETRA_TEST_ERR(A4.InvRowMaxs(xRange3),ierr);
  if (verbose1) cout << xRow3 << xRange3;

  EPETRA_TEST_ERR(A4.InvRowSums(xRow3),ierr);                      
  if (verbose1) cout << xRow3;
  EPETRA_TEST_ERR(A4.LeftScale(xRow3),ierr);
  float A4infNormFloat = A4.NormInf();
  if (verbose1) cout << A4;
  if (2.0!=A4infNormFloat && NumProc != 1) {
    if (verbose1) cout << "InfNorm should be = 2 (because one column is replicated on two processors and NormOne() does not handle replication), but InfNorm = " << A4infNormFloat <<endl;
    EPETRA_TEST_ERR(-63,ierr);
    InvSumsBroke = true;
  }
  else if (1.0!=A4infNormFloat && NumProc == 1) {
    if (verbose1) cout << "InfNorm should be = 1, but InfNorm = " << A4infNormFloat <<endl;
    EPETRA_TEST_ERR(-63,ierr);
    InvSumsBroke = true;
  }
  
  Epetra_Vector xCol3cm(ColMap3cm,false);
  Epetra_CrsMatrix A4cm(Copy, RowMap3cm, ColMap3cm, NumProc+1);
  //Use values from A3cm
  for (int ii=0; ii<NumMyElements3; ii++) {
    A4cm.InsertGlobalValues(ii,NumMyColumns3,Values3cm,Indices3cm);
  }
  EPETRA_TEST_ERR(A4cm.FillComplete(DomainMap3cm, RangeMap3cm,false),ierr);
  if (verbose1) cout << A4cm << endl;
  // The next two lines should be expanded into a verifiable test.
  EPETRA_TEST_ERR(A4cm.InvColMaxs(xCol3cm),ierr);
  EPETRA_TEST_ERR(A4cm.InvColMaxs(xDomain3cm),ierr);
  if (verbose1) cout << xCol3cm << xDomain3cm;

  EPETRA_TEST_ERR(A4cm.InvColSums(xCol3cm),ierr);
 
  if (verbose1) cout << xCol3cm << endl;
  EPETRA_TEST_ERR(A4cm.RightScale(xCol3cm),ierr);
  float A4cmOneNormFloat = A4cm.NormOne();
  if (verbose1) cout << A4cm << endl;
  if (2.0!=A4cmOneNormFloat && NumProc != 1) {
    if (verbose1) cout << "OneNorm should be = 2 (because one column is replicated on two processors and NormOne() does not handle replication), but OneNorm = " << A4cmOneNormFloat << endl;
    EPETRA_TEST_ERR(-64,ierr);
    InvSumsBroke = true;
  }
  else if (1.0!=A4cmOneNormFloat && NumProc == 1) {
    if (verbose1) cout << "OneNorm should be = 1, but OneNorm = " << A4infNormFloat <<endl;
    EPETRA_TEST_ERR(-64,ierr);
    InvSumsBroke = true;
  }

  if (verbose) cout << "End replicated testing" << endl;

  if (InvSumsBroke) {
    if (verbose) cout << endl << "InvRowSums tests FAILED" << endl << endl;
  }
  else
    if (verbose) cout << endl << "InvRowSums tests PASSED" << endl << endl;

  A3cm.PutScalar(2.0);
  long long nnz_A3cm = A3cm.Graph().NumGlobalNonzeros64();
  double check_frobnorm = sqrt(nnz_A3cm*4.0);
  double frobnorm = A3cm.NormFrobenius();

  bool frobnorm_test_failed = false;
  if (fabs(check_frobnorm-frobnorm) > 5.e-5) {
    frobnorm_test_failed = true;
  }

  if (frobnorm_test_failed) {
    if (verbose) std::cout << "Frobenius-norm test FAILED."<<std::endl;
    EPETRA_TEST_ERR(-65, ierr);
  }

  delete [] Values2;
  delete [] Indices2;
  delete [] myGlobalElements;
  delete [] Values3;
  delete [] Indices3;
  delete [] Values3cm;
  delete [] Indices3cm;
  delete [] RangeLeftScaleValues;
  delete [] RowLeftScaleValues;
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

int power_method(bool TransA, Epetra_CrsMatrix& A, Epetra_Vector& q, Epetra_Vector& z, 
								 Epetra_Vector& resid, double* lambda, int niters, double tolerance, bool verbose) 
{  
	
  // Fill z with random Numbers
  z.Random();
	
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

int check(Epetra_CrsMatrix& A, int NumMyRows1, long long NumGlobalRows1, int NumMyNonzeros1,
					long long NumGlobalNonzeros1, long long* MyGlobalElements, bool verbose) 
{  
  (void)MyGlobalElements;
  int ierr = 0, forierr = 0;
  int NumGlobalIndices;
  int NumMyIndices;
	int* MyViewIndices = 0;
	long long* GlobalViewIndices = 0;
  double* MyViewValues = 0;
	double* GlobalViewValues = 0;
  int MaxNumIndices = A.Graph().MaxNumIndices();
  int* MyCopyIndices = new int[MaxNumIndices];
  long long* GlobalCopyIndices = new long long[MaxNumIndices];
  double* MyCopyValues = new double[MaxNumIndices];
  double* GlobalCopyValues = new double[MaxNumIndices];

  // Test query functions

  int NumMyRows = A.NumMyRows();
  if (verbose) cout << "\n\nNumber of local Rows = " << NumMyRows << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyRows==NumMyRows1),ierr);

  int NumMyNonzeros = A.NumMyNonzeros();
  if (verbose) cout << "\n\nNumber of local Nonzero entries = " << NumMyNonzeros << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyNonzeros==NumMyNonzeros1),ierr);

  long long NumGlobalRows = A.NumGlobalRows64();
  if (verbose) cout << "\n\nNumber of global Rows = " << NumGlobalRows << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalRows==NumGlobalRows1),ierr);

  long long NumGlobalNonzeros = A.NumGlobalNonzeros64();
  if (verbose) cout << "\n\nNumber of global Nonzero entries = " << NumGlobalNonzeros << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalNonzeros==NumGlobalNonzeros1),ierr);

  // GlobalRowView should be illegal (since we have local indices)

  EPETRA_TEST_ERR(!(A.ExtractGlobalRowView(A.RowMap().MaxMyGID64(), NumGlobalIndices, GlobalViewValues, GlobalViewIndices)==-2),ierr);

  // Other binary tests

  EPETRA_TEST_ERR(A.NoDiagonal(),ierr);
  EPETRA_TEST_ERR(!(A.Filled()),ierr);
  EPETRA_TEST_ERR(!(A.MyGRID(A.RowMap().MaxMyGID64())),ierr);
  EPETRA_TEST_ERR(!(A.MyGRID(A.RowMap().MinMyGID64())),ierr);
  EPETRA_TEST_ERR(A.MyGRID(1+A.RowMap().MaxMyGID64()),ierr);
  EPETRA_TEST_ERR(A.MyGRID(-1+A.RowMap().MinMyGID64()),ierr);
  EPETRA_TEST_ERR(!(A.MyLRID(0)),ierr);
  EPETRA_TEST_ERR(!(A.MyLRID(NumMyRows-1)),ierr);
  EPETRA_TEST_ERR(A.MyLRID(-1),ierr);
  EPETRA_TEST_ERR(A.MyLRID(NumMyRows),ierr);

  forierr = 0;
  for (int i = 0; i < NumMyRows; i++) {
    long long Row = A.GRID64(i);
    A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumGlobalIndices, GlobalCopyValues, GlobalCopyIndices);
    A.ExtractMyRowView(i, NumMyIndices, MyViewValues, MyViewIndices); // this is where the problem comes from
    forierr += !(NumGlobalIndices == NumMyIndices);
    for(int j = 1; j < NumMyIndices; j++) {
			forierr += !(MyViewIndices[j-1] < MyViewIndices[j]); // this is where the test fails
		}
    for(int j = 0; j < NumGlobalIndices; j++) {
			forierr += !(GlobalCopyIndices[j] == A.GCID64(MyViewIndices[j]));
			forierr += !(A.LCID(GlobalCopyIndices[j]) == MyViewIndices[j]);
			forierr += !(GlobalCopyValues[j] == MyViewValues[j]);
    }
  }
  EPETRA_TEST_ERR(forierr,ierr);

  forierr = 0;
  for (int i = 0; i < NumMyRows; i++) {
    long long Row = A.GRID64(i);
    A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumGlobalIndices, GlobalCopyValues, GlobalCopyIndices);
    A.ExtractMyRowCopy(i, MaxNumIndices, NumMyIndices, MyCopyValues, MyCopyIndices);
    forierr += !(NumGlobalIndices == NumMyIndices);
    for (int j = 1; j < NumMyIndices; j++) 
			forierr += !(MyCopyIndices[j-1] < MyCopyIndices[j]);
    for (int j = 0; j < NumGlobalIndices; j++) {
			forierr += !(GlobalCopyIndices[j] == A.GCID64(MyCopyIndices[j]));
			forierr += !(A.LCID(GlobalCopyIndices[j]) == MyCopyIndices[j]);
			forierr += !(GlobalCopyValues[j] == MyCopyValues[j]);
    }

  }
  EPETRA_TEST_ERR(forierr,ierr);

  delete [] MyCopyIndices;
  delete [] GlobalCopyIndices;
  delete [] MyCopyValues;
  delete [] GlobalCopyValues;

  if (verbose) cout << "\n\nRows sorted check OK" << endl<< endl;

  return (ierr);
}

int check_graph_sharing(Epetra_Comm& Comm)
{
  int numLocalElems = 5;
  int localProc = Comm.MyPID();
  long long firstElem = localProc*numLocalElems;
  int err;
  Epetra_Map map((long long)-1, numLocalElems, 0LL, Comm);

  Epetra_CrsMatrix* A = new Epetra_CrsMatrix(Copy, map, 1);

  for (int i=0; i<numLocalElems; ++i) {
    long long row = firstElem+i;
    long long col = row;
    double val = 1.0;

    err = A->InsertGlobalValues(row, 1, &val, &col);
    if (err != 0) {
      cerr << "A->InsertGlobalValues("<<row<<") returned err="<<err<<endl;
      return(err);
    }
  }

  A->FillComplete(false);

  Epetra_CrsMatrix B(Copy, A->Graph());

  delete A;

  for (int i=0; i<numLocalElems; ++i) {
    long long row = firstElem+i;
    long long col = row;
    double val = 1.0;

    err = B.ReplaceGlobalValues(row, 1, &val, &col);
    if (err != 0) {
      cerr << "B.ReplaceGlobalValues("<<row<<") returned err="<<err<<endl;
      return(err);
    }
  }

  return(0);
}


