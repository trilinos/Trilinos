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

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_Flops.h"
#include "Epetra_VbrMatrix.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"

// prototypes

int CompareValues(double * A, int LDA, int NumRowsA, int NumColsA, 
		  double * B, int LDB, int NumRowsB, int NumColsB);

int check(Epetra_VbrMatrix& A, 
	  int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1, int NumGlobalNonzeros1, 
	  int NumMyBlockRows1, int NumGlobalBlockRows1, int NumMyBlockNonzeros1, int NumGlobalBlockNonzeros1, 
	  int * MyGlobalElements, bool verbose);

int power_method(bool TransA, Epetra_VbrMatrix& A, 
		 Epetra_MultiVector& q,
		 Epetra_MultiVector& z, 
		 Epetra_MultiVector& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose);

int checkMergeRedundantEntries(Epetra_Comm& comm, bool verbose);

int checkExtractMyRowCopy(Epetra_Comm& comm, bool verbose);

int main(int argc, char *argv[])
{
  int ierr = 0, i, j, forierr = 0;
  bool debug = false;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  //char tmp;
  //if (rank==0) cout << "Press any key to continue..."<< endl;
  //if (rank==0) cin >> tmp;
  //Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

//  int NumMyElements = 1000;
  int NumMyElements = 10;
  if (MyPID < 3) NumMyElements++;

  // Define pseudo-random block sizes using a Petra Vector of random numbers
  Epetra_Map randmap(-1, NumMyElements, 0, Comm);
  Epetra_Vector randvec(randmap);
  randvec.Random(); // Fill with random numbers
  int * ElementSizeList = new int[NumMyElements];
  int MinSize = 3;
  int MaxSize = 8;
  int SizeRange = MaxSize - MinSize + 1;
  
  for (i=0; i<NumMyElements; i++) ElementSizeList[i] = 3 + SizeRange * (int) (fabs(randvec[i]) * .99);

  // Construct a Map

  int *randMyGlobalElements = randmap.MyGlobalElements();

  Epetra_BlockMap Map (-1, NumMyElements, randMyGlobalElements, ElementSizeList, 0, Comm);
  
  // Get update list and number of local elements from newly created Map
  int NumGlobalElements = Map.NumGlobalElements();
  int * MyGlobalElements = Map.MyGlobalElements();
  bool DistributedGlobal = Map.DistributedGlobal();

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  int * NumNz = new int[NumMyElements];

  // We are building a block tridiagonal matrix

  for (i=0; i<NumMyElements; i++)
    if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalElements-1)
      NumNz[i] = 2;
    else
      NumNz[i] = 3;
  // Create a Epetra_Matrix

  Epetra_VbrMatrix A(Copy, Map, NumNz);
  EPETRA_TEST_ERR(A.IndicesAreGlobal(),ierr);
  EPETRA_TEST_ERR(A.IndicesAreLocal(),ierr);
  
  // Use an array of Epetra_SerialDenseMatrix objects to build VBR matrix

  Epetra_SerialDenseMatrix ** BlockEntries = new Epetra_SerialDenseMatrix*[SizeRange];

  // The array of dense matrices will increase in size from MinSize to MaxSize (defined above)
  for (int kr=0; kr<SizeRange; kr++) {
    BlockEntries[kr] = new Epetra_SerialDenseMatrix[SizeRange];
    int RowDim = ElementSizeList[kr];
    for (int kc = 0; kc<SizeRange; kc++) {
      int ColDim = ElementSizeList[kc];
      Epetra_SerialDenseMatrix * curmat = &(BlockEntries[kr][kc]);
      curmat->Shape(RowDim,ColDim);
    for (j=0; j < ColDim; j++)
      for (i=0; i < RowDim; i++) {
	BlockEntries[kr][kc][j][i] = -1.0;
	if (i==j && kr==kc) BlockEntries[kr][kc][j][i] = 9.0;
	else BlockEntries[kr][kc][j][i] = -1.0;
      }
    }
  }
  
    
  // Add  rows one-at-a-time

  int *Indices = new int[3];
  int *ColDims = new int[3];
  int NumEntries;
  int NumMyNonzeros = 0, NumMyEquations = 0;
  
  for (i=0; i<NumMyElements; i++) {
    int CurRow = MyGlobalElements[i];
    int RowDim = ElementSizeList[i]-MinSize;
    NumMyEquations += BlockEntries[RowDim][RowDim].M();
    
    if (CurRow==0)
      {
	Indices[0] = CurRow;
	Indices[1] = CurRow+1;
	NumEntries = 2;
	ColDims[0] = ElementSizeList[i] - MinSize;
	ColDims[1] = ElementSizeList[i+1] - MinSize; // Assumes linear global ordering and > 1 row/proc.
      }
    else if (CurRow == NumGlobalElements-1)
      {
	Indices[0] = CurRow-1;
	Indices[1] = CurRow;
	NumEntries = 2;
	ColDims[0] = ElementSizeList[i-1] - MinSize;
	ColDims[1] = ElementSizeList[i] - MinSize; // Assumes linear global ordering and > 1 row/proc.
      }
      else {
	Indices[0] = CurRow-1;
	Indices[1] = CurRow;
	Indices[2] = CurRow+1;
	NumEntries = 3;
	if (i==0) ColDims[0] = EPETRA_MAX(MinSize, EPETRA_MIN(MaxSize, MyPID-1)) - MinSize; // ElementSize on MyPID-1
	else ColDims[0] = ElementSizeList[i-1];
	ColDims[1] = ElementSizeList[i];
	// ElementSize on MyPID+1
	if (i==NumMyElements-1) ColDims[2] = EPETRA_MAX(MinSize, EPETRA_MIN(MaxSize, MyPID)) - MinSize;
	else ColDims[2] = ElementSizeList[i+1] - MinSize;
      }
    EPETRA_TEST_ERR(!(A.BeginInsertGlobalValues(CurRow, NumEntries, Indices)==0),ierr);
    forierr = 0;
    for (j=0; j < NumEntries; j++) {
      Epetra_SerialDenseMatrix * AD = &(BlockEntries[RowDim][ColDims[j]]);
      NumMyNonzeros += AD->M() * AD->N();	  
      forierr += !(A.SubmitBlockEntry(AD->A(), AD->LDA(), AD->M(), AD->N())==0);
    }
    EPETRA_TEST_ERR(forierr,ierr);

      A.EndSubmitEntries();
  }
  
  // Finish up
  EPETRA_TEST_ERR(!(A.IndicesAreGlobal()),ierr);
  EPETRA_TEST_ERR(!(A.FillComplete()==0),ierr);
  EPETRA_TEST_ERR(!(A.IndicesAreLocal()),ierr);
  EPETRA_TEST_ERR(A.StorageOptimized(),ierr);
  // A.OptimizeStorage();
  // EPETRA_TEST_ERR(!(A.StorageOptimized()),ierr);
  EPETRA_TEST_ERR(A.UpperTriangular(),ierr);
  EPETRA_TEST_ERR(A.LowerTriangular(),ierr);


  {for (int kr=0; kr<SizeRange; kr++) delete [] BlockEntries[kr];}
  delete [] BlockEntries;
  delete [] ColDims;
  delete [] Indices;
  delete [] ElementSizeList;
  

  int NumMyBlockEntries = 3*NumMyElements;
  if (A.LRID(0)>=0) NumMyBlockEntries--; // If I own first global row, then there is one less nonzero
  if (A.LRID(NumGlobalElements-1)>=0) NumMyBlockEntries--; // If I own last global row, then there is one less nonzero

  int NumGlobalBlockEntries = 3*NumGlobalElements-2;
  int NumGlobalNonzeros, NumGlobalEquations;
  Comm.SumAll(&NumMyNonzeros, &NumGlobalNonzeros, 1);
  Comm.SumAll(&NumMyEquations, &NumGlobalEquations, 1);
  

  EPETRA_TEST_ERR(!(check(A, NumMyEquations, NumGlobalEquations, NumMyNonzeros, NumGlobalNonzeros, 
	       NumMyElements, NumGlobalElements, NumMyBlockEntries, NumGlobalBlockEntries, 
	       MyGlobalElements, verbose)==0),ierr);
  forierr = 0;
  for (i=0; i<NumMyElements; i++) forierr += !(A.NumGlobalBlockEntries(MyGlobalElements[i])==NumNz[i]);
  EPETRA_TEST_ERR(forierr,ierr);
  forierr = 0;
  for (i=0; i<NumMyElements; i++) forierr += !(A.NumMyBlockEntries(i)==NumNz[i]);
  EPETRA_TEST_ERR(forierr,ierr);

  if (verbose) cout << "\n\nNumEntries function check OK" << endl<< endl;

  delete [] NumNz;


  // Create vectors for Power method

  Epetra_Vector q(Map);
  Epetra_Vector z(Map);
  Epetra_Vector z_initial(Map);
  Epetra_Vector resid(Map);

  
  // Fill z with random Numbers 
  z_initial.Random();

  // variable needed for iteration
  double lambda = 0.0;
  int niters = 100;
  // int niters = 200;
  double tolerance = 1.0e-3;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Iterate
  Epetra_Flops flopcounter;
  A.SetFlopCounter(flopcounter);
  q.SetFlopCounter(A);
  z.SetFlopCounter(A);
  resid.SetFlopCounter(A);
  z = z_initial;  // Start with common initial guess
  Epetra_Time timer(Comm);
  int ierr1 = power_method(false, A, q, z, resid, &lambda, niters, tolerance, verbose);
  double elapsed_time = timer.ElapsedTime();
  double total_flops = flopcounter.Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for first solve = " << MFLOPs << endl<< endl;
  if (verbose && ierr1==1) cout << "***** Power Method did not converge. *****" << endl << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << endl;
  // Iterate
  lambda = 0.0;
  z = z_initial;  // Start with common initial guess
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  ierr1 = power_method(true, A, q, z, resid, &lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();
  total_flops = flopcounter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for transpose solve = " << MFLOPs << endl<< endl;
  if (verbose && ierr1==1) cout << "***** Power Method did not converge. *****" << endl << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Increase diagonal dominance

  if (verbose) cout << "\n\nIncreasing the magnitude of first diagonal term and solving again\n\n"
		    << endl;

  
  if (A.MyGlobalBlockRow(0)) {
    int numvals = A.NumGlobalBlockEntries(0);
    Epetra_SerialDenseMatrix ** Rowvals;
    int* Rowinds = new int[numvals];
    int  RowDim;
    A.ExtractGlobalBlockRowPointers(0, numvals, RowDim, numvals, Rowinds, 
				    Rowvals); // Get A[0,:]

    for (i=0; i<numvals; i++) {
      if (Rowinds[i] == 0) {
	Rowvals[i]->A()[0] *= 10.0; // Multiply first diag value by 10.0
      }
    }
    delete [] Rowinds;
  }
  // Iterate (again)
  lambda = 0.0;
  z = z_initial;  // Start with common initial guess
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  ierr1 = power_method(false, A, q, z, resid, &lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();
  total_flops = flopcounter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for second solve = " << MFLOPs << endl<< endl;
  if (verbose && ierr1==1) cout << "***** Power Method did not converge. *****" << endl << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // Solve transpose problem

  if (verbose) cout << "\n\nUsing transpose of matrix and solving again (should give same result).\n\n"
		    << endl;

  // Iterate (again)
  lambda = 0.0;
  z = z_initial;  // Start with common initial guess
  flopcounter.ResetFlops();
  timer.ResetStartTime();
  ierr1 = power_method(true, A, q, z, resid, &lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();
  total_flops = flopcounter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for tranpose of second solve = " << MFLOPs << endl<< endl;
  if (verbose && ierr1==1) cout << "***** Power Method did not converge. *****" << endl << endl;


  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing copy constructor" << endl<< endl;

  Epetra_VbrMatrix B(A);

  EPETRA_TEST_ERR(!(check(B, NumMyEquations, NumGlobalEquations, NumMyNonzeros, NumGlobalNonzeros, 
	       NumMyElements, NumGlobalElements, NumMyBlockEntries, NumGlobalBlockEntries, 
	       MyGlobalElements, verbose)==0),ierr);


  if (debug) Comm.Barrier();

  if (verbose) cout << "\n\n*****Testing post construction modifications" << endl<< endl;

  int One = 1;
  if (B.MyGRID(0)) EPETRA_TEST_ERR(!(B.BeginInsertGlobalValues(0, 1, &One)==-2),ierr);

  Epetra_Vector checkDiag(B.RowMap());
  forierr = 0;
  int NumMyEquations1 = B.NumMyRows();
  double two1 = 2.0;

    // Test diagonal replacement and extraction methods

    forierr = 0;
    for (i=0; i<NumMyEquations1; i++) checkDiag[i]=two1;
    EPETRA_TEST_ERR(forierr,ierr);

    EPETRA_TEST_ERR(!(B.ReplaceDiagonalValues(checkDiag)==0),ierr);

    Epetra_Vector checkDiag1(B.RowMap());
    EPETRA_TEST_ERR(!(B.ExtractDiagonalCopy(checkDiag1)==0),ierr);

    forierr = 0;
    for (i=0; i<NumMyEquations1; i++) forierr += !(checkDiag[i]==checkDiag1[i]);
    EPETRA_TEST_ERR(forierr,ierr);

    if (verbose) cout << "\n\nDiagonal extraction and replacement OK.\n\n" << endl;

    double orignorm = B.NormOne();
    EPETRA_TEST_ERR(!(B.Scale(4.0)==0),ierr);
    EPETRA_TEST_ERR((B.NormOne()!=orignorm),ierr);
    
    if (verbose) cout << "\n\nMatrix scale OK.\n\n" << endl;


  /*
  if (verbose1) {
    // Test ostream << operator (if verbose1)
    // Construct a Map that puts 2 equations on each PE
    
    int NumMyElements1 = 2;
    int NumMyElements1 = NumMyElements1;
    int NumGlobalElements1 = NumMyElements1*NumProc;

    Epetra_Map Map1(-1, NumMyElements1, 0, Comm);
    
    // Get update list and number of local equations from newly created Map
    int * MyGlobalElements1 = new int[Map1.NumMyElements()];
    Map1.MyGlobalElements(MyGlobalElements1);
    
    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor
    
    int * NumNz1 = new int[NumMyElements1];
    
    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    
    for (i=0; i<NumMyElements1; i++)
      if (MyGlobalElements1[i]==0 || MyGlobalElements1[i] == NumGlobalElements1-1)
	NumNz1[i] = 1;
      else
	NumNz1[i] = 2;
    
    // Create a Epetra_Matrix
    
    Epetra_VbrMatrix A1(Copy, Map1, NumNz1);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    
    int *Indices1 = new int[2];
    double two1 = 2.0;
    int NumEntries1;

    forierr = 0;
    for (i=0; i<NumMyElements1; i++)
      {
	if (MyGlobalElements1[i]==0)
	  {
	    Indices1[0] = 1;
	    NumEntries1 = 1;
	  }
	else if (MyGlobalElements1[i] == NumGlobalElements1-1)
	  {
	    Indices1[0] = NumGlobalElements1-2;
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
    // Finish up
    EPETRA_TEST_ERR(!(A1.FillComplete()==0),ierr);
    
    if (verbose) cout << "\n\nPrint out tridiagonal matrix, each part on each processor.\n\n" << endl;
    cout << A1 << endl;
    
    // Release all objects
    delete [] NumNz1;
    delete [] Values1;
    delete [] Indices1;
    delete [] MyGlobalElements1;

  }
  */

  EPETRA_TEST_ERR( checkMergeRedundantEntries(Comm, verbose1), ierr);

  EPETRA_TEST_ERR( checkExtractMyRowCopy(Comm, verbose1), ierr);

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

int power_method(bool TransA, Epetra_VbrMatrix& A, 
		 Epetra_MultiVector& q,
		 Epetra_MultiVector& z, 
		 Epetra_MultiVector& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose) {  

  // variable needed for iteration
  double normz, residual;

  int ierr = 1;

  for (int iter = 0; iter < niters; iter++)
    {
      z.Norm2(&normz); // Compute 2-norm of z
      q.Scale(1.0/normz, z);
      A.Multiply(TransA, q, z); // Compute z = A*q
      q.Dot(z, lambda); // Approximate maximum eigenvaluE
      if (iter%100==0 || iter+1==niters)
	{
	  resid.Update(1.0, z, -(*lambda), q, 0.0); // Compute A*q - lambda*q
	  resid.Norm2(&residual);
	  if (verbose) cout << "Iter = " << iter << "  Lambda = " << *lambda 
			     << "  Residual of A*q - lambda*q = " << residual << endl;
	} 
      if (residual < tolerance) {
	ierr = 0;
	break;
      }
    }
  return(ierr);
}
int check(Epetra_VbrMatrix& A, 
	  int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1, int NumGlobalNonzeros1, 
	  int NumMyBlockRows1, int NumGlobalBlockRows1, int NumMyBlockNonzeros1, int NumGlobalBlockNonzeros1, 
	  int * MyGlobalElements, bool verbose) {

  int ierr = 0, forierr = 0;
  // Test query functions

  int NumMyRows = A.NumMyRows();
  if (verbose) cout << "\n\nNumber of local Rows = " << NumMyRows << endl<< endl;
  // TEMP
  if (verbose) cout << "\n\nNumber of local Rows should = " << NumMyRows1 << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyRows==NumMyRows1),ierr);

  int NumMyNonzeros = A.NumMyNonzeros();
  if (verbose) cout << "\n\nNumber of local Nonzero entries = " << NumMyNonzeros << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyNonzeros==NumMyNonzeros1),ierr);

  int NumGlobalRows = A.NumGlobalRows();
  if (verbose) cout << "\n\nNumber of global Rows = " << NumGlobalRows << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalRows==NumGlobalRows1),ierr);

  int NumGlobalNonzeros = A.NumGlobalNonzeros();
  if (verbose) cout << "\n\nNumber of global Nonzero entries = " << NumGlobalNonzeros << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalNonzeros==NumGlobalNonzeros1),ierr);

  int NumMyBlockRows = A.NumMyBlockRows();
  if (verbose) cout << "\n\nNumber of local Block Rows = " << NumMyBlockRows << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyBlockRows==NumMyBlockRows1),ierr);

  int NumMyBlockNonzeros = A.NumMyBlockEntries();
  if (verbose) cout << "\n\nNumber of local Nonzero Block entries = " << NumMyBlockNonzeros << endl<< endl;

  EPETRA_TEST_ERR(!(NumMyBlockNonzeros==NumMyBlockNonzeros1),ierr);

  int NumGlobalBlockRows = A.NumGlobalBlockRows();
  if (verbose) cout << "\n\nNumber of global Block Rows = " << NumGlobalBlockRows << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalBlockRows==NumGlobalBlockRows1),ierr);

  int NumGlobalBlockNonzeros = A.NumGlobalBlockEntries();
  if (verbose) cout << "\n\nNumber of global Nonzero Block entries = " << NumGlobalBlockNonzeros << endl<< endl;

  EPETRA_TEST_ERR(!(NumGlobalNonzeros==NumGlobalNonzeros1),ierr);

  
  // Test RowMatrix interface implementations
  int RowDim, NumBlockEntries, * BlockIndices;
  Epetra_SerialDenseMatrix ** Values;
  // Get View of last block row
  A.ExtractMyBlockRowView(NumMyBlockRows-1, RowDim, NumBlockEntries,
			  BlockIndices, Values);
  int NumMyEntries1 = 0;
  {for (int i=0; i < NumBlockEntries; i++) NumMyEntries1 += Values[i]->N();}
  int NumMyEntries;
  A.NumMyRowEntries(NumMyRows-1, NumMyEntries);
  if (verbose) {
    cout << "\n\nNumber of nonzero values in last row = "
	 << NumMyEntries << endl<< endl;
  }

  EPETRA_TEST_ERR(!(NumMyEntries==NumMyEntries1),ierr);
  
  // Other binary tests

  EPETRA_TEST_ERR(A.NoDiagonal(),ierr);
  EPETRA_TEST_ERR(!(A.Filled()),ierr);
  EPETRA_TEST_ERR(!(A.Sorted()),ierr);
  EPETRA_TEST_ERR(!(A.MyGRID(A.RowMap().MaxMyGID())),ierr);
  EPETRA_TEST_ERR(!(A.MyGRID(A.RowMap().MinMyGID())),ierr);
  EPETRA_TEST_ERR(A.MyGRID(1+A.RowMap().MaxMyGID()),ierr);
  EPETRA_TEST_ERR(A.MyGRID(-1+A.RowMap().MinMyGID()),ierr);
  EPETRA_TEST_ERR(!(A.MyLRID(0)),ierr);
  EPETRA_TEST_ERR(!(A.MyLRID(NumMyBlockRows-1)),ierr);
  EPETRA_TEST_ERR(A.MyLRID(-1),ierr);
  EPETRA_TEST_ERR(A.MyLRID(NumMyBlockRows),ierr);

    
  int i, j;
  int MaxNumBlockEntries = A.MaxNumBlockEntries();

  // Pointer Extraction approach

  //   Local index
  int MyPointersRowDim, MyPointersNumBlockEntries;
  int * MyPointersBlockIndices = new int[MaxNumBlockEntries];
  Epetra_SerialDenseMatrix **MyPointersValuesPointers;
  //   Global Index
  int GlobalPointersRowDim, GlobalPointersNumBlockEntries;
  int * GlobalPointersBlockIndices = new int[MaxNumBlockEntries];
  Epetra_SerialDenseMatrix **GlobalPointersValuesPointers;

  // Copy Extraction approach

  //   Local index
  int MyCopyRowDim, MyCopyNumBlockEntries;
  int * MyCopyBlockIndices = new int[MaxNumBlockEntries];
  int * MyCopyColDims = new int[MaxNumBlockEntries];
  int * MyCopyLDAs = new int[MaxNumBlockEntries];
  int MaxRowDim = A.MaxRowDim();
  int MaxColDim = A.MaxColDim();
  int MyCopySizeOfValues = MaxRowDim*MaxColDim;
  double ** MyCopyValuesPointers = new double*[MaxNumBlockEntries];
  for (i=0; i<MaxNumBlockEntries; i++)
    MyCopyValuesPointers[i] = new double[MaxRowDim*MaxColDim];

  //   Global Index
  int GlobalCopyRowDim, GlobalCopyNumBlockEntries;
  int * GlobalCopyBlockIndices = new int[MaxNumBlockEntries];
  int * GlobalCopyColDims = new int[MaxNumBlockEntries];
  int * GlobalCopyLDAs = new int[MaxNumBlockEntries];
  
  int GlobalMaxRowDim = A.GlobalMaxRowDim();
  int GlobalMaxColDim = A.GlobalMaxColDim();
  int GlobalCopySizeOfValues = GlobalMaxRowDim*GlobalMaxColDim;
  double ** GlobalCopyValuesPointers = new double*[MaxNumBlockEntries];
  for (i=0; i<MaxNumBlockEntries; i++)
    GlobalCopyValuesPointers[i] = new double[GlobalMaxRowDim*GlobalMaxColDim];

  // View Extraction approaches

  //   Local index (There is no global view available)
  int MyView1RowDim, MyView1NumBlockEntries;
  int * MyView1BlockIndices;
  Epetra_SerialDenseMatrix **MyView1ValuesPointers = new Epetra_SerialDenseMatrix*[MaxNumBlockEntries];

  //   Local index version 2 (There is no global view available)
  int MyView2RowDim, MyView2NumBlockEntries;
  int * MyView2BlockIndices;
  Epetra_SerialDenseMatrix **MyView2ValuesPointers;


  // For each row, test six approaches to extracting data from a given local index matrix
  forierr = 0;
  for (i=0; i<NumMyBlockRows; i++) {
    int MyRow = i;
    int GlobalRow = A.GRID(i);
    // Get a copy of block indices in local index space, pointers to everything else
    A.ExtractMyBlockRowPointers(MyRow, MaxNumBlockEntries, MyPointersRowDim, 
				MyPointersNumBlockEntries, MyPointersBlockIndices,
				MyPointersValuesPointers);
    // Get a copy of block indices in local index space, pointers to everything else
    A.ExtractGlobalBlockRowPointers(GlobalRow, MaxNumBlockEntries, GlobalPointersRowDim, 
				    GlobalPointersNumBlockEntries, GlobalPointersBlockIndices,
				    GlobalPointersValuesPointers);

    // Initiate a copy of block row in local index space.
    A.BeginExtractMyBlockRowCopy(MyRow, MaxNumBlockEntries, MyCopyRowDim, 
				 MyCopyNumBlockEntries, MyCopyBlockIndices,
				 MyCopyColDims);
    // Copy Values
    for (j=0; j<MyCopyNumBlockEntries; j++) {
      A.ExtractEntryCopy(MyCopySizeOfValues, MyCopyValuesPointers[j], MaxRowDim, false);
      MyCopyLDAs[j] = MaxRowDim;
    }

    // Initiate a copy of block row in global index space.
    A.BeginExtractGlobalBlockRowCopy(GlobalRow, MaxNumBlockEntries, GlobalCopyRowDim, 
				    GlobalCopyNumBlockEntries, GlobalCopyBlockIndices,
				    GlobalCopyColDims);
    // Copy Values
    for (j=0; j<GlobalCopyNumBlockEntries; j++) {
      A.ExtractEntryCopy(GlobalCopySizeOfValues, GlobalCopyValuesPointers[j], GlobalMaxRowDim, false);
      GlobalCopyLDAs[j] = GlobalMaxRowDim;
    }

    // Initiate a view of block row in local index space (Version 1)
    A.BeginExtractMyBlockRowView(MyRow, MyView1RowDim, 
				 MyView1NumBlockEntries, MyView1BlockIndices);
    // Set pointers to values
    for (j=0; j<MyView1NumBlockEntries; j++) 
      A.ExtractEntryView(MyView1ValuesPointers[j]);


    // Extract a view of block row in local index space (version 2)
    A.ExtractMyBlockRowView(MyRow, MyView2RowDim, 
			    MyView2NumBlockEntries, MyView2BlockIndices,
			    MyView2ValuesPointers);

    forierr += !(MyPointersNumBlockEntries==GlobalPointersNumBlockEntries);
    forierr += !(MyPointersNumBlockEntries==MyCopyNumBlockEntries);
    forierr += !(MyPointersNumBlockEntries==GlobalCopyNumBlockEntries);
    forierr += !(MyPointersNumBlockEntries==MyView1NumBlockEntries);
    forierr += !(MyPointersNumBlockEntries==MyView2NumBlockEntries);
    for (j=1; j<MyPointersNumBlockEntries; j++) {
      forierr += !(MyCopyBlockIndices[j-1]<MyCopyBlockIndices[j]);
      forierr += !(MyView1BlockIndices[j-1]<MyView1BlockIndices[j]);
      forierr += !(MyView2BlockIndices[j-1]<MyView2BlockIndices[j]);

      forierr += !(GlobalPointersBlockIndices[j]==A.GCID(MyPointersBlockIndices[j]));
      forierr += !(A.LCID(GlobalPointersBlockIndices[j])==MyPointersBlockIndices[j]);
      forierr += !(GlobalPointersBlockIndices[j]==GlobalCopyBlockIndices[j]);
      
      Epetra_SerialDenseMatrix* my = MyPointersValuesPointers[j];
      Epetra_SerialDenseMatrix* global = GlobalPointersValuesPointers[j];

      Epetra_SerialDenseMatrix* myview1 = MyView1ValuesPointers[j];
      Epetra_SerialDenseMatrix* myview2 = MyView2ValuesPointers[j];

      forierr += !(CompareValues(my->A(), my->LDA(), 
			   MyPointersRowDim, my->N(), 
			   global->A(), global->LDA(), 
			   GlobalPointersRowDim, global->N())==0);
      forierr += !(CompareValues(my->A(), my->LDA(), 
			   MyPointersRowDim, my->N(), 
			   MyCopyValuesPointers[j], MyCopyLDAs[j], 
			   MyCopyRowDim, MyCopyColDims[j])==0);
      forierr += !(CompareValues(my->A(), my->LDA(), 
			   MyPointersRowDim, my->N(), 
			   GlobalCopyValuesPointers[j], GlobalCopyLDAs[j], 
			   GlobalCopyRowDim, GlobalCopyColDims[j])==0);
      forierr += !(CompareValues(my->A(), my->LDA(), 
			   MyPointersRowDim, my->N(), 
			   myview1->A(), myview1->LDA(), 
			   MyView1RowDim, myview1->N())==0);
      forierr += !(CompareValues(my->A(), my->LDA(),
			   MyPointersRowDim, my->N(),
				 myview2->A(), myview2->LDA(),
			   MyView2RowDim, myview2->N())==0);
    }
  }
  EPETRA_TEST_ERR(forierr,ierr);

  // GlobalRowView should be illegal (since we have local indices)
  EPETRA_TEST_ERR(!(A.BeginExtractGlobalBlockRowView(A.GRID(0), MyView1RowDim, 
						     MyView1NumBlockEntries,
						     MyView1BlockIndices)==-2),ierr);
  
  // Extract a view of block row in local index space (version 2)
  EPETRA_TEST_ERR(!(A.ExtractGlobalBlockRowView(A.GRID(0), MyView2RowDim, 
				     MyView2NumBlockEntries, MyView2BlockIndices,
				     MyView2ValuesPointers)==-2),ierr);
  
  delete [] MyPointersBlockIndices;
  delete [] GlobalPointersBlockIndices;
  delete [] MyCopyBlockIndices;
  delete [] MyCopyColDims;
  delete [] MyCopyLDAs;
  for (i=0; i<MaxNumBlockEntries; i++) delete [] MyCopyValuesPointers[i];
  delete [] MyCopyValuesPointers;
  delete [] GlobalCopyBlockIndices;
  delete [] GlobalCopyColDims;
  delete [] GlobalCopyLDAs;
  for (i=0; i<MaxNumBlockEntries; i++) delete [] GlobalCopyValuesPointers[i];
  delete [] GlobalCopyValuesPointers;
  delete [] MyView1ValuesPointers;
  if (verbose) cout << "\n\nRows sorted check OK" << endl<< endl;
  
  return ierr;
}

//=============================================================================
int CompareValues(double * A, int LDA, int NumRowsA, int NumColsA, 
		  double * B, int LDB, int NumRowsB, int NumColsB) {
  
  int i, j, ierr = 0, forierr = 0;
  double * ptr1 = B;
  double * ptr2;
  
  if (NumRowsA!=NumRowsB) EPETRA_TEST_ERR(-2,ierr);
  if (NumColsA!=NumColsB) EPETRA_TEST_ERR(-3,ierr);
 

  forierr = 0;
  for (j=0; j<NumColsA; j++) {
    ptr1 = B + j*LDB;
    ptr2 = A + j*LDA;
    for (i=0; i<NumRowsA; i++) forierr += (*ptr1++ != *ptr2++);
  }
  EPETRA_TEST_ERR(forierr,ierr);
  return ierr;
}

int checkMergeRedundantEntries(Epetra_Comm& comm, bool verbose)
{
  int numProcs = comm.NumProc();
  int localProc = comm.MyPID();

  int myFirstRow = localProc*3;
  int myLastRow = myFirstRow+2;
  int numMyRows = myLastRow - myFirstRow + 1;
  int numGlobalRows = numProcs*numMyRows;
  int i,j, ierr;

  //We'll set up a matrix which is globally block-diagonal, i.e., on each
  //processor the list of columns == list of rows.
  //Set up a list of column-indices which is twice as long as it should be,
  //and its contents will be the list of local rows, repeated twice.
  int numCols = 2*numMyRows;
  int* myCols = new int[numCols];

  int col = myFirstRow;
  for(i=0; i<numCols; ++i) {
    myCols[i] = col++;
    if (col > myLastRow) col = myFirstRow;
  }

  int elemSize = 2;
  int indexBase = 0;

  Epetra_BlockMap map(numGlobalRows, numMyRows,
		      elemSize, indexBase, comm);

  Epetra_VbrMatrix A(Copy, map, numCols);

  Epetra_MultiVector x(map, 1), y(map, 1);
  x.PutScalar(1.0);

  Epetra_MultiVector x3(map, 3), y3(map, 3);
  x.PutScalar(1.0);

  double* coef = new double[elemSize*elemSize];
  for(i=0; i<elemSize*elemSize; ++i) {
    coef[i] = 0.5;
  }

  //we're going to insert each row twice, with coef values of 0.5. So after
  //FillComplete, which internally calls MergeRedundantEntries, the
  //matrix should contain 1.0 in each entry.

  for(i=myFirstRow; i<=myLastRow; ++i) {
    EPETRA_TEST_ERR( A.BeginInsertGlobalValues(i, numCols, myCols), ierr);

    for(j=0; j<numCols; ++j) {
      EPETRA_TEST_ERR( A.SubmitBlockEntry(coef, elemSize,
					  elemSize, elemSize), ierr);
    }

    EPETRA_TEST_ERR( A.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( A.FillComplete(), ierr);

  delete [] coef;

  if (verbose) cout << "Multiply x"<<endl;
  EPETRA_TEST_ERR( A.Multiply(false, x, y), ierr );
  if (verbose) cout << y <<endl;

  //Next we're going to extract pointers-to-block-rows and check values to make
  //sure that the internal method Epetra_VbrMatrix::mergeRedundantEntries()
  //worked correctly. 
  //At the same time, we're also going to create another VbrMatrix which will
  //be a View of the matrix we've already created. This will serve to provide
  //more test coverage of the VbrMatrix code.

  int numBlockEntries = 0;
  int RowDim;
  int** BlockIndices = new int*[numMyRows];
  Epetra_SerialDenseMatrix** Values;
  Epetra_VbrMatrix Aview(View, map, numMyRows);

  for(i=myFirstRow; i<=myLastRow; ++i) {
    BlockIndices[i-myFirstRow] = new int[numCols];
    EPETRA_TEST_ERR( A.ExtractGlobalBlockRowPointers(i, numCols,
						     RowDim, numBlockEntries,
						     BlockIndices[i-myFirstRow],
						     Values), ierr);

    EPETRA_TEST_ERR( Aview.BeginInsertGlobalValues(i, numBlockEntries,
					      BlockIndices[i-myFirstRow]), ierr);

    if (numMyRows != numBlockEntries) return(-1);
    if (RowDim != elemSize) return(-2);
    for(j=0; j<numBlockEntries; ++j) {
      if (Values[j]->A()[0] != 1.0) {
	cout << "Row " << i << " Values["<<j<<"][0]: "<< Values[j][0]
	     << " should be 1.0" << endl;
	return(-3); //comment-out this return to de-activate this test
      }

      EPETRA_TEST_ERR( Aview.SubmitBlockEntry(Values[j]->A(),
					      Values[j]->LDA(),
					      Values[j]->M(),
					      Values[j]->N()), ierr);
    }

    EPETRA_TEST_ERR( Aview.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( Aview.FillComplete(), ierr);

  //So the test appears to have passed for the original matrix A. Now check the
  //values of our second "view" of the matrix, 'Aview'.

  for(i=myFirstRow; i<=myLastRow; ++i) {
    EPETRA_TEST_ERR( Aview.ExtractGlobalBlockRowPointers(i, numMyRows,
							 RowDim, numBlockEntries,
							 BlockIndices[i-myFirstRow],
							 Values), ierr);

    if (numMyRows != numBlockEntries) return(-1);
    if (RowDim != elemSize) return(-2);
    for(j=0; j<numBlockEntries; ++j) {
      if (Values[j]->A()[0] != 1.0) {
	cout << "Aview: Row " << i << " Values["<<j<<"][0]: "<< Values[j][0]
	     << " should be 1.0" << endl;
	return(-3); //comment-out this return to de-activate this test
      }
    }
    
    delete [] BlockIndices[i-myFirstRow];
  }

  if (verbose&&localProc==0) {
    cout << "checkMergeRedundantEntries, A:" << endl;
  }

  if (verbose) {
  A.Print(cout);
  }

  delete [] BlockIndices;
  delete [] myCols;

  return(0);
}

int checkExtractMyRowCopy(Epetra_Comm& comm, bool verbose)
{
  int numProcs = comm.NumProc();
  int localProc = comm.MyPID();

  int myFirstRow = localProc*3;
  int myLastRow = myFirstRow+2;
  int numMyRows = myLastRow - myFirstRow + 1;
  int numGlobalRows = numProcs*numMyRows;
  int i,j, ierr;

  int numCols = numMyRows;
  int* myCols = new int[numCols];

  int col = myFirstRow;
  for(i=0; i<numCols; ++i) {
    myCols[i] = col++;
    if (col > myLastRow) col = myFirstRow;
  }

  int elemSize = 2;
  int indexBase = 0;

  Epetra_BlockMap map(numGlobalRows, numMyRows,
		      elemSize, indexBase, comm);

  Epetra_VbrMatrix A(Copy, map, numCols);

  double* coef = new double[elemSize*elemSize];

  for(i=myFirstRow; i<=myLastRow; ++i) {
    int myPointRow = i*elemSize;

    //The coefficients need to be laid out in column-major order. i.e., the
    //coefficients in a column occur contiguously.
    for(int ii=0; ii<elemSize; ++ii) {
      for(int jj=0; jj<elemSize; ++jj) {
	double val = (myPointRow+ii)*1.0;
	coef[ii+elemSize*jj] = val;
      }
    }

    EPETRA_TEST_ERR( A.BeginInsertGlobalValues(i, numCols, myCols), ierr);

    for(j=0; j<numCols; ++j) {
      EPETRA_TEST_ERR( A.SubmitBlockEntry(coef, elemSize,
					  elemSize, elemSize), ierr);
    }

    EPETRA_TEST_ERR( A.EndSubmitEntries(), ierr);
  }

  EPETRA_TEST_ERR( A.FillComplete(), ierr);

  delete [] coef;
  delete [] myCols;

  Epetra_SerialDenseMatrix** blockEntries;
  int len = elemSize*numCols, checkLen;
  double* values = new double[len];
  int* indices = new int[len];
  int RowDim, numBlockEntries;

  for(i=myFirstRow; i<=myLastRow; ++i) {
    EPETRA_TEST_ERR( A.ExtractGlobalBlockRowPointers(i, numMyRows,
						     RowDim, numBlockEntries,
						     indices,
						     blockEntries), ierr);
    if (numMyRows != numBlockEntries) return(-1);
    if (RowDim != elemSize) return(-2);

    int myPointRow = i*elemSize - myFirstRow*elemSize;
    int ii,jj;
    for(ii=0; ii<elemSize; ++ii) {
      EPETRA_TEST_ERR( A.ExtractMyRowCopy(myPointRow+ii, len,
					  checkLen, values, indices), ierr);
      if (len != checkLen) return(-3);

      double val = (i*elemSize+ii)*1.0;
      double blockvalue = blockEntries[0]->A()[ii];

      for(jj=0; jj<len; ++jj) {
	if (values[jj] != val) return(-4);
	if (values[jj] != blockvalue) return(-5);
      }
    }
  }

  delete [] values;
  delete [] indices;

  return(0);
}
