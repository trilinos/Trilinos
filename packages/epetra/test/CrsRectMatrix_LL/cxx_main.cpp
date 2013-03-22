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


#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Flops.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "../epetra_test_err.h"
#include "Epetra_Version.h"
 
int main(int argc, char *argv[])
{
  int ierr = 0, i, j, forierr = 0;
  int ntrials = 1;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;



  // char tmp; if (Comm.MyPID()==0) { cout << "Press any key to continue..."<< endl;  cin >> tmp;} Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if (verbose && Comm.MyPID()==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && Comm.MyPID()!=0) verbose = false;

  int NumMyEquations = 10000;

  long long NumGlobalEquations = NumMyEquations*NumProc;
  long long NumGlobalVariables = 2 * NumGlobalEquations+1;

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_Map RowMap(NumGlobalEquations, 0LL, Comm);
  Epetra_Map XMap(NumGlobalVariables, 0LL, Comm);
  Epetra_Map& YMap = RowMap;
  
  // Get update list and number of local equations from newly created Map
  long long * MyGlobalElements = new long long[RowMap.NumMyElements()];
  RowMap.MyGlobalElements(MyGlobalElements);

  // Get update list and number of local equations from newly created XMap
  long long * XGlobalElements = new long long[XMap.NumMyElements()];
  XMap.MyGlobalElements(XGlobalElements);

  // Get update list and number of local variables from newly created YMap
  long long * YGlobalElements = new long long[YMap.NumMyElements()];
  YMap.MyGlobalElements(YGlobalElements);

  // We need vectors to compute:
  // X = A^T*Y
  // AATY = A*A^T*Y = A*X
  //  and 
  // BY = B*Y

  Epetra_Vector Y(YMap);
  Epetra_Vector X(XMap);
  Epetra_Vector AATY(YMap);
  Epetra_Vector BY(YMap);


  // Fill Y Vector
  Y.Random();
  //Y.PutScalar(1.0);

  // To create A^T explicitly we need an assembly map that is two elements longer than
  // the XMap, because each processor will be making contributions to two rows beyond what
  // it will own.
  long long ATAssemblyNumMyElements = 2*MyGlobalElements[NumMyEquations-1] + 2 - 2*MyGlobalElements[0] + 1;
  long long * ATAssemblyGlobalElements = new long long[ATAssemblyNumMyElements];

  for (i=0; i<ATAssemblyNumMyElements; i++) ATAssemblyGlobalElements[i] = 2*MyGlobalElements[0] + i;
  Epetra_Map ATAssemblyMap((long long)-1, ATAssemblyNumMyElements, ATAssemblyGlobalElements, 0LL, Comm);

  // Create a Epetra_Matrix with the values of A
  // A is a simple 1D weighted average operator that mimics a restriction operator
  // that might be found in a multigrid code.
  // Also create A^T explicitly

  Epetra_CrsMatrix A(Copy, RowMap, 3);
  Epetra_CrsMatrix ATAssembly(Copy, ATAssemblyMap, 2);
  Epetra_CrsMatrix AT(Copy, XMap, 2);
  
  //cout << "ATAssemblyMap = "<< endl<< ATAssemblyMap << endl
  //     << "XMap = " << endl << XMap << endl
  //     << "RowMap = " << endl << RowMap << endl;
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[3];
  long long *Indices = new long long[3];
  int NumEntries;
  /*
  Values[0] = 0.25;
  Values[1] = 0.5;
  Values[2] = 0.25;
  */
  Values[0] = 0.5;
  Values[1] = 0.25;
  Values[2] = 0.25;
  forierr = 0;
  for (i=0; i<NumMyEquations; i++)
    {
  /*
      Indices[0] = 2*MyGlobalElements[i];
      Indices[1] = 2*MyGlobalElements[i]+1;
      Indices[2] = 2*MyGlobalElements[i]+2;
   */
      Indices[0] = 2*MyGlobalElements[i]+1;
      Indices[1] = 2*MyGlobalElements[i]+2;
      Indices[2] = 2*MyGlobalElements[i];
      NumEntries = 3;
      forierr += !(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
      for (j=0; j<3; j++)
	forierr += !(ATAssembly.InsertGlobalValues(Indices[j],1, &(Values[j]), &(MyGlobalElements[i]))>=0);
    }
  EPETRA_TEST_ERR(forierr,ierr);


  EPETRA_TEST_ERR(!(ATAssembly.FillComplete()==0),ierr);
  // Gather AT values from ATAssembly matrix
  Epetra_Export Exporter(ATAssemblyMap, XMap);
  EPETRA_TEST_ERR(!(AT.Export(ATAssembly, Exporter, Add)==0),ierr);

  // Finish up
  EPETRA_TEST_ERR(!(A.FillComplete(XMap, YMap)==0),ierr);
  EPETRA_TEST_ERR(!(AT.FillComplete(YMap, XMap)==0),ierr);


  if (verbose1 && NumGlobalEquations<20) { 
    if (verbose) cout << "\n\n Matrix A\n" << endl;
    cout << A << endl;
    if (verbose) cout << " \n\n Matrix A Transpose\n" << endl;
    cout <<  AT << endl;
  }


  // Create a Epetra_Matrix containing B = A*A^T.
  // This matrix will be a square tridiagonal matrix.  We will use it to compare the results
  // of A*(A^T*X) using two methods: (1) with two calls to Multiply using A^T and then A and
  // (2) using B directly.

  Epetra_CrsMatrix B(Copy, RowMap, 3);

  Values[0] = 1.0/16.0;
  Values[1] = 3.0/8.0;
  Values[2] = 1.0/16.0;
  int Valstart;
  forierr = 0;
  for (i=0; i<NumMyEquations; i++)
    {
      if (MyGlobalElements[i] == 0) {
      Indices[0] = MyGlobalElements[i];
      Indices[1] = MyGlobalElements[i]+1;
      NumEntries = 2;
      Valstart = 1;
      }
      else {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i];
	Indices[2] = MyGlobalElements[i]+1;
	NumEntries = 3;
	Valstart = 0;
      }
      if (MyGlobalElements[i] == NumGlobalEquations-1) NumEntries--;
      forierr += !(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values+Valstart, Indices)==0);
    }
  EPETRA_TEST_ERR(forierr,ierr);

  // Finish up
  EPETRA_TEST_ERR(!(B.FillComplete()==0),ierr);
  if (verbose && NumGlobalEquations<20) cout << "\n\nMatrix B \n" << endl;
  if (verbose1 && NumGlobalEquations<20) cout << B << endl;


  Epetra_Flops counter;
  A.SetFlopCounter(counter);
  B.SetFlopCounter(A);
  Epetra_Time timer(Comm);
  for (i=0; i<ntrials; i++) {
    EPETRA_TEST_ERR(!(B.Multiply(false, Y, BY)==0),ierr); // Compute BY = B*Y
  }
  double elapsed_time = timer.ElapsedTime();
  double total_flops = B.Flops();
  counter.ResetFlops();
  
  double MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for B*Y = " << MFLOPs << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = B*Y \n";
  if (verbose1 && NumGlobalEquations<20) cout << BY << endl;
 

  /////////////////////////////////////////////////////////////////////////////////////////////////

  timer.ResetStartTime();
  for (i=0; i<ntrials; i++) {
    EPETRA_TEST_ERR(!(A.Multiply(true, Y, X)==0),ierr); // Compute X = A^T*Y
  }
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops();
  counter.ResetFlops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for A^T*Y using A and trans=true = " << MFLOPs << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = AT*Y \n";
  if (verbose1 && NumGlobalEquations<20) cout << X << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  timer.ResetStartTime();
  EPETRA_TEST_ERR(!(A.Multiply(false, X, AATY)==0),ierr); // Compute AATY = A*X
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  counter.ResetFlops();
  Epetra_Vector resid(YMap);
  resid.Update(1.0, BY, -1.0, AATY, 0.0);
  double residual;
  resid.Norm2(&residual);

  if (verbose) cout << "\n\nTotal MFLOPs for A*X using A and trans=false = " << MFLOPs << endl<< endl;
  if (verbose) cout << "Residual = " << residual << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = A*ATY \n";
  if (verbose1 && NumGlobalEquations<20) cout << AATY << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  AT.SetFlopCounter(counter);
  timer.ResetStartTime();
  for (i=0; i<ntrials; i++) {
    EPETRA_TEST_ERR(!(AT.Multiply(false, Y, X)==0),ierr); // Compute X = A^T*Y
  }
  elapsed_time = timer.ElapsedTime();
  total_flops = AT.Flops();
  counter.ResetFlops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for A^T*Y using AT and trans=false = " << MFLOPs << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = AT*Y \n";
  if (verbose1 && NumGlobalEquations<20) cout << X << endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  timer.ResetStartTime();
  for (i=0; i<ntrials; i++) {
    EPETRA_TEST_ERR(!(AT.Multiply(true, X, AATY)==0),ierr); // Compute AATY = A*X
  }
  elapsed_time = timer.ElapsedTime();
  total_flops = AT.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  counter.ResetFlops();
  resid.Update(1.0, BY, -1.0, AATY, 0.0);
  resid.Norm2(&residual);

  if (verbose) cout << "\n\nTotal MFLOPs for A*X using AT and trans=true = " << MFLOPs << endl<< endl;
  if (verbose) cout << "Residual = " << residual << endl<< endl;
  if (verbose && NumGlobalEquations<20) cout << "\n\nVector Z = A*ATY \n";
  if (verbose1 && NumGlobalEquations<20) cout <<AATY << endl;


  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Now test use of Epetra_LocalMap vectors: First test case of local replicated range vector

  {
    Epetra_CrsMatrix AL(Copy, RowMap, 3);
    for (i=0; i<NumMyEquations; i++)
      {
	forierr += !(A.ExtractGlobalRowCopy(MyGlobalElements[i], 3, NumEntries, Values, Indices)==0);
	forierr += !(AL.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
      }
    EPETRA_TEST_ERR(forierr,ierr);

    Epetra_LocalMap YLMap(NumGlobalEquations, 0LL, Comm);
    EPETRA_TEST_ERR(!(AL.FillComplete(XMap, YLMap)==0),ierr);
    AL.SetFlopCounter(A);
    Epetra_Vector YL(YLMap);
    Epetra_Vector ALX(YLMap);
  
    timer.ResetStartTime();
    for (i=0; i<ntrials; i++) {
      EPETRA_TEST_ERR(!(A.Multiply(false, X, Y)==0),ierr); // Compute Y= A*X
    }
    elapsed_time = timer.ElapsedTime();
    total_flops = A.Flops();
    counter.ResetFlops();
    MFLOPs = total_flops/elapsed_time/1000000.0;



    if (verbose) cout << "\n\nTotal MFLOPs for Y=A*X using global distributed Y = " << MFLOPs << endl<< endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector Y = A*X using distributed Y \n";
    if (verbose1 && NumGlobalEquations<20) cout << Y << endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nA using dist Y range map\n";
    if (verbose1 && NumGlobalEquations<20) cout << A << endl;

    timer.ResetStartTime();
    for (i=0; i<ntrials; i++) {
      EPETRA_TEST_ERR(!(AL.Multiply(false, X, ALX)==0),ierr); // Compute YL= A*X
    }
    elapsed_time = timer.ElapsedTime();
    total_flops = AL.Flops();
    counter.ResetFlops();
    MFLOPs = total_flops/elapsed_time/1000000.0;

    if (verbose) cout << "\n\nTotal MFLOPs for Y=A*X using Local replicated Y = " << MFLOPs << endl<< endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector YL = AL*X using local replicated Y \n";
    if (verbose1 && NumGlobalEquations<20) cout << ALX << endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nA using local Y range map\n";
    if (verbose1 && NumGlobalEquations<20) cout << AL << endl;

    // Now gather Y values from the distributed Y and compare them to the local replicated Y values
    Epetra_Import g2limporter(YLMap, YMap); // Import from distributed Y map to local Y map
    EPETRA_TEST_ERR(!(YL.Import(Y, g2limporter, Insert)==0),ierr);
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector YL = imported from distributed Y \n";
    if (verbose1 && NumGlobalEquations<20) cout << YL << endl;
    EPETRA_TEST_ERR(!(YL.Update(-1.0, ALX, 1.0)==0),ierr);
    EPETRA_TEST_ERR(!(YL.Norm2(&residual)==0),ierr);
    if (verbose) cout << "Residual = " << residual << endl<< endl;
			

    // 
    // Multiply by transpose
    //

    timer.ResetStartTime();
    for (i=0; i<ntrials; i++) {
      EPETRA_TEST_ERR(!(A.Multiply(true, Y, X)==0),ierr); // Compute X = A^TY
    }
    elapsed_time = timer.ElapsedTime();
    total_flops = A.Flops();
    counter.ResetFlops();
    MFLOPs = total_flops/elapsed_time/1000000.0;



    if (verbose) cout << "\n\nTotal MFLOPs for X=A^TY using global distributed Y = " << MFLOPs << endl<< endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector X using distributed Y \n";
    if (verbose1 && NumGlobalEquations<20) cout << Y << endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nA using dist Y range map\n";
    if (verbose1 && NumGlobalEquations<20) cout << A << endl;

    Epetra_Vector X1(XMap);

    timer.ResetStartTime();
    for (i=0; i<ntrials; i++) {
      EPETRA_TEST_ERR(!(AL.Multiply(true, ALX, X1)==0),ierr); // Compute X1 = AL^T*Y
    }
    elapsed_time = timer.ElapsedTime();
    total_flops = AL.Flops();
    counter.ResetFlops();
    MFLOPs = total_flops/elapsed_time/1000000.0;

    if (verbose) cout << "\n\nTotal MFLOPs for X1=A^T*Y using Local replicated Y = " << MFLOPs << endl<< endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector X1 using local replicated Y \n";
    if (verbose1 && NumGlobalEquations<20) cout << X1 << endl;

    EPETRA_TEST_ERR(!(X1.Update(-1.0, X, 1.0)==0),ierr);
    EPETRA_TEST_ERR(!(X1.Norm2(&residual)==0),ierr);
    if (verbose) cout << "Residual = " << residual << endl<< endl;
  }
			
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Finally test use of Epetra_LocalMap vectors using local replicated domain vector

  {
    Epetra_CrsMatrix AL(Copy, RowMap, 3);
    for (i=0; i<NumMyEquations; i++)
      {
	forierr += !(A.ExtractGlobalRowCopy(MyGlobalElements[i], 3, NumEntries, Values, Indices)==0);
	forierr += !(AL.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
      }
    EPETRA_TEST_ERR(forierr,ierr);

    Epetra_LocalMap XLMap(NumGlobalVariables, 0LL, Comm);
    EPETRA_TEST_ERR(!(AL.FillComplete(XLMap, YMap)==0),ierr);
    AL.SetFlopCounter(A);
    Epetra_Vector XL(XLMap);
    Epetra_Vector ALX(XLMap);
  
    timer.ResetStartTime();
    for (i=0; i<ntrials; i++) {
      EPETRA_TEST_ERR(!(A.Multiply(false, X, Y)==0),ierr); // Compute Y= A*X
    }
    elapsed_time = timer.ElapsedTime();
    total_flops = A.Flops();
    counter.ResetFlops();
    MFLOPs = total_flops/elapsed_time/1000000.0;



    if (verbose) cout << "\n\nTotal MFLOPs for Y=A*X using global distributed X = " << MFLOPs << endl<< endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector Y = A*X using distributed X \n";
    if (verbose1 && NumGlobalEquations<20) cout << Y << endl;
    //if (verbose && NumGlobalEquations<20) cout << "\n\nA using dist X range map\n";
    //if (verbose1 && NumGlobalEquations<20) cout << A << endl;

    // Now gather X values from the distributed X 
    Epetra_Import g2limporter(XLMap, XMap); // Import from distributed X map to local X map
    EPETRA_TEST_ERR(!(XL.Import(X, g2limporter, Insert)==0),ierr);
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector XL = imported from distributed X \n";
    if (verbose1 && NumGlobalEquations<20) cout << XL << endl;
    Epetra_Vector Y1(Y);
    timer.ResetStartTime();
    for (i=0; i<ntrials; i++) {
      EPETRA_TEST_ERR(!(AL.Multiply(false, XL, Y1)==0),ierr); // Compute Y1= AL*XL
    }
    elapsed_time = timer.ElapsedTime();
    total_flops = AL.Flops();
    counter.ResetFlops();
    MFLOPs = total_flops/elapsed_time/1000000.0;

    if (verbose) cout << "\n\nTotal MFLOPs for Y1=A*XL using Local replicated X = " << MFLOPs << endl<< endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector Y1 = AL*XL using local replicated X \n";
    if (verbose1 && NumGlobalEquations<20) cout << Y1 << endl;
    //if (verbose && NumGlobalEquations<20) cout << "\n\nA using local X domain map\n";
    //if (verbose1 && NumGlobalEquations<20) cout << AL << endl;

    EPETRA_TEST_ERR(!(Y1.Update(-1.0, Y, 1.0)==0),ierr);
    EPETRA_TEST_ERR(!(Y1.Norm2(&residual)==0),ierr);
    if (verbose) cout << "Residual = " << residual << endl<< endl;
			

    // 
    // Multiply by transpose
    //

    timer.ResetStartTime();
    for (i=0; i<ntrials; i++) {
      EPETRA_TEST_ERR(!(A.Multiply(true, Y, X)==0),ierr); // Compute X = A^TY
    }
    elapsed_time = timer.ElapsedTime();
    total_flops = A.Flops();
    counter.ResetFlops();
    MFLOPs = total_flops/elapsed_time/1000000.0;



    if (verbose) cout << "\n\nTotal MFLOPs for X=A^TY using global distributed X = " << MFLOPs << endl<< endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector X using distributed X \n";
    if (verbose1 && NumGlobalEquations<20) cout << X << endl;
    //if (verbose && NumGlobalEquations<20) cout << "\n\nA using dist X domain map\n";
    //if (verbose1 && NumGlobalEquations<20) cout << A << endl;

    timer.ResetStartTime();
    for (i=0; i<ntrials; i++) {
      EPETRA_TEST_ERR(!(AL.Multiply(true, Y, XL)==0),ierr); // Compute XL = AL^T*Y
    }
    elapsed_time = timer.ElapsedTime();
    total_flops = AL.Flops();
    counter.ResetFlops();
    MFLOPs = total_flops/elapsed_time/1000000.0;

    if (verbose) cout << "\n\nTotal MFLOPs for XL=A^T*Y1 using Local replicated X = " << MFLOPs << endl<< endl;
    if (verbose && NumGlobalEquations<20) cout << "\n\nVector XL using local replicated X \n";
    if (verbose1 && NumGlobalEquations<20) cout << XL << endl;

    Epetra_Vector XL1(XLMap);
    EPETRA_TEST_ERR(!(XL1.Import(X, g2limporter, Insert)==0),ierr);
    EPETRA_TEST_ERR(!(XL1.Update(-1.0, XL, 1.0)==0),ierr);
    EPETRA_TEST_ERR(!(XL1.Norm2(&residual)==0),ierr);
    if (verbose) cout << "Residual = " << residual << endl<< endl;
  }			
  // Release all objects
  delete [] Values;
  delete [] Indices;
  delete [] MyGlobalElements;
  delete [] XGlobalElements;
  delete [] YGlobalElements;
  delete [] ATAssemblyGlobalElements;


#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

