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


#include "Epetra_BLAS.h"
#include "ExecuteTestProblems.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"

int MultiVectorTests(const Epetra_BlockMap & Map, int NumVectors, bool verbose)
{
  (void)NumVectors;
  const Epetra_Comm & Comm = Map.Comm();
  int ierr = 0;
  
  /* get number of processors and the name of this processor */
  
  // int NumProc = Comm.getNumProc();
  int MyPID   = Comm.MyPID();
  
  // Construct FEVector
  
  if (verbose&&MyPID==0) cout << "constructing Epetra_FEVector" << endl;

  Epetra_FEVector A(Map, 1);
 
  //For an extreme test, we'll have each processor sum-in a 1.0 for All
  //global ids.

  int minGID = Map.MinAllGID();
  int numGlobalIDs = Map.NumGlobalElements();

  //For now we're going to have just one point associated with
  //each GID (element).

  int* ptIndices = new int[numGlobalIDs];
  double* ptCoefs = new double[numGlobalIDs];

  Epetra_IntSerialDenseVector epetra_indices(View, ptIndices, numGlobalIDs);
  Epetra_SerialDenseVector epetra_coefs(View, ptCoefs, numGlobalIDs);

  {for(int i=0; i<numGlobalIDs; ++i) {
    ptIndices[i] = minGID+i;
    ptCoefs[i] = 1.0;
  }}

  if (verbose&&MyPID==0) {
    cout << "calling A.SumIntoGlobalValues with " << numGlobalIDs << " values"<<endl;
  }
  EPETRA_TEST_ERR( A.SumIntoGlobalValues(numGlobalIDs, ptIndices, ptCoefs), ierr);

  if (verbose&&MyPID==0) {
    cout << "calling A.SumIntoGlobalValues with " << numGlobalIDs << " values"<<endl;
  }
  EPETRA_TEST_ERR( A.SumIntoGlobalValues(epetra_indices, epetra_coefs), ierr);

  if (verbose&&MyPID==0) {
    cout << "calling A.GlobalAssemble()" << endl;
  }

  EPETRA_TEST_ERR( A.GlobalAssemble(), ierr );

  if (verbose&&MyPID==0) {
  cout << "after globalAssemble"<<endl;
  }
  if (verbose) {
  A.Print(cout);
  }

  //now do a quick test of the copy constructor
  Epetra_FEVector B(A);

  double nrm2a, nrm2b;
  A.Norm2(&nrm2a);
  B.Norm2(&nrm2b);

  if (nrm2a != nrm2b) {
    cerr << "copy-constructor test failed, norm of copy doesn't equal"
         << " norm of original."<<endl;
    return(-1);
  }

  delete [] ptIndices;
  delete [] ptCoefs;

  return(ierr);
}

int fevec0(Epetra_Comm& Comm, bool verbose)
{
  int ierr = 0;
  int NumGlobalRows = 4;
  int indexBase = 0;
  Epetra_Map Map(NumGlobalRows, indexBase, Comm);

  int Numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if (Numprocs != 2) return(0);


  int NumCols = 3;
  int* Indices = new int[NumCols];

  double* Values = new double[NumCols];

// Create vectors

  Epetra_FEVector b(Map, 1);
  Epetra_FEVector x0(Map, 1);

// source terms
  NumCols = 2;

  if(MyPID==0)  // indices corresponding to element 0 on processor 0
  {
    Indices[0] = 0;
    Indices[1] = 3;

    Values[0] = 1./2.;
    Values[1] = 1./2.;

   }
   else
   {
    Indices[0] = 1;
    Indices[1] = 2;

    Values[0] = 0;
    Values[1] = 0;
   }

  EPETRA_TEST_ERR( b.SumIntoGlobalValues(NumCols, Indices, Values),
                   ierr);

  EPETRA_TEST_ERR( b.GlobalAssemble(), ierr);

  if (verbose&&MyPID==0) {
    cout << "b:"<<endl;
  }

  if (verbose) {
  b.Print(cout);
  }

  x0 = b;

  if (verbose&&MyPID==0) {
    cout << "x:"<<endl;
  }

  if (verbose) {
  x0.Print(cout);
  }

  delete [] Values;
  delete [] Indices;

  return(0);
}

int fevec1(Epetra_Comm& Comm, bool verbose)
{
  int Numprocs = Comm.NumProc();

  if (Numprocs != 2) return(0);
  int MyPID   = Comm.MyPID();

  int ierr = 0;
  int NumGlobalRows = 6;
  const int NumVectors = 4;
  int indexBase = 0;
  Epetra_Map Map(NumGlobalRows, indexBase, Comm);

  const int Num = 4;
  int Indices[Num];
 
  double Values[Num];
 
// Create vectors

  Epetra_FEVector b(Map, NumVectors);
  Epetra_FEVector x0(Map, NumVectors);
 
// source terms
 
  if(MyPID==0)  // indices corresponding to element 0 on processor 0
  {
    Indices[0] = 0;
    Indices[1] = 1;
    Indices[2] = 4;
    Indices[3] = 5;
 
    Values[0] = 1./2.;
    Values[1] = 1./2.;
    Values[2] = 1./2.;
    Values[3] = 1./2.;

   }
   else
   {
    Indices[0] = 1;
    Indices[1] = 2;
    Indices[2] = 3;
    Indices[3] = 4;
 
    Values[0] = 0;
    Values[1] = 0;
    Values[2] = 0;
    Values[3] = 0;
   }

  for(int i=0; i<NumVectors; ++i) {
    EPETRA_TEST_ERR( b.SumIntoGlobalValues(Num, Indices, Values, i),
		   ierr);
  }

  EPETRA_TEST_ERR( b.GlobalAssemble(), ierr);

  double nrm2[NumVectors];

  b.Norm2(nrm2);

  for(int i=1; i<NumVectors; ++i) {
    if (fabs(nrm2[i]-nrm2[0]) > 1.e-12) {
      EPETRA_TEST_ERR(-1, ierr);
      return(-1);
    }
  }


  //now sum-in again, to make sure the previous call to GlobalAssemble
  //didn't do something nasty to internal non-local data structures.
  //(This is a specific case that has bitten me. Hence this test...)
  for(int i=0; i<NumVectors; ++i) {
    EPETRA_TEST_ERR( b.SumIntoGlobalValues(Num, Indices, Values, i),
                   ierr);
  }

  //and now GlobalAssemble again...
  EPETRA_TEST_ERR( b.GlobalAssemble(), ierr);


  if (verbose&&MyPID==0) {
    cout << "b:"<<endl;
  }

  if (verbose) {
    b.Print(cout);
  }

  x0 = b;

  if (verbose&&MyPID==0) {
    cout << "x:"<<endl;
  }

  if (verbose) {
    x0.Print(cout);
  }

  return(0);
}

int fevec2(Epetra_Comm& Comm, bool verbose)
{
  int ierr = 0;
  int NumGlobalElems = 4;
  int elemSize = 3;
  int indexBase = 0;
  Epetra_BlockMap Map(NumGlobalElems, elemSize, indexBase, Comm);

  int Numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if (Numprocs != 2) return(0);

  int NumCols = 3;
  int* Indices = new int[NumCols];
  int* numValuesPerID = new int[NumCols];
  for(int i=0; i<NumCols; ++i) {
    numValuesPerID[i] = elemSize;
  }
 
  double* Values = new double[NumCols*elemSize];
 
// Create vectors

  Epetra_FEVector b(Map, 1);
  Epetra_FEVector x0(Map, 1);
 
// source terms
  NumCols = 2;
 
  if(MyPID==0)  // indices corresponding to element 0 on processor 0
  {
    Indices[0] = 0;
    Indices[1] = 3;
 
    Values[0] = 1./2.;
    Values[1] = 1./2.;
    Values[2] = 1./2.;
    Values[3] = 1./2.;
    Values[4] = 1./2.;
    Values[5] = 1./2.;

   }
   else
   {
    Indices[0] = 1;
    Indices[1] = 2;
 
    Values[0] = 0;
    Values[1] = 0;
    Values[2] = 0;
    Values[3] = 0;
    Values[4] = 0;
    Values[5] = 0;
   }

  EPETRA_TEST_ERR( b.SumIntoGlobalValues(NumCols, Indices,
					 numValuesPerID, Values),
		   ierr);

  EPETRA_TEST_ERR( b.GlobalAssemble(), ierr);

  if (verbose&&MyPID==0) {
    cout << "b:"<<endl;
  }

  if (verbose) {
  b.Print(cout);
  }

  x0 = b;

  if (verbose&&MyPID==0) {
    cout << "x:"<<endl;
  }

  if (verbose) {
  x0.Print(cout);
  }

  delete [] Values;
  delete [] Indices;
  delete [] numValuesPerID;

  return(0);
}

int fevec3(Epetra_Comm& Comm, bool verbose)
{
  int ierr = 0;
  int NumGlobalElems = 4;
  int elemSize = 40;
  int indexBase = 0;
  Epetra_BlockMap Map(NumGlobalElems, elemSize, indexBase, Comm);

  int Numprocs = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if (Numprocs != 2) return(0);

  int NumCols = 3;
  int* Indices = new int[NumCols];
  int* numValuesPerID = new int[NumCols];
  for(int i=0; i<NumCols; ++i) {
    numValuesPerID[i] = elemSize;
  }
 
  double* Values = new double[NumCols*elemSize];
 
// Create vectors

  Epetra_FEVector b(Map, 1);
  Epetra_FEVector x0(Map, 1);
 
// source terms
  NumCols = 2;
 
  if(MyPID==0)  // indices corresponding to element 0 on processor 0
  {
    Indices[0] = 0;
    Indices[1] = 3;
 
    for(int ii=0; ii<NumCols*elemSize; ++ii) {
      Values[ii] = 1./2.;
    }

  }
  else
  {
    Indices[0] = 1;
    Indices[1] = 2;
 
    for(int ii=0; ii<NumCols*elemSize; ++ii) {
      Values[ii] = 0.;
    }

  }

  EPETRA_TEST_ERR( b.SumIntoGlobalValues(NumCols, Indices,
					 numValuesPerID, Values),
		   ierr);

  EPETRA_TEST_ERR( b.GlobalAssemble(), ierr);

  if (verbose&&MyPID==0) {
    cout << "b:"<<endl;
  }

  if (verbose) {
  b.Print(cout);
  }

  x0 = b;

  if (verbose&&MyPID==0) {
    cout << "x:"<<endl;
  }

  if (verbose) {
  x0.Print(cout);
  }

  delete [] Values;
  delete [] Indices;
  delete [] numValuesPerID;

  return(0);
}

int fevec4(Epetra_Comm& Comm, bool verbose)
{
  int NumElements = 4;
  Epetra_Map     Map(NumElements, 0, Comm);
  Epetra_FEVector x1(Map);
  const double value = 1.;
  x1.PutScalar (value);
				// replace one element by itself. processor 0
				// does not own this element
  const int GID = 3;
  x1.ReplaceGlobalValues(1, &GID, &value);
  x1.GlobalAssemble (Insert);

  if (Map.MyGID(3)) {
    //insist that the value for GID==3 is 1:
    if (std::abs(x1.Values()[Map.LID(3)] - 1) > 1.e-9) return -1;
  }

  std::cout << x1;

  Comm.Barrier();

				// re-apply GlobalAssemble. Nothing should
				// happen
  x1.GlobalAssemble (Insert);
  std::cout << x1;
  if (Map.MyGID(3)) {
    //insist that the value for GID==3 is 1:
    if (std::abs(x1.Values()[Map.LID(3)] - 1) > 1.e-9) return -1;
  }

  return 0;
}

int fevec5(Epetra_Comm& Comm, bool verbose)
{
  int NumElements = 4;
  Epetra_Map     Map(NumElements, 0, Comm);
  Epetra_FEVector x1(Map);
  x1.PutScalar (0);

        // let all processors set global entry 0 to 1
  const int GID = 0;
  const double value = 1;
  x1.ReplaceGlobalValues(1, &GID, &value);
  x1.GlobalAssemble (Insert);
  if (Comm.MyPID()==0)
    std::cout << "Entry " << GID << " after construct & set: " 
        << x1[0][0] << std::endl;

        // copy vector
  Epetra_FEVector x2 (x1);

  x2.PutScalar(0);

        // re-apply 1 to the vector, but only on the
        // owning processor. should be enough to set
        // the value (as non-local data in x1 should
        // have been eliminated after calling
        // GlobalAssemble).
  if (Comm.MyPID()==0)
    x2.ReplaceGlobalValues(1, &GID, &value);
  x2.GlobalAssemble (Insert);

  if (Comm.MyPID()==0)
    std::cout << "Entry " << GID << " after copy & set:      " 
        << x2[0][0] << std::endl;

  return 0;
}

int fevec6(Epetra_Comm& Comm, bool verbose)
{
  int NumElements = 4;
  Epetra_Map     Map(NumElements, 0, Comm);
  Epetra_FEVector x1(Map);
  x1.PutScalar (0);

        // let all processors set global entry 0 to 1
  const int GID = 0;
  const double value = 1;
  x1.ReplaceGlobalValues(1, &GID, &value);
  x1.GlobalAssemble (Insert);
  if (Comm.MyPID()==0)
    std::cout << "Entry " << GID << " after construct & set: " 
        << x1[0][0] << std::endl;

  x1.PutScalar(0);

        // re-apply 1 to the vector, but only on the
        // owning processor. should be enough to set
        // the value (as non-local data in x1 should
        // have been eliminated after calling
        // GlobalAssemble).
  if (Comm.MyPID()==0)
    x1.ReplaceGlobalValues(1, &GID, &value);
  x1.GlobalAssemble (Insert);

  if (Comm.MyPID()==0)
    std::cout << "Entry " << GID << " after PutScalar & set:      " 
        << x1[0][0] << std::endl;

  return 0;
}

