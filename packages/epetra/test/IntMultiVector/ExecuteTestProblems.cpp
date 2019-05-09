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
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"

int IntMultiVectorTests(const Epetra_BlockMap & Map, int NumVectors, bool verbose)
{
  int ierr = 0;

  Epetra_BLAS BLAS;

  //*******************************************************************
  // Post-construction modification tests
  //*******************************************************************

  if (verbose) cout <<  "\n\nXXXXX Testing Post-construction modification of a multivector"
		    <<endl<<endl;

  int err = 0;

  Epetra_IntMultiVector X(Map, NumVectors);

  // Pick middle range values for GID, LID and Vector Index
  int testGID = Map.NumGlobalElements()/2;
  int testVecIndex = NumVectors/2;

  int GIDSize = 1;
  int LIDOfGID = 0;
  int FirstEntryOfGID = 0;

  if (Map.MyGID(testGID)) {
    LIDOfGID = Map.LID(testGID);
    GIDSize = Map.ElementSize(LIDOfGID);
    FirstEntryOfGID = Map.FirstPointInElement(LIDOfGID);
  }

  // ========================================================================
  // Test int ReplaceGlobalValue (int GlobalRow, int VectorIndex, int OrdinalValue)
  // ========================================================================

  int newGIDValue = 4;
  int locerr = 0;
  locerr = X.ReplaceGlobalValue(testGID, testVecIndex, newGIDValue);

  if (Map.MyGID(testGID)) {
    if (X[testVecIndex][FirstEntryOfGID]!=newGIDValue) err++;
    if (verbose) cout << "X["<<testVecIndex<<"]["<<FirstEntryOfGID<<"] = "
		      <<  X[testVecIndex][FirstEntryOfGID]
		      << " should = " << newGIDValue << endl;
  }
  else
    if (locerr!=0) err++; // Test for GID out of range error (=1)

  // ========================================================================
  // Test int ReplaceGlobalValue (int GlobalRow, intBlockRowOffset, int VectorIndex, int OrdinalValue)
  // ========================================================================
  newGIDValue = 8;
  locerr = X.ReplaceGlobalValue(testGID, GIDSize-1, testVecIndex, newGIDValue);

  if (Map.MyGID(testGID)) {
    if (X[testVecIndex][FirstEntryOfGID+GIDSize-1]!=newGIDValue) err++;
    if (verbose) cout << "X["<<testVecIndex<<"]["<<FirstEntryOfGID+GIDSize-1<<"] = "
		      <<  X[testVecIndex][FirstEntryOfGID+GIDSize-1]
		      << " should = " << newGIDValue << endl;
  }
  else
    if (locerr!=0) err++; // Test for GID out of range error (=1)

  // ========================================================================
  // Test int SumIntoGlobalValue (int GlobalRow, int VectorIndex, int OrdinalValue)
  // ========================================================================

  newGIDValue = 1;
  locerr = X.ReplaceGlobalValue(testGID, testVecIndex, newGIDValue);
  locerr = X.SumIntoGlobalValue(testGID, testVecIndex, newGIDValue);
  if (Map.MyGID(testGID)) {
    if (X[testVecIndex][FirstEntryOfGID]!=(newGIDValue+newGIDValue)) err++;
    if (verbose) cout << "X["<<testVecIndex<<"]["<<FirstEntryOfGID<<"] = "
		      <<  X[testVecIndex][FirstEntryOfGID]
		      << " should = " << newGIDValue+newGIDValue << endl;
  }
  else
    if (locerr!=0) err++; // Test for GID out of range error (=1)

  // ========================================================================
  // Test int SumIntoGlobalValue (int GlobalRow, intBlockRowOffset, int VectorIndex, int OrdinalValue)
  // ========================================================================

  newGIDValue = 1;
  locerr = X.ReplaceGlobalValue(testGID, GIDSize-1, testVecIndex, newGIDValue);
  locerr = X.SumIntoGlobalValue(testGID, GIDSize-1, testVecIndex, newGIDValue);

  if (Map.MyGID(testGID)) {
    if (X[testVecIndex][FirstEntryOfGID+GIDSize-1]!=(newGIDValue+newGIDValue)) err++;
    if (verbose) cout << "X["<<testVecIndex<<"]["<<FirstEntryOfGID+GIDSize-1<<"] = "
		      <<  X[testVecIndex][FirstEntryOfGID+GIDSize-1]
		      << " should = " << newGIDValue+newGIDValue << endl;
  }
  else
    if (locerr!=0) err++; // Test for GID out of range error (=1)

  // ========================================================================
  // Test Local "My" versions of same routine (less complicated)
  // ========================================================================

  // Pick middle range values for LID
  int testLID = Map.NumMyElements()/2;

  int LIDSize = Map.ElementSize(testLID);
  int FirstEntryOfLID = Map.FirstPointInElement(testLID);


  int newLIDValue = 4;
  locerr = X.ReplaceMyValue(testLID, testVecIndex, newLIDValue);

  if (X[testVecIndex][FirstEntryOfLID]!=newLIDValue) err++;
  if (verbose) cout << "X["<<testVecIndex<<"]["<<FirstEntryOfLID<<"] = "
		    <<  X[testVecIndex][FirstEntryOfLID]
		    << " should = " << newLIDValue << endl;

  newLIDValue = 8;
  locerr = X.ReplaceMyValue(testLID, LIDSize-1, testVecIndex, newLIDValue);
  if (X[testVecIndex][FirstEntryOfLID+LIDSize-1]!=newLIDValue) err++;
  if (verbose) cout << "X["<<testVecIndex<<"]["<<FirstEntryOfLID+LIDSize-1<<"] = "
		    <<  X[testVecIndex][FirstEntryOfLID+LIDSize-1]
		    << " should = " << newLIDValue << endl;
  newLIDValue = 1;
  locerr = X.ReplaceMyValue(testLID, testVecIndex, newLIDValue);
  locerr = X.SumIntoMyValue(testLID, testVecIndex, newLIDValue);
  if (X[testVecIndex][FirstEntryOfLID]!=(newLIDValue+newLIDValue)) err++;
  if (verbose) cout << "X["<<testVecIndex<<"]["<<FirstEntryOfLID<<"] = "
		    <<  X[testVecIndex][FirstEntryOfLID]
		    << " should = " << newLIDValue+newLIDValue << endl;
  newLIDValue = 2;
  locerr = X.ReplaceMyValue(testLID, LIDSize-1, testVecIndex, newLIDValue);
  locerr = X.SumIntoMyValue(testLID, LIDSize-1, testVecIndex, newLIDValue);
  if (X[testVecIndex][FirstEntryOfLID+LIDSize-1]!=(newLIDValue+newLIDValue)) err++;
  if (verbose) cout << "X["<<testVecIndex<<"]["<<FirstEntryOfLID+LIDSize-1<<"] = "
		    <<  X[testVecIndex][FirstEntryOfLID+LIDSize-1]
		    << " should = " << newLIDValue+newLIDValue << endl;

  ierr += err;

  return(ierr);
}
