
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

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
//=============================================================================
Epetra_Vector::Epetra_Vector(const Epetra_BlockMap& map, bool zeroOut)
  : Epetra_MultiVector(map,1,zeroOut) // Vector is just special case of MultiVector
{
  SetLabel("Epetra::Vector");
}
//=============================================================================
Epetra_Vector::Epetra_Vector(const Epetra_Vector& Source)
  : Epetra_MultiVector(Source) // Vector is just special case of MultiVector
{
}
//=============================================================================
Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_BlockMap& map, double *V)
  : Epetra_MultiVector(CV, map, V, map.NumMyPoints(), 1) // Vector is just special case of MultiVector
{
  SetLabel("Epetra::Vector");
}
//=============================================================================
Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_MultiVector& Source, int Index)
  : Epetra_MultiVector(CV, Source, Index, 1) // Vector is just special case of MultiVector
{
  SetLabel("Epetra::Vector");
}
//=========================================================================
Epetra_Vector::~Epetra_Vector(){}

//=============================================================================
int Epetra_Vector::ExtractCopy(double *V) const {
  return(Epetra_MultiVector::ExtractCopy(V, 1));
}

//=============================================================================
int Epetra_Vector::ExtractView(double **V) const {
  int junk;
  return(Epetra_MultiVector::ExtractView(V, &junk));
}

/*
//=========================================================================
double& Epetra_Vector::operator [] (int Index)  {

   return(Values_[Index]);
}

//=========================================================================
const double& Epetra_Vector::operator [] (int Index) const  {

   return(Values_[Index]);
}
*/

//=========================================================================
int Epetra_Vector::ReplaceGlobalValues(int NumEntries, const double * values, const int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, 0, values, Indices, true, false));
  return(0);
}
//=========================================================================
int Epetra_Vector::ReplaceMyValues(int NumEntries, const double * values, const int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, 0, values, Indices, false, false));
  return(0);
}
//=========================================================================
int Epetra_Vector::SumIntoGlobalValues(int NumEntries, const double * values, const int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, 0, values, Indices, true, true));
  return(0);
}
//=========================================================================
int Epetra_Vector::SumIntoMyValues(int NumEntries, const double * values, const int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, 0, values, Indices, false, true));
  return(0);
}
//=========================================================================
int Epetra_Vector::ReplaceGlobalValues(int NumEntries, int BlockOffset, const double * values, const int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, BlockOffset, values, Indices, true, false));
  return(0);
}
//=========================================================================
int Epetra_Vector::ReplaceMyValues(int NumEntries, int BlockOffset, const double * values, const int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, BlockOffset, values, Indices, false, false));
  return(0);
}
//=========================================================================
int Epetra_Vector::SumIntoGlobalValues(int NumEntries, int BlockOffset, const double * values, const int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, BlockOffset, values, Indices, true, true));
  return(0);
}
//=========================================================================
int Epetra_Vector::SumIntoMyValues(int NumEntries, int BlockOffset, const double * values, const int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, BlockOffset, values, Indices, false, true));
  return(0);
}
//=========================================================================
int Epetra_Vector::ChangeValues(int NumEntries, int BlockOffset, const double * values, const int * Indices,
				bool IndicesGlobal, bool SumInto) {

  int cur_index;
  int ierr = 0;
  if (BlockOffset<0) EPETRA_CHK_ERR(-1); // Offset is out-of-range

  for (int i=0; i<NumEntries; i++) {
    if (IndicesGlobal) 
      cur_index = Map().LID(Indices[i]);
    else
      cur_index = Indices[i];
    
    if (Map().MyLID(cur_index)) {
      if (BlockOffset>=Map().ElementSize(cur_index)) EPETRA_CHK_ERR(-1); // Offset is out-of-range
      int entry = Map().FirstPointInElement(cur_index);

      if (SumInto)
	Values_[entry+BlockOffset] += values[i];
      else
	Values_[entry+BlockOffset] = values[i];
    }
    else ierr = 1;
  }

  EPETRA_CHK_ERR(ierr);
  return(0);
}
