
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

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
//=============================================================================
Epetra_Vector::Epetra_Vector(const Epetra_BlockMap& Map, bool zeroOut)
  : Epetra_MultiVector(Map,1,zeroOut) // Vector is just special case of MultiVector
{
  SetLabel("Epetra::Vector");
}
//=============================================================================
Epetra_Vector::Epetra_Vector(const Epetra_Vector& Source)
  : Epetra_MultiVector(Source) // Vector is just special case of MultiVector
{
}
//=============================================================================
Epetra_Vector::Epetra_Vector(Epetra_DataAccess CV, const Epetra_BlockMap& Map, double *V)
  : Epetra_MultiVector(CV, Map, V, Map.NumMyPoints(), 1) // Vector is just special case of MultiVector
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
int Epetra_Vector::ReplaceGlobalValues(int NumEntries, double * Values, int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, 0, Values, Indices, true, false));
  return(0);
}
//=========================================================================
int Epetra_Vector::ReplaceMyValues(int NumEntries, double * Values, int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, 0, Values, Indices, false, false));
  return(0);
}
//=========================================================================
int Epetra_Vector::SumIntoGlobalValues(int NumEntries, double * Values, int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, 0, Values, Indices, true, true));
  return(0);
}
//=========================================================================
int Epetra_Vector::SumIntoMyValues(int NumEntries, double * Values, int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, 0, Values, Indices, false, true));
  return(0);
}
//=========================================================================
int Epetra_Vector::ReplaceGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, BlockOffset, Values, Indices, true, false));
  return(0);
}
//=========================================================================
int Epetra_Vector::ReplaceMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, BlockOffset, Values, Indices, false, false));
  return(0);
}
//=========================================================================
int Epetra_Vector::SumIntoGlobalValues(int NumEntries, int BlockOffset, double * Values, int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, BlockOffset, Values, Indices, true, true));
  return(0);
}
//=========================================================================
int Epetra_Vector::SumIntoMyValues(int NumEntries, int BlockOffset, double * Values, int * Indices) {
  // Use the more general method below
  EPETRA_CHK_ERR(ChangeValues(NumEntries, BlockOffset, Values, Indices, false, true));
  return(0);
}
//=========================================================================
int Epetra_Vector::ChangeValues(int NumEntries, int BlockOffset, double * Values, int * Indices,
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
	Values_[entry+BlockOffset] += Values[i];
      else
	Values_[entry+BlockOffset] = Values[i];
    }
    else ierr = 1;
  }

  EPETRA_CHK_ERR(ierr);
  return(0);
}
