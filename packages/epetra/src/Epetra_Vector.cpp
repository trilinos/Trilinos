
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
//=============================================================================
Epetra_Vector::Epetra_Vector(const Epetra_BlockMap& Map)
  : Epetra_MultiVector(Map,1) // Vector is just special case of MultiVector
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
int Epetra_Vector::ExtractCopy(double *V)
{
  return(Epetra_MultiVector::ExtractCopy(V, 1));
}

//=============================================================================
int Epetra_Vector::ExtractView(double **V)
{
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
