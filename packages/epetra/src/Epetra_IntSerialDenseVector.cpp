
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


#include "Epetra_IntSerialDenseVector.h"
//=============================================================================
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(void)
  : Epetra_IntSerialDenseMatrix()
{
}

//=============================================================================
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(Epetra_DataAccess CV, int *Values, int Length)
  : Epetra_IntSerialDenseMatrix(CV, Values, Length, Length, 1)
{}
//=============================================================================
Epetra_IntSerialDenseVector::Epetra_IntSerialDenseVector(const Epetra_IntSerialDenseVector& Source)
  : Epetra_IntSerialDenseMatrix(Source)
{}
//=============================================================================
Epetra_IntSerialDenseVector::~Epetra_IntSerialDenseVector()
{}
//=========================================================================
int& Epetra_IntSerialDenseVector::operator() (int Index)  {

  if (Index>=M_) throw ReportError("Index = " +toString(Index) + 
				   " Out of Range 0 - " + toString(M_-1),-1);
  return(A_[Index]);
}

//=========================================================================
const int& Epetra_IntSerialDenseVector::operator() (int Index) const  {

  if (Index>=M_) throw ReportError("Index = " +toString(Index) + 
				   " Out of Range 0 - " + toString(M_-1),-1);
   return(A_[Index]);
}
//=========================================================================
const int& Epetra_IntSerialDenseVector::operator [] (int Index) const  {
   return(A_[Index]);
}
//=========================================================================
int& Epetra_IntSerialDenseVector::operator [] (int Index)  {

  if (Index>=M_) throw ReportError("Index = " +toString(Index) + 
				   " Out of Range 0 - " + toString(M_-1),-1);
   return(A_[Index]);
}
