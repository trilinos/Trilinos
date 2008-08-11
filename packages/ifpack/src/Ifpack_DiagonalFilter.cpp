/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_DiagonalFilter.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//==============================================================================
Ifpack_DiagonalFilter::Ifpack_DiagonalFilter(const Teuchos::RefCountPtr<Epetra_RowMatrix>& Matrix,
					     double AbsoluteThreshold,
					     double RelativeThreshold) :
  A_(Matrix),
  AbsoluteThreshold_(AbsoluteThreshold),
  RelativeThreshold_(RelativeThreshold)
{
  Epetra_Time Time(Comm());
  
  pos_.resize(NumMyRows());
  val_.resize(NumMyRows());
  
  vector<int> Indices(MaxNumEntries());
  vector<double> Values(MaxNumEntries());
  int NumEntries;
  
  for (int MyRow = 0 ; MyRow < NumMyRows() ; ++MyRow) {
    
    pos_[MyRow] = -1;
    val_[MyRow] = 0.0;
    int ierr = A_->ExtractMyRowCopy(MyRow, MaxNumEntries(), NumEntries,
				    &Values[0], &Indices[0]);
    assert (ierr == 0);
    
    for (int i = 0 ; i < NumEntries ; ++i) {
      if (Indices[i] == MyRow) {
	pos_[MyRow] = i;
	val_[MyRow] = Values[i] * (RelativeThreshold_ - 1) +
	  AbsoluteThreshold_ * EPETRA_SGN(Values[i]);
      }
      break;
    }
  }
  cout << "TIME = " << Time.ElapsedTime() << endl;
}

//==============================================================================
int Ifpack_DiagonalFilter::
ExtractMyRowCopy(int MyRow, int Length, int& NumEntries, 
		 double* Values, int* Indices) const
{

  IFPACK_CHK_ERR(A_->ExtractMyRowCopy(MyRow, Length, NumEntries,
				     Values,Indices));

  if (pos_[MyRow] != -1)
    Values[pos_[MyRow]] += val_[MyRow];

  return(0);
}

//==============================================================================
int Ifpack_DiagonalFilter::
Multiply(bool TransA, const Epetra_MultiVector& X, 
	 Epetra_MultiVector& Y) const
{

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  IFPACK_CHK_ERR(A_->Multiply(TransA, X, Y));

  for (int v = 0 ; v < X.NumVectors() ; ++v)
    for (int i = 0 ; i < NumMyRows() ; ++i)
      Y[v][i] += val_[i] * X[v][i];
      

  return(0);
}
