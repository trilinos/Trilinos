/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
  
  std::vector<int> Indices(MaxNumEntries());
  std::vector<double> Values(MaxNumEntries());
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
