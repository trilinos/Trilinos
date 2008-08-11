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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_MultiVector.h"
#include "Ifpack_Graph.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_AMDReordering.h"

extern "C" {
#include <amesos_amd.h>
}

//==============================================================================
Ifpack_AMDReordering::
Ifpack_AMDReordering() :
  NumMyRows_(0),
  IsComputed_(false)
{
}

//==============================================================================
Ifpack_AMDReordering::
Ifpack_AMDReordering(const Ifpack_AMDReordering& RHS) :
  NumMyRows_(RHS.NumMyRows()),
  IsComputed_(RHS.IsComputed())
{
  Reorder_.resize(NumMyRows());
  InvReorder_.resize(NumMyRows());
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    Reorder_[i] = RHS.Reorder(i);
    InvReorder_[i] = RHS.InvReorder(i);
  }
}

//==============================================================================
Ifpack_AMDReordering& Ifpack_AMDReordering::
operator=(const Ifpack_AMDReordering& RHS)
{
  if (this == &RHS) {
    return (*this);
  }

  NumMyRows_ = RHS.NumMyRows(); // set number of local rows
  IsComputed_ = RHS.IsComputed();
  // resize vectors, and copy values from RHS
  Reorder_.resize(NumMyRows()); 
  InvReorder_.resize(NumMyRows());
  if (IsComputed()) {
    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      Reorder_[i] = RHS.Reorder(i);
      InvReorder_[i] = RHS.InvReorder(i);
    }
  }
  return (*this);
}

//==============================================================================
int Ifpack_AMDReordering::
SetParameter(const string Name, const int Value)
{
  return(0);
}

//==============================================================================
int Ifpack_AMDReordering::
SetParameter(const string Name, const double Value)
{
  return(0);
}

//==============================================================================
int Ifpack_AMDReordering::
SetParameters(Teuchos::ParameterList& List)
{
  return(0);
}

//==============================================================================
int Ifpack_AMDReordering::Compute(const Epetra_RowMatrix& Matrix)
{
  Ifpack_Graph_Epetra_RowMatrix Graph(Teuchos::rcp(&Matrix,false));

  IFPACK_CHK_ERR(Compute(Graph));

  return(0);
}

//==============================================================================
int Ifpack_AMDReordering::Compute(const Ifpack_Graph& Graph)
{
  IsComputed_ = false;
  NumMyRows_ = Graph.NumMyRows();
  int NumNz = Graph.NumMyNonzeros();
  
  if (NumMyRows_ == 0)
    IFPACK_CHK_ERR(-1); // strange graph this one
  
  // Extract CRS format
  vector<int> ia(NumMyRows_+1,0);
  vector<int> ja(NumNz);
  int cnt;
  for( int i = 0; i < NumMyRows_; ++i )
  {
    int * tmpP = &ja[ia[i]];
    Graph.ExtractMyRowCopy( i, NumNz-ia[i], cnt, tmpP );
    ia[i+1] = ia[i] + cnt;
  }

  // Trim down to local only
  vector<int> iat(NumMyRows_+1);
  vector<int> jat(NumNz);
  int loc = 0;
  for( int i = 0; i < NumMyRows_; ++i )
  {
    iat[i] = loc;
    for( int j = ia[i]; j < ia[i+1]; ++j )
    {
      if( ja[j] < NumMyRows_ )
        jat[loc++] = ja[j];
      else
	break;
    }
  }
  iat[NumMyRows_] = loc;

  // Compute AMD permutation
  Reorder_.resize(NumMyRows_);
  vector<double> info(AMD_INFO);

  amesos_amd_order( NumMyRows_, &iat[0], &jat[0], &Reorder_[0], NULL, &info[0] );

  if( info[AMD_STATUS] == AMD_INVALID )
    cout << "AMD ORDERING: Invalid!!!!\n";

  // Build inverse reorder (will be used by ExtractMyRowCopy() 
  InvReorder_.resize(NumMyRows_);

  for (int i = 0 ; i < NumMyRows_ ; ++i)
    InvReorder_[i] = -1;

  for (int i = 0 ; i < NumMyRows_ ; ++i)
    InvReorder_[Reorder_[i]] = i;

  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    if (InvReorder_[i] == -1)
      IFPACK_CHK_ERR(-1);
  }

  IsComputed_ = true;
  return(0);
}

//==============================================================================
int Ifpack_AMDReordering::Reorder(const int i) const
{
#ifdef IFPACK_ABC
  if (!IsComputed())
    IFPACK_CHK_ERR(-1);
  if ((i < 0) || (i >= NumMyRows_))
    IFPACK_CHK_ERR(-1);
#endif

  return(Reorder_[i]);
}

//==============================================================================
int Ifpack_AMDReordering::InvReorder(const int i) const
{
#ifdef IFPACK_ABC
  if (!IsComputed())
    IFPACK_CHK_ERR(-1);
  if ((i < 0) || (i >= NumMyRows_))
    IFPACK_CHK_ERR(-1);
#endif

  return(InvReorder_[i]);
}
//==============================================================================
int Ifpack_AMDReordering::P(const Epetra_MultiVector& Xorig,
			    Epetra_MultiVector& X) const
{  
  int NumVectors = X.NumVectors();

  for (int j = 0 ; j < NumVectors ; ++j) {
    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      int np = Reorder_[i];
      X[j][np] = Xorig[j][i];
    }
  }

  return(0);
}

//==============================================================================
int Ifpack_AMDReordering::Pinv(const Epetra_MultiVector& Xorig,
			       Epetra_MultiVector& X) const
{
  int NumVectors = X.NumVectors();

  for (int j = 0 ; j < NumVectors ; ++j) {
    for (int i = 0 ; i < NumMyRows_ ; ++i) {
      int np = Reorder_[i];
      X[j][i] = Xorig[j][np];
    }
  }

  return(0);
}

//==============================================================================
ostream& Ifpack_AMDReordering::Print(std::ostream& os) const
{
  os << "*** Ifpack_AMDReordering" << endl << endl;
  if (!IsComputed())
    os << "*** Reordering not yet computed." << endl;
  
  os << "*** Number of local rows = " << NumMyRows_ << endl;
  os << endl;
  os << "Local Row\tReorder[i]\tInvReorder[i]" << endl;
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    os << '\t' << i << "\t\t" << Reorder_[i] << "\t\t" << InvReorder_[i] << endl;
  }
   
  return(os);
}
