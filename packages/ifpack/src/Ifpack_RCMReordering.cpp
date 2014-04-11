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
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_MultiVector.h"
#include "Ifpack_Graph.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_RCMReordering.h"

//==============================================================================
Ifpack_RCMReordering::
Ifpack_RCMReordering() :
  RootNode_(0),
  NumMyRows_(0),
  IsComputed_(false)
{
}

//==============================================================================
Ifpack_RCMReordering::
Ifpack_RCMReordering(const Ifpack_RCMReordering& RHS) :
  RootNode_(RHS.RootNode()),
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
Ifpack_RCMReordering& Ifpack_RCMReordering::
operator=(const Ifpack_RCMReordering& RHS)
{
  if (this == &RHS) {
    return (*this);
  }

  NumMyRows_ = RHS.NumMyRows(); // set number of local rows
  RootNode_ = RHS.RootNode(); // set root node
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
int Ifpack_RCMReordering::
SetParameter(const string Name, const int Value)
{
  if (Name == "reorder: root node")
    RootNode_ = Value;
  return(0);
}

//==============================================================================
int Ifpack_RCMReordering::
SetParameter(const string Name, const double Value)
{
  return(0);
}

//==============================================================================
int Ifpack_RCMReordering::
SetParameters(Teuchos::ParameterList& List)
{
  RootNode_ = List.get("reorder: root node", RootNode_);
  return(0);
}

//==============================================================================
int Ifpack_RCMReordering::Compute(const Epetra_RowMatrix& Matrix)
{
  Ifpack_Graph_Epetra_RowMatrix Graph(Teuchos::rcp(&Matrix,false));

  IFPACK_CHK_ERR(Compute(Graph));

  return(0);
}

//==============================================================================
int Ifpack_RCMReordering::Compute(const Ifpack_Graph& Graph)
{
  IsComputed_ = false;
  NumMyRows_ = Graph.NumMyRows();
  
  if ((RootNode_ < 0) || (RootNode_ >= NumMyRows_))
    RootNode_ = 0;
    
  Reorder_.resize(NumMyRows_);

  // the case where one processor holds no chunk of the graph happens...
  if (!NumMyRows_) {
    InvReorder_.resize(NumMyRows_);
    IsComputed_ = true;
    return(0);
  }

  for (int i = 0 ; i < NumMyRows_ ; ++i)
    Reorder_[i] = -1;

  std::vector<int> tmp;
  tmp.push_back(RootNode_);

  int count = NumMyRows_ - 1;
  int Length = Graph.MaxMyNumEntries();
  std::vector<int> Indices(Length);
  
  Reorder_[RootNode_] = count;
  count--;

  // stop when no nodes have been added in the previous level

  while (tmp.size()) {

    std::vector<int> tmp2;

    // for each node in the previous level, look for non-marked
    // neighbors. 
    for (int i = 0 ; i < (int)tmp.size() ; ++i) {
      int NumEntries;
      IFPACK_CHK_ERR(Graph.ExtractMyRowCopy(tmp[i], Length,
					     NumEntries, &Indices[0]));

      if (Length > 1)
	std::sort(Indices.begin(), Indices.begin() + Length);

      for (int j = 0 ; j < NumEntries ; ++j) {
	int col = Indices[j];
	if (col >= NumMyRows_) 
	  continue;

	if (Reorder_[col] == -1) {
	  Reorder_[col] = count;
	  count--;
	  if (col != tmp[i]) {
	    tmp2.push_back(col);
	  }
	}
      }
    }

    // if no nodes have been found but we still have
    // rows to walk through, to localize the next -1 
    // and restart.
    // FIXME: I can replace with STL
    if ((tmp2.size() == 0) && (count != -1)) {
      for (int i = 0 ; i < NumMyRows_ ; ++i)
	if (Reorder_[i] == -1) {
	  tmp2.push_back(i);
	  Reorder_[i] = count--;
	  break;
	}
    }

    // prepare for the next level
    tmp = tmp2;
  }

  // check nothing went wrong
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    if (Reorder_[i] == -1)
      IFPACK_CHK_ERR(-1);
  }
  
  // build inverse reorder (will be used by ExtractMyRowCopy() 
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
int Ifpack_RCMReordering::Reorder(const int i) const
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
int Ifpack_RCMReordering::InvReorder(const int i) const
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
int Ifpack_RCMReordering::P(const Epetra_MultiVector& Xorig,
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
int Ifpack_RCMReordering::Pinv(const Epetra_MultiVector& Xorig,
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
ostream& Ifpack_RCMReordering::Print(std::ostream& os) const
{
  os << "*** Ifpack_RCMReordering" << endl << endl;
  if (!IsComputed())
    os << "*** Reordering not yet computed." << endl;
  
  os << "*** Number of local rows = " << NumMyRows_ << endl;
  os << "*** Root node = " << RootNode_ << endl;
  os << endl;
  os << "Local Row\tReorder[i]\tInvReorder[i]" << endl;
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    os << '\t' << i << "\t\t" << Reorder_[i] << "\t\t" << InvReorder_[i] << endl;
  }
   
  return(os);
}
