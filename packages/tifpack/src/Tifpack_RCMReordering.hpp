/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef TIFPACK_RCMREORDERING_HPP
#define TIFPACK_RCMREORDERING_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Reordering.hpp"

namespace Teuchos {
  class ParameterList;
}
class Tifpack_Graph;
class Tpetra_MultiVector;
class Tpetra_RowMatrix;

//! Tifpack_RCMReordering: reverse Cuthill-McKee reordering.

class Tifpack_RCMReordering : public Tifpack_Reordering {

public:

  //! Constructor for Tifpack_Graph's.
  Tifpack_RCMReordering();

  //! Copy Constructor.
  Tifpack_RCMReordering(const Tifpack_RCMReordering& RHS);

  //! Assignment operator.
  Tifpack_RCMReordering& operator=(const Tifpack_RCMReordering& RHS);

  //! Destructor.
  virtual ~Tifpack_RCMReordering() {};
  
  //! Sets integer parameters `Name'.
  virtual int SetParameter(const string Name, const int Value);

  //! Sets double parameters `Name'.
  virtual int SetParameter(const string Name, const double Value);
  
  //! Sets all parameters.
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Tifpack_Graph& Graph);

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Tpetra_RowMatrix& Matrix);

  //! Returns \c true is the reordering object has been successfully initialized, false otherwise.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Returns the reordered index of row \c i.
  virtual inline int Reorder(const int i) const;

  //! Returns the inverse reordered index of row \c i.
  virtual inline int InvReorder(const int i) const;

  //! Applies reordering to multivector X, whose local length equals the number of local rows.
  virtual int P(const Tpetra_MultiVector& Xorig,
		Tpetra_MultiVector& Xreord) const;

  //! Applies inverse reordering to multivector X, whose local length equals the number of local rows.
  virtual int Pinv(const Tpetra_MultiVector& Xorig,
		   Tpetra_MultiVector& Xinvreord) const;

  
  //! Prints basic information on iostream. This function is used by operator<<.
  virtual ostream& Print(std::ostream& os) const;

  //! Returns the number of local rows.
  virtual int NumMyRows() const 
  {
    return(NumMyRows_);
  }

  //! Returns the root node.
  virtual int RootNode() const 
  {
    return(RootNode_);
  }

private:
  //! Defines the root node (defaulted to 0).
  int RootNode_;
  //! Number of local rows in the graph.
  int NumMyRows_;
  //! If \c true, the reordering has been successfully computed.
  bool IsComputed_;
  //! Contains the reordering.
  std::vector<int> Reorder_;
  //! Contains the inverse reordering.
  std::vector<int> InvReorder_;
}; 

#endif
