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

#ifndef TIFPACK_METISREORDERING_HPP
#define TIFPACK_METISREORDERING_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Reordering.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Tpetra_Comm;
class Tifpack_Graph;
class Tpetra_Map;
class Tpetra_BlockMap;
class Tpetra_Import;

//! Tifpack_METISReordering: A class to reorder a graph using METIS.

class Tifpack_METISReordering : public Tifpack_Reordering {

public:

  //! Constructor.
  Tifpack_METISReordering();

  //! Destructor.
  virtual ~Tifpack_METISReordering() {};

  //! Sets integer parameters `Name'.
  virtual int SetParameter(const string Name, const int Value)
  {
    if (Name == "partitioner: use symmetric graph")
      UseSymmetricGraph_ = (bool)Value;
    return(0);
  }
 
  //! Sets double parameters `Name'.
  virtual int SetParameter(const string Name, const double Value)
  {
    return(0);
  };

  //! Sets all the parameters for the partitioner (none at moment).
  virtual int SetParameters(Teuchos::ParameterList& List)
  {
    UseSymmetricGraph_ = List.get("partitioner: use symmetric graph", 
				  UseSymmetricGraph_);

    return(0);
  }

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
  virtual int Reorder(const int i) const;

  //! Returns the inverse reordered index of row \c i.
  virtual int InvReorder(const int i) const;

  //! Applies reordering to multivector Xorig, whose local length equals the number of local rows, stores result in X.
  virtual int P(const Tpetra_MultiVector& Xorig,
		Tpetra_MultiVector& X) const;

  //! Applies inverse reordering to multivector Xorig, whose local length equals the number of local rows, stores result in X.
  virtual int Pinv(const Tpetra_MultiVector& Xorig,
		   Tpetra_MultiVector& X) const;

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual ostream& Print(std::ostream& os) const;

private:
  //! If \c true, the graph has to be symmetrized before calling METIS.
  bool UseSymmetricGraph_;
  //! Number of local rows in the graph.
  int NumMyRows_;
  //! If \c true, the reordering has been successfully computed.
  bool IsComputed_;
  //! Contains the reordering.
  std::vector<int> Reorder_;
  //! Contains the inverse reordering.
  std::vector<int> InvReorder_;

}; // class Tifpack_METISReordering

#endif // TIFPACK_METISREORDERING_HPP
