//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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

#ifndef KOKKOS_DEFAULTSPARSEMULTIPLY_H
#define KOKKOS_DEFAULTSPARSEMULTIPLY_H

#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_MultiVector.hpp"

namespace Kokkos {

  template<class Scalar, class Ordinal, class Node>
  class DefaultSparseMultiply {
  public:

    //! @name Constructors/Destructor

    //@{

    //! DefaultSparseMultiply constuctor with variable number of indices per row.
    DefaultSparseMultiply();

    //! DefaultSparseMultiply Destructor
    ~DefaultSparseMultiply();

    //@}

    //! @name Initialization of structure

    //@{

    //! Initialize structure of matrix
    template <class Matrix>
    int initializeStructure(const Matrix& A);

    //! Initialize values of matrix
    template <class Matrix>
    int initializeValues(const Matrix& A);

    //@}

    //! @name Computational methods

    //@{

    //! Applies the matrix to a MultiVector.
    int apply(const MultiVector<Scalar,Ordinal,Node> &X, MultiVector<Scalar,Ordinal,Node> &Y) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultSparseMultiply(const DefaultSparseMultiply& source);

    Ordinal numRows_;
    Ordinal numCols_;
    Ordinal numEntries_;

    typename Node::template buffer<const Ordinal>::buffer_t offsets_, indices_;
    typename Node::template buffer<const Scalar>::buffer_t values_;
  };

  template <class Scalar, class Ordinal, class Node>
  DefaultSparseMultiply::DefaultSparseMultiply()
  {
  }

  template <class Scalar, class Ordinal, class Node>
  DefaultSparseMultiply::~DefaultSparseMultiply()
  {
  }

  template <class Scalar, class Ordinal, class Node>
  template <class Matrix>
  int DefaultSparseMultiply::initializeStructure(const Matrix& A)
  {
    return 0;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class Matrix>
  int DefaultSparseMultiply::initializeValues(const Matrix& A)
  {
    return 0;
  }

  template <class Scalar, class Ordinal, class Node>
  int DefaultSparseMultiply::apply(const MultiVector<Scalar,Ordinal,Node> &X, MultiVector<Scalar,Ordinal,Node> &Y) const
  {
    return 0;
  }

} // namespace Kokkos
#endif /* KOKKOS_DEFAULTSPARSEMULTIPLY_H */
