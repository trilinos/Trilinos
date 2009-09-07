//@HEADER
// ************************************************************************
// 
//                Kokkos: A Fast Kernel Package
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

#ifndef KOKKOS_CRSMATRIX_H
#define KOKKOS_CRSMATRIX_H

#include <Teuchos_RCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_CrsGraph.hpp"

namespace Kokkos {

//! Kokkos::CrsMatrix: Kokkos compressed index sparse matrix class.

template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
class CrsMatrix {
public:
  typedef Scalar  ScalarType;
  typedef Ordinal OrdinalType;
  typedef Node    NodeType;

  //! @name Constructors/Destructor
  //@{

  //! Default CrsMatrix constuctor.
  CrsMatrix(size_t numRows, const Teuchos::RCP<Node> &node = DefaultNode::getDefaultNode());

  //! CrsMatrix Destructor
  ~CrsMatrix();

  //@}

  //! @name Accessor routines.
  //@{ 
  
  //! Node accessor.
  Teuchos::RCP<Node> getNode() const;

  //! Graph accessor.
  const CrsGraph<Ordinal,Node> & getGraph() const;

  //@}

  //! @name Data entry and accessor methods.
  //@{

  //! Submit the values for a 1D storage.
  void set1DValues(const Teuchos::ArrayRCP<const Scalar> &allvals);

  //! Submit the values for one row of 2D storage.
  void set2DValues(size_t row, const Teuchos::ArrayRCP<const Scalar> &rowvals);

  //! Retrieve the values for a 1D storage.
  Teuchos::ArrayRCP<const Scalar> get1DValues() const;

  //! Retrieve the values for one row of 2D storage.
  Teuchos::ArrayRCP<const Scalar> get2DValues(size_t row) const;

  //! Release data associated with this matrix.
  void clear();

  //@}

protected:

  //! Copy constructor (protected and not implemented)
  CrsMatrix(const CrsMatrix& source);

  Teuchos::RCP<Node> node_;
  size_t numRows_;
  bool isInitialized_;

  Teuchos::ArrayRCP<Scalar>                      pbuf_values1D_;
  Teuchos::ArrayRCP< Teuchos::ArrayRCP<Scalar> > pbuf_values2D_;
};


//==============================================================================
template <class Scalar, class Ordinal, class Node>
CrsMatrix<Scalar,Ordinal,Node>::CrsMatrix(size_t numRows, const Teuchos::RCP<Node> &node)
: node_(node)
, numRows_(numRows)
, isInitialized_(false) {
}

//==============================================================================
template <class Scalar, class Ordinal, class Node>
CrsMatrix<Scalar,Ordinal,Node>::~CrsMatrix() {
}

//==============================================================================
template <class Scalar, class Ordinal, class Node>
Teuchos::RCP<Node> CrsMatrix<Scalar,Ordinal,Node>::getNode() const { 
  return node_; 
}

//==============================================================================
template <class Scalar, class Ordinal, class Node>
void CrsMatrix<Scalar,Ordinal,Node>::clear() { 
  pbuf_values1D_ = Teuchos::null;
  pbuf_values2D_ = Teuchos::null;
  isInitialized_ = false;
}

} // namespace Kokkos


#endif /* KOKKOS_CRSMATRIX_H */
