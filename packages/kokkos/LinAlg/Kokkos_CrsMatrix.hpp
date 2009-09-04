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
    CrsMatrix(const Teuchos::RCP<Node> &node = DefaultNode::getDefaultNode());

    //! CrsMatrix Destructor
    ~CrsMatrix();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    Teuchos::RCP<Node> getNode() const;

    //@}

    //! @name Matrix entry methods
    //@{

    //@}

  protected:

    //! Copy constructor (protected and not implemented)
    CrsMatrix(const CrsMatrix& source);

    Teuchos::RCP<Node> node_;
  };


  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  CrsMatrix<Scalar,Ordinal,Node>::CrsMatrix(const Teuchos::RCP<Node> &node)
  : node_(node) {
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

} // namespace Kokkos


#endif /* KOKKOS_CRSMATRIX_H */
