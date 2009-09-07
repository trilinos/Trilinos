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

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_CrsGraph.hpp" 
#include "Kokkos_MultiVector.hpp"

#ifndef KERNEL_PREFIX
  #define KERNEL_PREFIX
#endif

namespace Kokkos {

  template <class Scalar, class Ordinal, class Node>
  struct DefaultSparseMultiplyOp {
    // mat data
    const size_t  *offsets;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    Scalar        alpha, beta;
    // mv data
    const Scalar  *x;
    Scalar        *y;

    inline KERNEL_PREFIX void execute(size_t i) {
      Scalar tmp = 0;
      for (size_t c=offsets[i]; c != offsets[i+1]; ++c) {
        tmp += vals[c] * x[inds[c]];
      }
      Scalar tmp2 = beta*y[i];
      y[i] = alpha*tmp + tmp2;
    }
  };


  // default implementation
  template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  class DefaultSparseMultiply {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! DefaultSparseMultiply constuctor with variable number of indices per row.
    DefaultSparseMultiply(const Teuchos::RCP<Node> &node = DefaultNode::getDefaultNode());

    //! DefaultSparseMultiply Destructor
    ~DefaultSparseMultiply();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    Teuchos::RCP<Node> getNode() const;

    //@}

    //! @name Initialization of structure

    //@{

    //! Initialize structure of matrix
    template <class GRAPH>
    Teuchos::DataAccess initializeStructure(GRAPH &graph, Teuchos::DataAccess cv);

    //! Initialize values of matrix
    template <class MATRIX>
    Teuchos::DataAccess initializeValues(MATRIX &matrix, Teuchos::DataAccess cv);

    //! Initialize structure of matrix, using Kokkos::CrsGraph
    Teuchos::DataAccess initializeStructure(CrsGraph<Ordinal,Node> &graph, Teuchos::DataAccess cv);

    //! Initialize values of matrix, using Kokkos::CrsMatrix
    Teuchos::DataAccess initializeValues(CrsMatrix<Scalar,Node> &matrix, Teuchos::DataAccess cv);

    //! Clear all matrix structure and values.
    void clear();

    //@}

    //! @name Computational methods

    //@{

    //! Applies the matrix to a MultiVector.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, Scalar alpha, const MultiVector<DomainScalar,Node> &X, Scalar beta, MultiVector<RangeScalar,Node> &Y) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultSparseMultiply(const DefaultSparseMultiply& source);

    Teuchos::RCP<Node> node_;
    mutable DefaultSparseMultiplyOp<Scalar,Ordinal,Node> op_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseMultiply<Scalar,Ordinal,Node>::DefaultSparseMultiply(const Teuchos::RCP<Node> &node)
  : node_(node) {
    op_.alpha = Teuchos::ScalarTraits<Scalar>::one();
    op_.beta =  Teuchos::ScalarTraits<Scalar>::zero();
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseMultiply<Scalar,Ordinal,Node>::~DefaultSparseMultiply() {
  }

  template<class Scalar, class Ordinal, class Node>
  template <class GRAPH>
  Teuchos::DataAccess DefaultSparseMultiply<Scalar,Ordinal,Node>::initializeStructure(GRAPH &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse graphs
    TEST_FOR_EXCEPT(true);
  }

  template<class Scalar, class Ordinal, class Node>
  template <class MATRIX>
  Teuchos::DataAccess DefaultSparseMultiply<Scalar,Ordinal,Node>::initializeValues(MATRIX &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse matrices
    TEST_FOR_EXCEPT(true);
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::DataAccess DefaultSparseMultiply<Scalar,Ordinal,Node>::initializeStructure(CrsGraph<Ordinal,Node> &graph, Teuchos::DataAccess cv) {
    // FINISH
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::DataAccess DefaultSparseMultiply<Scalar,Ordinal,Node>::initializeValues(CrsMatrix<Scalar,Node> &graph, Teuchos::DataAccess cv) {
    // FINISH
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::RCP<Node> DefaultSparseMultiply<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultSparseMultiply<Scalar,Ordinal,Node>::clear() { 
    // FINISH
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultSparseMultiply<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                Scalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                Scalar beta, MultiVector<RangeScalar,Node> &Y) const {
    // TEST_FOR_EXCEPT(true); // FINISH
  }

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEMULTIPLY_H */
