//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

  //! Kokkos::CrsMatrix: Kokkos compressed index sparse matrix class.

  template <class Scalar, class Node = DefaultNode::DefaultNodeType>
  class CrsMatrix {
  public:

    typedef typename Node::size_t size_t;
    typedef Scalar  ScalarType;
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

    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Submit the values for a 1D storage.
    void setPackedValues(const Teuchos::ArrayRCP<const Scalar> &allvals);

    //! Submit the values for one row of 2D storage.
    void set2DValues(size_t row, const Teuchos::ArrayRCP<const Scalar> &rowvals);

    //! Retrieve the values for a 1D storage.
    Teuchos::ArrayRCP<const Scalar> getPackedValues() const;

    //! Retrieve the values for one row of 2D storage.
    Teuchos::ArrayRCP<const Scalar> get2DValues(size_t row) const;

    //! Indicates whether or not the matrix entries are packed.
    bool isPacked() const;
  
    //! Indicates that the matrix is initialized, but empty.
    bool isEmpty() const;

    //! Return the number of rows in the matrix.
    size_t getNumRows() const;

    //! Release data associated with this matrix.
    void clear();

    //@}

  private:
    //! Copy constructor (protected and not implemented)
    CrsMatrix(const CrsMatrix& source);

    Teuchos::RCP<Node> node_;
    size_t numRows_;
    bool isInitialized_, isPacked_, isEmpty_;

    Teuchos::ArrayRCP<const Scalar>                      pbuf_values1D_;
    Teuchos::ArrayRCP< Teuchos::ArrayRCP<const Scalar> > pbuf_values2D_;
  };


  //==============================================================================
  template <class Scalar, class Node>
  CrsMatrix<Scalar,Node>::CrsMatrix(typename Node::size_t numRows, const Teuchos::RCP<Node> &node)
  : node_(node)
  , numRows_(numRows)
  , isInitialized_(false)
  , isPacked_(false)
  , isEmpty_(true) {
  }

  //==============================================================================
  template <class Scalar, class Node>
  CrsMatrix<Scalar,Node>::~CrsMatrix() {
  }

  //==============================================================================
  template <class Scalar, class Node>
  Teuchos::RCP<Node> CrsMatrix<Scalar,Node>::getNode() const {
    return node_;
  }

  //==============================================================================
  template <class Scalar, class Node>
  void CrsMatrix<Scalar,Node>::clear() { 
    pbuf_values1D_ = Teuchos::null;
    pbuf_values2D_ = Teuchos::null;
    isInitialized_ = false;
    isEmpty_       = true;
    isPacked_      = false;
  }

  //==============================================================================
  template <class Scalar, class Node>
  void CrsMatrix<Scalar,Node>::setPackedValues(
                        const Teuchos::ArrayRCP<const Scalar> &allvals) {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(isInitialized_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::setPackedValues(): matrix is already initialized. Call clear() before reinitializing.");
#endif
    isEmpty_ = (allvals == Teuchos::null);
    pbuf_values1D_ = allvals;
    isInitialized_ = true;
    isPacked_ = (allvals != Teuchos::null);
  }

  //==============================================================================
  template <class Scalar, class Node>
  void CrsMatrix<Scalar,Node>::set2DValues(
                              typename Node::size_t row,
                              const Teuchos::ArrayRCP<const Scalar> &rowvals) {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(isPacked_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::set2DValues(): matrix is already initialized with 1D structure. Call clear() before reinitializing.");
#endif
    if (isInitialized_ == false) {
      pbuf_values2D_ = Teuchos::arcp<Teuchos::ArrayRCP<const Scalar> >(numRows_);
      isInitialized_ = true;
    }
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION((row < 1 && row != 0) || row > numRows_, std::runtime_error,
        Teuchos::typeName(*this) << ":;set2DValues(): specified row is invalid.");
#endif
    isEmpty_ = isEmpty_ && (rowvals == Teuchos::null);
    pbuf_values2D_[row] = rowvals;
  }

  //==============================================================================
  template <class Scalar, class Node>
  Teuchos::ArrayRCP<const Scalar> 
  CrsMatrix<Scalar,Node>::getPackedValues() const {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(isPacked_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getPackedValues(): matrix is uninitialized or not packed.");
#endif
    return pbuf_values1D_;
  }

  //==============================================================================
  template <class Scalar, class Node>
  Teuchos::ArrayRCP<const Scalar>
  CrsMatrix<Scalar,Node>::get2DValues(typename Node::size_t row) const {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(isInitialized_ == false || isPacked_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::get2DValues(): matrix is uninitialized or initialized packed.");
    TEST_FOR_EXCEPTION((row < 1 && row != 0) || row > numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::get2DValues(): row number is invalid.");
#endif
    return pbuf_values2D_[row];
  }

  //==============================================================================
  template <class Scalar, class Node>
  bool CrsMatrix<Scalar,Node>::isPacked() const {
    return isPacked_;
  }

  //==============================================================================
  template <class Scalar, class Node>
  bool CrsMatrix<Scalar,Node>::isEmpty() const {
    return isEmpty_;
  }

  //==============================================================================
  template <class Scalar, class Node>
  typename Node::size_t CrsMatrix<Scalar,Node>::getNumRows() const {
    return numRows_;
  }

} // namespace Kokkos


#endif /* KOKKOS_CRSMATRIX_H */
