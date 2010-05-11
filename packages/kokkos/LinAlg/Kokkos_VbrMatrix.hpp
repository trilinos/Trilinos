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

#ifndef KOKKOS_VBRMATRIX_H
#define KOKKOS_VBRMATRIX_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

  //! Kokkos::VbrMatrix: Kokkos variable block row matrix class.

  template <class Scalar, class Node = DefaultNode::DefaultNodeType>
  class VbrMatrix {
  public:

    typedef typename Node::size_t size_t;
    typedef Scalar  ScalarType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor
    //@{

    //! 'Default' VbrMatrix constuctor.
    VbrMatrix(size_t numBlockRows, const Teuchos::RCP<Node> &node = DefaultNode::getDefaultNode());

    //! VbrMatrix Destructor
    ~VbrMatrix();

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

    //! Submit the values for one row of block-entries.
    void setBlockRow(size_t row, const Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > &blockEntry);

    //! Retrieve the values for 1D storage.
    Teuchos::ArrayRCP<const Scalar> getPackedValues() const;

    //! Retrieve the values for one block row of 2D storage.
    /** Each row is an array of block-entries, where each block-entry is an array.
     * That is why this method returns an array-of-arrays...
     */
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > getBlockRow(size_t row) const;

    //! Indicates whether or not the matrix entries are packed.
    bool isPacked() const;
  
    //! Indicates that the matrix is initialized, but empty.
    bool isEmpty() const;

    //! Return the number of block rows in the matrix.
    size_t getNumBlockRows() const;

    //! Release data associated with this matrix.
    void clear();

    //@}

  private:
    //! Copy constructor (private and not implemented)
    VbrMatrix(const VbrMatrix& source);

    Teuchos::RCP<Node> node_;
    size_t numRows_;
    bool isInitialized_, isPacked_, isEmpty_;

    Teuchos::ArrayRCP<const Scalar>                      pbuf_values1D_;
    
    Teuchos::ArrayRCP< Teuchos::ArrayRCP< Teuchos::ArrayRCP<const Scalar> > > pbuf_values2D_;
    //We use the name 'pbuf_values2D_' even though it's a 3D structure...
    //It's a 3D structure because each row is a collection of block-entries.
    //The outer-most array is an array of rows, each row is an array of block-entries,
    //and each block-entry is an array of scalars.
    //Logically speaking, it is a 2D structure of block-entries.
    //Logically/mathematically, each block-entry is a dense matrix (a rectangular
    //array).
    //In keeping with the tradition of Aztec's DVBR and Epetra_VbrMatrix, each block-entry
    //is assumed to be stored in column-major order. I.e., the scalars for a given
    //column of the block-entry are stored consecutively (in contiguous memory).
  };


  //==============================================================================
  template <class Scalar, class Node>
  VbrMatrix<Scalar,Node>::VbrMatrix(typename Node::size_t numRows, const Teuchos::RCP<Node> &node)
  : node_(node)
  , numRows_(numRows)
  , isInitialized_(false)
  , isPacked_(false)
  , isEmpty_(true) {
  }

  //==============================================================================
  template <class Scalar, class Node>
  VbrMatrix<Scalar,Node>::~VbrMatrix() {
  }

  //==============================================================================
  template <class Scalar, class Node>
  Teuchos::RCP<Node> VbrMatrix<Scalar,Node>::getNode() const {
    return node_;
  }

  //==============================================================================
  template <class Scalar, class Node>
  void VbrMatrix<Scalar,Node>::clear() { 
    pbuf_values1D_ = Teuchos::null;
    pbuf_values2D_ = Teuchos::null;
    isInitialized_ = false;
    isEmpty_       = true;
    isPacked_      = false;
  }

  //==============================================================================
  template <class Scalar, class Node>
  void VbrMatrix<Scalar,Node>::setPackedValues(
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
  void VbrMatrix<Scalar,Node>::setBlockRow(
                              typename Node::size_t row,
                              const Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > &rowvals) {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(isPacked_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::setBlockRow(): matrix is already initialized with 1D (packed) structure. Call clear() before reinitializing.");
#endif
    if (isInitialized_ == false) {
      pbuf_values2D_ = Teuchos::arcp<Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > >(numRows_);
      isInitialized_ = true;
    }
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION((row < 1 && row != 0) || row > numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::setBlockRow(): specified row is invalid.");
#endif
    isEmpty_ = isEmpty_ && (rowvals == Teuchos::null);
    pbuf_values2D_[row] = rowvals;
  }

  //==============================================================================
  template <class Scalar, class Node>
  Teuchos::ArrayRCP<const Scalar> 
  VbrMatrix<Scalar,Node>::getPackedValues() const {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(isPacked_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getPackedValues(): matrix is uninitialized or not packed.");
#endif
    return pbuf_values1D_;
  }

  //==============================================================================
  template <class Scalar, class Node>
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> >
  VbrMatrix<Scalar,Node>::getBlockRow(typename Node::size_t row) const {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(isInitialized_ == false || isPacked_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getBlockRow(): matrix is uninitialized or initialized packed.");
    TEST_FOR_EXCEPTION((row < 1 && row != 0) || row > numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::getBlockRow(): row number is invalid.");
#endif
    return pbuf_values2D_[row];
  }

  //==============================================================================
  template <class Scalar, class Node>
  bool VbrMatrix<Scalar,Node>::isPacked() const {
    return isPacked_;
  }

  //==============================================================================
  template <class Scalar, class Node>
  bool VbrMatrix<Scalar,Node>::isEmpty() const {
    return isEmpty_;
  }

  //==============================================================================
  template <class Scalar, class Node>
  typename Node::size_t VbrMatrix<Scalar,Node>::getNumBlockRows() const {
    return numRows_;
  }

} // namespace Kokkos


#endif /* KOKKOS_VBRMATRIX_H */
