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

#ifndef KOKKOS_VBRMATRIX_HPP
#define KOKKOS_VBRMATRIX_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

  //! Kokkos::VbrMatrix: Kokkos variable block row matrix class.

  template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  class VbrMatrix {
  public:

    typedef typename Node::size_t size_t;
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
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

    //! Submit the values for 1D storage.
    void setPackedValues(const Teuchos::ArrayRCP<const Scalar>& allvals,
                         const Teuchos::ArrayRCP<const Ordinal>& rptr,
                         const Teuchos::ArrayRCP<const Ordinal>& cptr,
                         const Teuchos::ArrayRCP<const size_t>& bptr,
                         const Teuchos::ArrayRCP<const Ordinal>& bindx,
                         const Teuchos::ArrayRCP<const Ordinal>& indx);

    //! Submit the values for one row of block-entries.
    void setBlockRow(size_t row, const Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > &blockEntry);

    //! Indicates whether or not the matrix entries are packed.
    bool isPacked() const;
  
    //! Indicates that the matrix is initialized, but empty.
    bool isEmpty() const;

    //! Return the number of block rows in the matrix.
    size_t getNumBlockRows() const;

    //! Release data associated with this matrix.
    void clear();

    const Teuchos::ArrayRCP<const Scalar>& get_values() const { return pbuf_values1D_; }
    const Teuchos::ArrayRCP<const Ordinal>& get_rptr() const { return pbuf_rptr_; }
    const Teuchos::ArrayRCP<const Ordinal>& get_cptr() const { return pbuf_cptr_; }
    const Teuchos::ArrayRCP<const size_t>& get_bptr() const { return pbuf_bptr_; }
    const Teuchos::ArrayRCP<const Ordinal>& get_bindx() const { return pbuf_bindx_; }
    const Teuchos::ArrayRCP<const Ordinal>& get_indx() const { return pbuf_indx_; }

    //@}

  private:
    //! Copy constructor (private and not implemented)
    VbrMatrix(const VbrMatrix& source);

    Teuchos::RCP<Node> node_;
    size_t numRows_;
    bool isInitialized_, isPacked_, isEmpty_;

    Teuchos::ArrayRCP<const Scalar> pbuf_values1D_;
    Teuchos::ArrayRCP<const Ordinal> pbuf_rptr_;
    Teuchos::ArrayRCP<const Ordinal> pbuf_cptr_;
    Teuchos::ArrayRCP<const size_t> pbuf_bptr_;
    Teuchos::ArrayRCP<const Ordinal> pbuf_bindx_;
    Teuchos::ArrayRCP<const Ordinal> pbuf_indx_;
    
    //Logically/mathematically, each block-entry is a dense matrix (a rectangular
    //array).
    //In keeping with the tradition of Aztec's DVBR and Epetra_VbrMatrix, each block-entry
    //is assumed to be stored in column-major order. I.e., the scalars for a given
    //column of the block-entry are stored consecutively (in contiguous memory).
  };


  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  VbrMatrix<Scalar,Ordinal,Node>::VbrMatrix(typename Node::size_t numRows, const Teuchos::RCP<Node> &node)
  : node_(node)
  , numRows_(numRows)
  , isInitialized_(false)
  , isPacked_(false)
  , isEmpty_(true) {
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  VbrMatrix<Scalar,Ordinal,Node>::~VbrMatrix() {
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  Teuchos::RCP<Node> VbrMatrix<Scalar,Ordinal,Node>::getNode() const {
    return node_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void VbrMatrix<Scalar,Ordinal,Node>::clear() { 
    pbuf_values1D_ = Teuchos::null;
    pbuf_rptr_ = Teuchos::null;
    pbuf_cptr_ = Teuchos::null;
    pbuf_bptr_ = Teuchos::null;
    pbuf_bindx_ = Teuchos::null;
    pbuf_indx_ = Teuchos::null;
    isInitialized_ = false;
    isEmpty_       = true;
    isPacked_      = false;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void VbrMatrix<Scalar,Ordinal,Node>::setPackedValues(
                        const Teuchos::ArrayRCP<const Scalar> &allvals,
                         const Teuchos::ArrayRCP<const Ordinal>& rptr,
                         const Teuchos::ArrayRCP<const Ordinal>& cptr,
                         const Teuchos::ArrayRCP<const size_t>& bptr,
                         const Teuchos::ArrayRCP<const Ordinal>& bindx,
                         const Teuchos::ArrayRCP<const Ordinal>& indx)
{
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION(isInitialized_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::setPackedValues(): matrix is already initialized. Call clear() before reinitializing.");
#endif
    isEmpty_ = (allvals == Teuchos::null);
    pbuf_values1D_ = allvals;
    pbuf_rptr_ = rptr;
    pbuf_cptr_ = cptr;
    pbuf_bptr_ = bptr;
    pbuf_bindx_ = bindx;
    pbuf_indx_ = indx;
    isInitialized_ = true;
    isPacked_ = (allvals != Teuchos::null);
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  bool VbrMatrix<Scalar,Ordinal,Node>::isPacked() const {
    return isPacked_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  bool VbrMatrix<Scalar,Ordinal,Node>::isEmpty() const {
    return isEmpty_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  typename Node::size_t VbrMatrix<Scalar,Ordinal,Node>::getNumBlockRows() const {
    return numRows_;
  }

} // namespace Kokkos


#endif /* KOKKOS_VBRMATRIX_H */
