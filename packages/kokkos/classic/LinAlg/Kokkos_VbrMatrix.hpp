//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER

#ifndef KOKKOS_VBRMATRIX_HPP
#define KOKKOS_VBRMATRIX_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

namespace Kokkos {

  //! Kokkos variable block row matrix class.

  template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  class VbrMatrix {
  public:

    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor
    //@{

    //! 'Default' VbrMatrix constuctor.
    VbrMatrix(size_t numBlockRows, const RCP<Node> &node = DefaultNode::getDefaultNode());

    //! VbrMatrix Destructor
    ~VbrMatrix();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    RCP<Node> getNode() const;

    //@}

    //! @name Data entry and accessor methods.
    //@{

    //! Submit the values for 1D storage.
    void setPackedValues(const ArrayRCP<const Scalar>& allvals,
                         const ArrayRCP<const Ordinal>& rptr,
                         const ArrayRCP<const Ordinal>& cptr,
                         const ArrayRCP<const size_t>& bptr,
                         const ArrayRCP<const Ordinal>& bindx,
                         const ArrayRCP<const Ordinal>& indx);

    //! Submit the values for one row of block-entries.
    void setBlockRow(size_t row, const ArrayRCP<ArrayRCP<const Scalar> > &blockEntry);

    //! Indicates whether or not the matrix entries are packed.
    bool isPacked() const;
  
    //! Indicates that the matrix is initialized, but empty.
    bool isEmpty() const;

    //! Return the number of block rows in the matrix.
    size_t getNumBlockRows() const;

    //! Release data associated with this matrix.
    void clear();

    const ArrayRCP<const Scalar>& get_values() const { return pbuf_values1D_; }
    const ArrayRCP<const Ordinal>& get_rptr() const { return pbuf_rptr_; }
    const ArrayRCP<const Ordinal>& get_cptr() const { return pbuf_cptr_; }
    const ArrayRCP<const size_t>& get_bptr() const { return pbuf_bptr_; }
    const ArrayRCP<const Ordinal>& get_bindx() const { return pbuf_bindx_; }
    const ArrayRCP<const Ordinal>& get_indx() const { return pbuf_indx_; }

    //@}

  private:
    //! Copy constructor (private and not implemented)
    VbrMatrix(const VbrMatrix& source);

    RCP<Node> node_;
    size_t numRows_;
    bool isInitialized_, isPacked_, isEmpty_;

    ArrayRCP<const Scalar> pbuf_values1D_;
    ArrayRCP<const Ordinal> pbuf_rptr_;
    ArrayRCP<const Ordinal> pbuf_cptr_;
    ArrayRCP<const size_t> pbuf_bptr_;
    ArrayRCP<const Ordinal> pbuf_bindx_;
    ArrayRCP<const Ordinal> pbuf_indx_;
    
    //Logically/mathematically, each block-entry is a dense matrix (a rectangular
    //array).
    //In keeping with the tradition of Aztec's DVBR and Epetra_VbrMatrix, each block-entry
    //is assumed to be stored in column-major order. I.e., the scalars for a given
    //column of the block-entry are stored consecutively (in contiguous memory).
  };


  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  VbrMatrix<Scalar,Ordinal,Node>::VbrMatrix(size_t numRows, const RCP<Node> &node)
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
  RCP<Node> VbrMatrix<Scalar,Ordinal,Node>::getNode() const {
    return node_;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void VbrMatrix<Scalar,Ordinal,Node>::clear() { 
    pbuf_values1D_ = null;
    pbuf_rptr_ = null;
    pbuf_cptr_ = null;
    pbuf_bptr_ = null;
    pbuf_bindx_ = null;
    pbuf_indx_ = null;
    isInitialized_ = false;
    isEmpty_       = true;
    isPacked_      = false;
  }

  //==============================================================================
  template <class Scalar, class Ordinal, class Node>
  void VbrMatrix<Scalar,Ordinal,Node>::setPackedValues(
                        const ArrayRCP<const Scalar> &allvals,
                         const ArrayRCP<const Ordinal>& rptr,
                         const ArrayRCP<const Ordinal>& cptr,
                         const ArrayRCP<const size_t>& bptr,
                         const ArrayRCP<const Ordinal>& bindx,
                         const ArrayRCP<const Ordinal>& indx)
{
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(isInitialized_ == true, std::runtime_error,
        Teuchos::typeName(*this) << "::setPackedValues(): matrix is already initialized. Call clear() before reinitializing.");
#endif
    isEmpty_ = (allvals == null);
    pbuf_values1D_ = allvals;
    pbuf_rptr_ = rptr;
    pbuf_cptr_ = cptr;
    pbuf_bptr_ = bptr;
    pbuf_bindx_ = bindx;
    pbuf_indx_ = indx;
    isInitialized_ = true;
    isPacked_ = (allvals != null);
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
  size_t VbrMatrix<Scalar,Ordinal,Node>::getNumBlockRows() const {
    return numRows_;
  }

} // namespace Kokkos


#endif /* KOKKOS_VBRMATRIX_H */
