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

#ifndef KOKKOS_DEFAULTBLOCKSPARSEMULTIPLY_HPP
#define KOKKOS_DEFAULTBLOCKSPARSEMULTIPLY_HPP

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_VbrMatrix.hpp" 
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultBlockSparseMultiplyKernelOps.hpp"

namespace Kokkos {

  // default implementation
  template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  /*! \class DefaultBlockSparseMultiply
      \brief DefaultBlockSparseMultiply
  */
  class DefaultBlockSparseMultiply {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! DefaultBlockSparseMultiply constuctor with variable number of indices per row.
    DefaultBlockSparseMultiply(const Teuchos::RCP<Node> &node = DefaultNode::getDefaultNode());

    //! DefaultBlockSparseMultiply Destructor
    ~DefaultBlockSparseMultiply();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    Teuchos::RCP<Node> getNode() const;

    //@}

    //! @name Initialization of structure

    //@{

    //! Initialize values of matrix, using Kokkos::VbrMatrix
    void initializeValues(const VbrMatrix<Scalar,Ordinal,Node> &matrix);

    //! Clear all matrix structure and values.
    void clear();

    //@}

    //! @name Computational methods

    //@{

    //! Applies the matrix to a MultiVector, overwriting Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, MultiVector<RangeScalar,Node> &Y) const;

    //! Applies the matrix to a MultiVector, accumulating into Y.
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultBlockSparseMultiply(const DefaultBlockSparseMultiply& source);

    Teuchos::RCP<Node> node_;

    Teuchos::ArrayRCP<const Ordinal> pbuf_rptr_;
    Teuchos::ArrayRCP<const Ordinal> pbuf_cptr_;
    Teuchos::ArrayRCP<const size_t> pbuf_bptr_;
    Teuchos::ArrayRCP<const Ordinal> pbuf_bindx_;
    Teuchos::ArrayRCP<const Ordinal> pbuf_indx_;
    Teuchos::ArrayRCP<const Scalar>  pbuf_vals1D_;

    size_t numBlockRows_;
    bool valsInit_, isPacked_, isEmpty_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultBlockSparseMultiply<Scalar,Ordinal,Node>::DefaultBlockSparseMultiply(const Teuchos::RCP<Node> &node)
  : node_(node)
  , valsInit_(false)
  , isPacked_(false)
  , isEmpty_(false) {
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultBlockSparseMultiply<Scalar,Ordinal,Node>::~DefaultBlockSparseMultiply() {
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultBlockSparseMultiply<Scalar,Ordinal,Node>::initializeValues(const VbrMatrix<Scalar,Ordinal,Node> &matrix) {
    using Teuchos::ArrayRCP;
    isEmpty_ = false;
    pbuf_vals1D_ = matrix.get_values();
    pbuf_rptr_ = matrix.get_rptr();
    pbuf_cptr_ = matrix.get_cptr();
    pbuf_bptr_ = matrix.get_bptr();
    pbuf_bindx_ = matrix.get_bindx();
    pbuf_indx_ = matrix.get_indx();
    numBlockRows_ = matrix.getNumBlockRows();
    valsInit_ = true;
    isPacked_ = true;
  }


  template <class Scalar, class Ordinal, class Node>
  Teuchos::RCP<Node> DefaultBlockSparseMultiply<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }


  template <class Scalar, class Ordinal, class Node>
  void DefaultBlockSparseMultiply<Scalar,Ordinal,Node>::clear() {
    pbuf_vals1D_  = Teuchos::null;
    pbuf_rptr_    = Teuchos::null;
    pbuf_cptr_    = Teuchos::null;
    pbuf_bptr_    = Teuchos::null;
    pbuf_bindx_   = Teuchos::null;
    pbuf_indx_    = Teuchos::null;
    valsInit_ = false;
    isPacked_ = false;
    isEmpty_  = false;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultBlockSparseMultiply<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                MultiVector<RangeScalar,Node> &Y) const {
    typedef DefaultBlockSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1;
    typedef DefaultBlockSparseMultiplyOp1Transpose<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1T;
    TEST_FOR_EXCEPTION(valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else if (isPacked_ == true) {
      if (trans == Teuchos::NO_TRANS) {
        Op1 wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numBlockRows = numBlockRows_;
        wdp.vals = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.rptr = rbh.template addConstBuffer<Ordinal>(pbuf_rptr_);
        wdp.cptr = rbh.template addConstBuffer<Ordinal>(pbuf_cptr_);
        wdp.bptr = rbh.template addConstBuffer<size_t>(pbuf_bptr_);
        wdp.bindx= rbh.template addConstBuffer<Ordinal>(pbuf_bindx_);
        wdp.indx = rbh.template addConstBuffer<Ordinal>(pbuf_indx_);
        wdp.y    = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x    = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1>(0,numBlockRows_*numRHS,wdp);
      }
      else {
        //start by initializing Y = beta*Y
        DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
        Op1T wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.numBlockRows = numBlockRows_;
        wdp.vals = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.rptr = rbh.template addConstBuffer<Ordinal>(pbuf_rptr_);
        wdp.cptr = rbh.template addConstBuffer<Ordinal>(pbuf_cptr_);
        wdp.bptr = rbh.template addConstBuffer<size_t>(pbuf_bptr_);
        wdp.bindx= rbh.template addConstBuffer<Ordinal>(pbuf_bindx_);
        wdp.indx = rbh.template addConstBuffer<Ordinal>(pbuf_indx_);
        wdp.y    = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x    = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1T>(0,numRHS,wdp);
      }
    }
    else {
      throw std::runtime_error("DefaultBlockSparseMultiply ERROR, not implemented for non-packed '2D' storage.");
    }
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultBlockSparseMultiply<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const {
    typedef DefaultBlockSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1;
    typedef DefaultBlockSparseMultiplyOp1Transpose<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1T;
    TEST_FOR_EXCEPTION(valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else if (isPacked_ == true) {
      if (trans == Teuchos::NO_TRANS) {
        Op1 wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numBlockRows = numBlockRows_;
        wdp.vals = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.rptr = rbh.template addConstBuffer<Ordinal>(pbuf_rptr_);
        wdp.cptr = rbh.template addConstBuffer<Ordinal>(pbuf_cptr_);
        wdp.bptr = rbh.template addConstBuffer<size_t>(pbuf_bptr_);
        wdp.bindx= rbh.template addConstBuffer<Ordinal>(pbuf_bindx_);
        wdp.indx = rbh.template addConstBuffer<Ordinal>(pbuf_indx_);
        wdp.y    = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x    = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1>(0,numBlockRows_*numRHS,wdp);
      }
      else {
        //start by initializing Y = beta*Y
        if (beta == Teuchos::ScalarTraits<RangeScalar>::zero()) {
          DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
        }
        else {
          DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta,Y);
        }
        Op1T wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.numBlockRows = numBlockRows_;
        wdp.vals = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.rptr = rbh.template addConstBuffer<Ordinal>(pbuf_rptr_);
        wdp.cptr = rbh.template addConstBuffer<Ordinal>(pbuf_cptr_);
        wdp.bptr = rbh.template addConstBuffer<size_t>(pbuf_bptr_);
        wdp.bindx= rbh.template addConstBuffer<Ordinal>(pbuf_bindx_);
        wdp.indx = rbh.template addConstBuffer<Ordinal>(pbuf_indx_);
        wdp.y    = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x    = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1T>(0,numRHS,wdp);
      }
    }
    else {
      throw std::runtime_error("DefaultBlockSparseMultiply ERROR, not implemented for non-packed '2D' storage.");
    }
  }

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEMULTIPLY_HPP */
