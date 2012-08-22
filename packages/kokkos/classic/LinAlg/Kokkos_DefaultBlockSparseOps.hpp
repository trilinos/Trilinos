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

#ifndef KOKKOS_DEFAULTBLOCKSPARSEOPS_HPP
#define KOKKOS_DEFAULTBLOCKSPARSEOPS_HPP

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_VbrMatrix.hpp" 
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultBlockSparseKernelOps.hpp"

namespace Kokkos {

  // default implementation
  template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  /*! \class DefaultBlockSparseOps
      \brief DefaultBlockSparseOps
  */
  class DefaultBlockSparseOps {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! DefaultBlockSparseOps constuctor with variable number of indices per row.
    DefaultBlockSparseOps(const RCP<Node> &node = DefaultNode::getDefaultNode());

    //! DefaultBlockSparseOps Destructor
    ~DefaultBlockSparseOps();

    //@}

    //! @name Accessor routines.
    //@{ 
    
    //! Node accessor.
    RCP<Node> getNode() const;

    //@}

    //! @name Initialization of structure

    //@{

    //! Initialize values of matrix, using VbrMatrix
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

    //! Triangular solve: find x such that A*x=y, only if A is triangular.
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo triang, Teuchos::EDiag diag,
            const MultiVector<DomainScalar,Node>& Y, MultiVector<RangeScalar,Node>& X) const;
    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultBlockSparseOps(const DefaultBlockSparseOps& source);

    RCP<Node> node_;

    ArrayRCP<const Ordinal> pbuf_rptr_;
    ArrayRCP<const Ordinal> pbuf_cptr_;
    ArrayRCP<const size_t> pbuf_bptr_;
    ArrayRCP<const Ordinal> pbuf_bindx_;
    ArrayRCP<const Ordinal> pbuf_indx_;
    ArrayRCP<const Scalar>  pbuf_vals1D_;

    size_t numBlockRows_;
    bool valsInit_, isPacked_, isEmpty_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultBlockSparseOps<Scalar,Ordinal,Node>::DefaultBlockSparseOps(const RCP<Node> &node)
  : node_(node)
  , valsInit_(false)
  , isPacked_(false)
  , isEmpty_(false) {
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultBlockSparseOps<Scalar,Ordinal,Node>::~DefaultBlockSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultBlockSparseOps<Scalar,Ordinal,Node>::initializeValues(const VbrMatrix<Scalar,Ordinal,Node> &matrix) {
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
  RCP<Node> DefaultBlockSparseOps<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }


  template <class Scalar, class Ordinal, class Node>
  void DefaultBlockSparseOps<Scalar,Ordinal,Node>::clear() {
    pbuf_vals1D_  = null;
    pbuf_rptr_    = null;
    pbuf_cptr_    = null;
    pbuf_bptr_    = null;
    pbuf_bindx_   = null;
    pbuf_indx_    = null;
    valsInit_ = false;
    isPacked_ = false;
    isEmpty_  = false;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultBlockSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                MultiVector<RangeScalar,Node> &Y) const {
    typedef DefaultBlockSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1;
    typedef DefaultBlockSparseMultiplyOp1Transpose<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1T;
    TEUCHOS_TEST_FOR_EXCEPTION(valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
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
      throw std::runtime_error("DefaultBlockSparseOps ERROR, not implemented for non-packed '2D' storage.");
    }
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultBlockSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const {
    typedef DefaultBlockSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1;
    typedef DefaultBlockSparseMultiplyOp1Transpose<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1T;
    TEUCHOS_TEST_FOR_EXCEPTION(valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
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
      throw std::runtime_error("DefaultBlockSparseOps::Multiply ERROR, not implemented for non-packed '2D' storage.");
    }
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultBlockSparseOps<Scalar,Ordinal,Node>::solve(
                                Teuchos::ETransp trans, 
                                Teuchos::EUplo triang, 
                                Teuchos::EDiag diag, 
                                const MultiVector<DomainScalar,Node> &Y, 
                                MultiVector<RangeScalar,Node> &X) const {
    typedef DefaultBlockSparseSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op;
    typedef DefaultBlockSparseTransposeSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  OpT;
    TEUCHOS_TEST_FOR_EXCEPTION(valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // X <= Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else if (isPacked_ == true) {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.upper   = (triang == Teuchos::UPPER_TRI);
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG);
        wdp.numBlockRows = numBlockRows_;
        wdp.vals = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.rptr = rbh.template addConstBuffer<Ordinal>(pbuf_rptr_);
        wdp.cptr = rbh.template addConstBuffer<Ordinal>(pbuf_cptr_);
        wdp.bptr = rbh.template addConstBuffer<size_t>(pbuf_bptr_);
        wdp.bindx= rbh.template addConstBuffer<Ordinal>(pbuf_bindx_);
        wdp.indx = rbh.template addConstBuffer<Ordinal>(pbuf_indx_);
        wdp.x    = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y    = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op>(0,numRHS,wdp);
      }
      else {
        OpT wdp;
        rbh.begin();
        wdp.upper   = (triang == Teuchos::UPPER_TRI);
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG);
        wdp.numBlockRows = numBlockRows_;
        wdp.vals = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.rptr = rbh.template addConstBuffer<Ordinal>(pbuf_rptr_);
        wdp.cptr = rbh.template addConstBuffer<Ordinal>(pbuf_cptr_);
        wdp.bptr = rbh.template addConstBuffer<size_t>(pbuf_bptr_);
        wdp.bindx= rbh.template addConstBuffer<Ordinal>(pbuf_bindx_);
        wdp.indx = rbh.template addConstBuffer<Ordinal>(pbuf_indx_);
        wdp.x    = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y    = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<OpT>(0,numRHS,wdp);
      }
    }
    else {
      throw std::runtime_error("DefaultBlockSparseOps::solve ERROR, not implemented for non-packed '2D' storage.");
    }
  }

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */
