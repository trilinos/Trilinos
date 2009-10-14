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
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#ifndef KERNEL_PREFIX
  #define KERNEL_PREFIX
#endif

namespace Kokkos {

  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseMultiplyOp1 {
    // mat data
    const size_t  *offsets;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    Scalar        alpha, beta;
    size_t numRows;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      const size_t rhs = (i - row) / numRows;
      Scalar tmp = 0;
      const DomainScalar *xj = x + rhs * xstride;
      RangeScalar        *yj = y + rhs * ystride;
      for (size_t c=offsets[row]; c != offsets[row+1]; ++c) {
        tmp += vals[c] * xj[inds[c]];
      }
      Scalar tmp2 = beta * yj[row];
      yj[row] = (RangeScalar)(alpha * tmp + tmp2);
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseTransposeMultiplyOp1 {
    // mat data
    const size_t  *offsets;
    const Ordinal *inds;
    const Scalar  *vals;
    // matvec params
    Scalar        alpha, beta;
    size_t numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // multiply entire matrix for rhs i
      const size_t rhs = i;
      const DomainScalar *xj = x + rhs * xstride;
      RangeScalar        *yj = y + rhs * ystride;
      for (size_t row=0; row < numCols; ++row) {
        yj[row] = (RangeScalar)(yj[row] * beta);
      }
      for (size_t row=0; row < numRows; ++row) {
        for (size_t c=offsets[row]; c != offsets[row+1]; ++c) {
          yj[inds[c]] += (RangeScalar)(alpha * vals[c] * xj[row]);
        }
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseMultiplyOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    // matvec params
    Scalar        alpha, beta;
    size_t numRows;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      const size_t row = i % numRows;
      const size_t rhs = (i - row) / numRows;
      Scalar tmp = 0;
      const DomainScalar *xj = x + rhs * xstride;
      RangeScalar        *yj = y + rhs * ystride;
      const Scalar  *curval = vals_beg[row];
      const Ordinal *curind = inds_beg[row];
      for (size_t j=0; j != numEntries[row]; ++j) {
        tmp += (curval[j]) * xj[curind[j]];
      }
      Scalar tmp2 = beta * yj[row];
      yj[row] = (Scalar)(alpha * tmp + tmp2);
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseTransposeMultiplyOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    // matvec params
    Scalar        alpha, beta;
    size_t numRows, numCols;
    // mv data
    const DomainScalar  *x;
    RangeScalar         *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // multiply entire matrix for rhs i
      const size_t rhs = i;
      const DomainScalar *xj = x + rhs * xstride;
      RangeScalar        *yj = y + rhs * ystride;
      for (size_t row=0; row < numCols; ++row) {
        yj[row] = (RangeScalar)(yj[row] * beta);
      }
      for (size_t row=0; row < numRows; ++row) {
        const Scalar  *rowval = vals_beg[row];
        const Ordinal *rowind = inds_beg[row];
        for (size_t j=0; j != numEntries[row]; ++j) {
          yj[rowind[j]] += (RangeScalar)(alpha * rowval[j] * xj[row]);
        }
      }
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

    // we do this one of two ways: 
    // 1D/packed: array of offsets, pointer for ordinals, pointer for values. obviously the smallest footprint.
    Teuchos::ArrayRCP<const Ordinal> pbuf_inds1D_;
    Teuchos::ArrayRCP<const size_t>  pbuf_offsets1D_;
    Teuchos::ArrayRCP<const Scalar>  pbuf_vals1D_;
    // 2D: array of pointers
    Teuchos::ArrayRCP<const Ordinal *> pbuf_inds2D_;
    Teuchos::ArrayRCP<const Scalar  *> pbuf_vals2D_;
    Teuchos::ArrayRCP<size_t>          pbuf_numEntries_;
    
    size_t numRows_;
    bool indsInit_, valsInit_, isPacked_, isEmpty_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseMultiply<Scalar,Ordinal,Node>::DefaultSparseMultiply(const Teuchos::RCP<Node> &node)
  : node_(node)
  , indsInit_(false)
  , valsInit_(false)
  , isPacked_(false) 
  , isEmpty_(false) {
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
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(cv != Teuchos::View, std::runtime_error,
        Teuchos::typeName(*this) << "::initializeStructure(): requires View access.");
    TEST_FOR_EXCEPTION(indsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (graph.isPacked()) {
      isEmpty_ = false;
      isPacked_ = true;
      pbuf_inds1D_    = graph.get1DIndices();
      pbuf_offsets1D_ = graph.get1DOffsets();
    }
    else {
      isEmpty_ = false;
      isPacked_ = false;
      pbuf_inds2D_     = node_->template allocBuffer<const Ordinal *>(numRows_);
      pbuf_numEntries_ = node_->template allocBuffer<size_t>(numRows_);
      ArrayRCP<const Ordinal *> inds2Dview = node_->template viewBufferNonConst<const Ordinal *>(WriteOnly, numRows_, pbuf_inds2D_);
      ArrayRCP<         size_t> numEntview = node_->template viewBufferNonConst<         size_t>(WriteOnly, numRows_, pbuf_numEntries_);
      for (size_t r=0; r < numRows_; ++r) {
        ArrayRCP<const Ordinal> rowinds = graph.get2DIndices(r);
        if (rowinds != Teuchos::null) {
          inds2Dview[r] = rowinds.getRawPtr();
          numEntview[r] = rowinds.size();
        }
        else {
          inds2Dview[r] = NULL;
          numEntview[r] = 0;
        }
      }
    }
    indsInit_ = true;
    return Teuchos::View;
  }


  template <class Scalar, class Ordinal, class Node>
  Teuchos::DataAccess DefaultSparseMultiply<Scalar,Ordinal,Node>::initializeValues(CrsMatrix<Scalar,Node> &matrix, Teuchos::DataAccess cv) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(cv != Teuchos::View, std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): requires View access.");
    TEST_FOR_EXCEPTION(valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): values already initialized.");
    TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || (!isEmpty_ && isPacked_ != matrix.isPacked()), std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (isEmpty_ || matrix.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (matrix.isPacked()) {
      isEmpty_ = false;
      pbuf_vals1D_ = matrix.get1DValues();
    }
    else {
      isEmpty_ = false;
      pbuf_vals2D_ = node_->template allocBuffer<const Scalar *>(numRows_);
      ArrayRCP<const Scalar *> vals2Dview = node_->template viewBufferNonConst<const Scalar *>(WriteOnly, numRows_, pbuf_vals2D_);
      for (size_t r=0; r < numRows_; ++r) {
        ArrayRCP<const Scalar> rowvals = matrix.get2DValues(r);
        if (rowvals != Teuchos::null) {
          vals2Dview[r] = rowvals.getRawPtr();
        }
        else {
          vals2Dview[r] = NULL;
        }
      }
    }
    valsInit_ = true;
    return Teuchos::View;
  }


  template <class Scalar, class Ordinal, class Node>
  Teuchos::RCP<Node> DefaultSparseMultiply<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }


  template <class Scalar, class Ordinal, class Node>
  void DefaultSparseMultiply<Scalar,Ordinal,Node>::clear() {
    pbuf_inds1D_     = Teuchos::null;
    pbuf_offsets1D_  = Teuchos::null;
    pbuf_vals1D_     = Teuchos::null;
    pbuf_inds2D_     = Teuchos::null;
    pbuf_vals2D_     = Teuchos::null;
    pbuf_numEntries_ = Teuchos::null;
    indsInit_ = false;
    valsInit_ = false;
    isPacked_ = false;
    isEmpty_  = false;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultSparseMultiply<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                Scalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                Scalar beta, MultiVector<RangeScalar,Node> &Y) const {
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1D;
    typedef DefaultSparseMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar>  Op2D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar> TOp1D;
    typedef DefaultSparseTransposeMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar> TOp2D;
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= alpha * 0 * X + beta * Y 
      //   <= beta * Y
      // TODO: this neglects NaNs in X, which don't satisfy 0*NaN == 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else if (isPacked_ == true) {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.offsets = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.offsets = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op2D wdp;
        rbh.begin();
        wdp.numRows = numRows_;
        wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(pbuf_inds2D_);
        wdp.vals_beg   = rbh.template addConstBuffer<const Scalar *>(pbuf_vals2D_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        rbh.end();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op2D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp2D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(pbuf_inds2D_);
        wdp.vals_beg   = rbh.template addConstBuffer<const Scalar *>(pbuf_vals2D_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp2D>(0,numRHS,wdp);
      }
    }
    return;
  }

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEMULTIPLY_H */
