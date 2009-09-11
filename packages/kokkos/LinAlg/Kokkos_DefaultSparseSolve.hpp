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

#ifndef KOKKOS_BASESPARSESOLVE_H
#define KOKKOS_BASESPARSESOLVE_H

#include <Teuchos_BLAS_types.hpp>
#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_CrsGraph.hpp" 
#include "Kokkos_MultiVector.hpp"

#ifndef KERNEL_PREFIX
  #define KERNEL_PREFIX
#endif

namespace Kokkos {

  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseSolveOp1 {
    // mat data
    const size_t  *offsets;
    const Ordinal *inds;
    const Scalar  *vals;
    size_t numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar  *x;
    const RangeScalar *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // solve rhs i for lhs i
      const size_t rhs = i;
      DomainScalar      *xj = x + rhs * xstride;
      const RangeScalar *yj = y + rhs * ystride;
      // unitDiag indictes whether we neglect the diagonal row entry and scale by it
      // or utilize all row entries and implicitly scale by a unit diagonal (i.e., don't scale)
      // upper (versus lower) will determine the ordering of the solve and the location of the diagonal
      // upper -> diagonal is first entry on row
      // lower -> diagonal is last entry on row
      // 
      // upper triangular requires backwards substition, solving in reverse order
      // must unroll the last iteration, because decrement results in wrap-around
      if (upper && unitDiag) {
        // upper + unit
        xj[numRows-1] = yj[numRows-1];
        for (size_t row=numRows-2; row != 0; --row) {
          xj[row] = yj[row];
          for (size_t c=offsets[row]; c != offsets[row+1]; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
        }
        xj[0] = yj[0];
        for (size_t c=offsets[0]; c != offsets[1]; ++c) {
          xj[0] -= vals[c] * xj[inds[c]];
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit
        xj[numRows-1] = yj[numRows-1] / vals[offsets[numRows-1]];
        for (size_t row=numRows-2; row != 0; --row) {
          xj[row] = yj[row];
          size_t d = offsets[row], e = offsets[row+1];
          Scalar dval = vals[d];
          for (size_t c=d+1; c != e; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
          xj[row] /= dval;
        }
        xj[0] = yj[0];
        size_t d = offsets[0], e = offsets[1];
        Scalar dval = vals[d];
        for (size_t c=d+1; c != e; ++c) {
          xj[0] -= vals[c] * xj[inds[c]];
        }
        xj[0] /= dval;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        xj[0] = yj[0];
        for (size_t row=1; row < numRows; ++row) {
          xj[row] = yj[row];
          for (size_t c=offsets[row]; c != offsets[row+1]; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit
        xj[0] = yj[0] / vals[0];
        for (size_t row=1; row < numRows; ++row) {
          xj[row] = yj[row];
          size_t b = offsets[row], d = offsets[row+1]-1;
          for (size_t c=b; c != d; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
          xj[row] /= vals[d];
        }
      }
    }
  };


  // default implementation
  template <class Scalar, class Ordinal, class Node = DefaultNode::DefaultNodeType>
  class DefaultSparseSolve {
  public:
    typedef Scalar  ScalarType;
    typedef Ordinal OrdinalType;
    typedef Node    NodeType;

    //! @name Constructors/Destructor

    //@{

    //! DefaultSparseSolve constuctor with variable number of indices per row.
    DefaultSparseSolve(const Teuchos::RCP<Node> &node = DefaultNode::getDefaultNode());

    //! DefaultSparseSolve Destructor
    ~DefaultSparseSolve();

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
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const MultiVector<DomainScalar,Node> &Y, MultiVector<RangeScalar,Node> &X) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultSparseSolve(const DefaultSparseSolve& source);

    Teuchos::RCP<Node> node_;

    // we do this one of two ways: 
    // 1D/packed: array of offsets, pointer for ordinals, pointer for values. obviously the smallest footprint.
    Teuchos::ArrayRCP<const Ordinal> pbuf_inds1D_;
    Teuchos::ArrayRCP<const size_t>  pbuf_offsets1D_;
    Teuchos::ArrayRCP<const Scalar>  pbuf_vals1D_;

    size_t numRows_;
    bool indsInit_, valsInit_, isPacked_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseSolve<Scalar,Ordinal,Node>::DefaultSparseSolve(const Teuchos::RCP<Node> &node)
  : node_(node)
  , indsInit_(false)
  , valsInit_(false)
  , isPacked_(false) {
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseSolve<Scalar,Ordinal,Node>::~DefaultSparseSolve() {
  }

  template<class Scalar, class Ordinal, class Node>
  template <class GRAPH>
  Teuchos::DataAccess DefaultSparseSolve<Scalar,Ordinal,Node>::initializeStructure(GRAPH &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse graphs
    TEST_FOR_EXCEPT(true);
  }

  template<class Scalar, class Ordinal, class Node>
  template <class MATRIX>
  Teuchos::DataAccess DefaultSparseSolve<Scalar,Ordinal,Node>::initializeValues(MATRIX &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse matrices
    TEST_FOR_EXCEPT(true);
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::DataAccess DefaultSparseSolve<Scalar,Ordinal,Node>::initializeStructure(CrsGraph<Ordinal,Node> &graph, Teuchos::DataAccess cv) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(cv != Teuchos::View, std::runtime_error,
        Teuchos::typeName(*this) << "::initializeStructure(): requires View access.");
    TEST_FOR_EXCEPTION(indsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isPacked()) {
      isPacked_ = true;
      pbuf_inds1D_    = graph.get1DIndices();
      pbuf_offsets1D_ = graph.get1DOffsets();
    }
    else {
      // TEST_FOR_EXCEPT(true);
    }
    indsInit_ = true;
    return Teuchos::View;
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::DataAccess DefaultSparseSolve<Scalar,Ordinal,Node>::initializeValues(CrsMatrix<Scalar,Node> &matrix, Teuchos::DataAccess cv) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(cv != Teuchos::View, std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): requires View access.");
    TEST_FOR_EXCEPTION(valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): values already initialized.");
    TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isPacked_ != matrix.isPacked(), std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (matrix.isPacked()) {
      isPacked_ = true;
      pbuf_vals1D_ = matrix.get1DValues();
    }
    else {
      // TEST_FOR_EXCEPT(true);
    }
    valsInit_ = true;
    return Teuchos::View;
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::RCP<Node> DefaultSparseSolve<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultSparseSolve<Scalar,Ordinal,Node>::clear() { 
    pbuf_inds1D_     = Teuchos::null;
    pbuf_offsets1D_  = Teuchos::null;
    pbuf_vals1D_     = Teuchos::null;
    indsInit_ = false;
    valsInit_ = false;
    isPacked_ = false;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultSparseSolve<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                      MultiVector<RangeScalar,Node> &X) const {
    typedef DefaultSparseSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1D;
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): operation not fully initialized.");
    TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isPacked_ == true) {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.offsets = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRHS,wdp);
      }
      else {
        TEST_FOR_EXCEPT(true);
      }
    }
    else {
      TEST_FOR_EXCEPT(true);
    }
    return;
  }

} // namespace Kokkos

#endif /* KOKKOS_BASESPARSESOLVE_H */
