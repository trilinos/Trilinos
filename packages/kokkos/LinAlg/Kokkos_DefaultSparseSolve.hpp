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

#ifndef KOKKOS_BASESPARSESOLVE_H
#define KOKKOS_BASESPARSESOLVE_H

#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_BLAS_types.hpp>
#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_CrsGraph.hpp" 
#include "Kokkos_MultiVector.hpp"

#ifndef KERNEL_PREFIX
  #define KERNEL_PREFIX
#endif


namespace Kokkos {

  // 
  // Matrix formatting and mat-vec options
  // Applies to all four operations below
  // 
  // unitDiag indicates whether we neglect the diagonal row entry and scale by it
  // or utilize all row entries and implicitly scale by a unit diagonal (i.e., don't need to scale)
  // upper (versus lower) will determine the ordering of the solve and the location of the diagonal
  // 
  // upper -> diagonal is first entry on row
  // lower -> diagonal is last entry on row
  // 

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
      // 
      // upper triangular requires backwards substition, solving in reverse order
      // must unroll the last iteration, because decrement results in wrap-around
      // 
      if (upper && unitDiag) {
        // upper + unit
        xj[numRows-1] = yj[numRows-1];
        for (size_t r=2; r < numRows+1; ++r) {
          const size_t row = numRows - r; // for row=numRows-2 to 0 step -1
          const size_t begin = offsets[row], end = offsets[row+1];
          xj[row] = yj[row];
          for (size_t c=begin; c != end; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit
        xj[numRows-1] = yj[numRows-1] / vals[offsets[numRows-1]];
        for (size_t r=2; r < numRows+1; ++r) {
          const size_t row = numRows - r; // for row=numRows-2 to 0 step -1
          const size_t diag = offsets[row], end = offsets[row+1];
          const Scalar dval = vals[diag];
          xj[row] = yj[row];
          for (size_t c=diag+1; c != end; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
          xj[row] /= dval;
        }
      }
      else if (!upper && unitDiag) {
        // lower + unit
        xj[0] = yj[0];
        for (size_t row=1; row < numRows; ++row) {
          const size_t begin = offsets[row], end = offsets[row+1];
          xj[row] = yj[row];
          for (size_t c=begin; c != end; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit
        xj[0] = yj[0] / vals[0];
        for (size_t row=1; row < numRows; ++row) {
          const size_t begin = offsets[row], diag = offsets[row+1]-1;
          const Scalar dval = vals[diag];
          xj[row] = yj[row];
          for (size_t c=begin; c != diag; ++c) {
            xj[row] -= vals[c] * xj[inds[c]];
          }
          xj[row] /= dval;
        }
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseSolveOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    size_t numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar      *x;
    const RangeScalar *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // solve rhs i for lhs i
      const size_t rhs = i;
      DomainScalar      *xj = x + rhs * xstride;
      const RangeScalar *yj = y + rhs * ystride;
      const Scalar  *rowvals;
      const Ordinal *rowinds;
      Scalar dval;
      size_t nE;
      // 
      // upper triangular requires backwards substition, solving in reverse order
      // must unroll the last iteration, because decrement results in wrap-around
      // 
      if (upper && unitDiag) {
        // upper + unit
        xj[numRows-1] = yj[numRows-1];
        for (size_t row=numRows-2; row != 0; --row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          xj[row] = yj[row];
          for (size_t j=0; j != nE; ++j) {
            xj[row] -= rowvals[j] * xj[rowinds[j]];
          }
        }
        nE = numEntries[0];
        rowvals = vals_beg[0];
        rowinds = inds_beg[0];
        xj[0] = yj[0];
        for (size_t j=0; j != nE; ++j) {
          xj[0] -= rowvals[j] * xj[rowinds[j]];
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit: diagonal is first entry
        dval = vals_beg[numRows-1][0];
        xj[numRows-1] = yj[numRows-1] / dval;
        for (size_t row=numRows-2; row != 0; --row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          xj[row] = yj[row];
          Scalar dval = rowvals[0];
          for (size_t j=1; j < nE; ++j) {
            xj[row] -= rowvals[j] * xj[rowinds[j]];
          }
          xj[row] /= dval;
        }
        nE = numEntries[0];
        rowvals = vals_beg[0];
        rowinds = inds_beg[0];
        xj[0] = yj[0];
        Scalar dval = rowvals[0];
        for (size_t j=1; j < nE; ++j) {
          xj[0] -= rowvals[j] * xj[rowinds[j]];
        }
        xj[0] /= dval;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        xj[0] = yj[0];
        for (size_t row=1; row < numRows; ++row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          xj[row] = yj[row];
          for (size_t j=0; j < nE; ++j) {
            xj[row] -= rowvals[j] * xj[rowinds[j]];
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit; diagonal is last entry
        nE = numEntries[0];
        rowvals = vals_beg[0];
        dval = rowvals[0];
        xj[0] = yj[0];
        for (size_t row=1; row < numRows; ++row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          dval = rowvals[nE-1];
          xj[row] = yj[row];
          if (nE > 1) {
            for (size_t j=0; j < nE-1; ++j) {
              xj[row] -= rowvals[j] * xj[rowinds[j]];
            }
          }
          xj[row] /= dval;
        }
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseTransposeSolveOp1 {
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
      // 
      // put y into x and solve system in-situ
      // this is copy-safe, in the scenario that x and y point to the same location.
      //
      for (size_t row=0; row < numRows; ++row) {
        xj[row] = yj[row];
      }
      // 
      if (upper && unitDiag) {
        // upper + unit
        size_t beg, endplusone;
        for (size_t row=0; row < numRows-1; ++row) {
          beg = offsets[row]; 
          endplusone = offsets[row+1];
          for (size_t j=beg; j < endplusone; ++j) {
            xj[inds[j]] -= (vals[j] * xj[row]);
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit; diag is first element in row
        size_t diag, endplusone;
        Scalar dval;
        for (size_t row=0; row < numRows-1; ++row) {
          diag = offsets[row]; 
          endplusone = offsets[row+1];
          dval = vals[diag];
          xj[row] /= dval;
          for (size_t j=diag+1; j < endplusone; ++j) {
            xj[inds[j]] -= (vals[j] * xj[row]);
          }
        }
        diag = offsets[numRows-1];
        dval = vals[diag];
        xj[numRows-1] /= dval;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        for (size_t row=numRows-1; row > 0; --row) {
          size_t beg = offsets[row], endplusone = offsets[row+1];
          for (size_t j=beg; j < endplusone; ++j) {
            xj[inds[j]] -= (vals[j] * xj[row]);
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit; diag is last element in row
        Scalar dval;
        for (size_t row=numRows-1; row > 0; --row) {
          size_t beg = offsets[row], diag = offsets[row+1]-1;
          dval = vals[diag];
          xj[row] /= dval;
          for (size_t j=beg; j < diag; ++j) {
            xj[inds[j]] -= (vals[j] * xj[row]);
          }
        }
        // last row
        dval = offsets[0];
        xj[0] /= dval;
      }
    }
  };


  template <class Scalar, class Ordinal, class DomainScalar, class RangeScalar>
  struct DefaultSparseTransposeSolveOp2 {
    // mat data
    const Ordinal * const * inds_beg;
    const Scalar  * const * vals_beg;
    const size_t  *         numEntries;
    size_t numRows;
    // matvec params
    bool unitDiag, upper;
    // mv data
    DomainScalar      *x;
    const RangeScalar *y;
    size_t xstride, ystride;

    inline KERNEL_PREFIX void execute(size_t i) {
      // solve rhs i for lhs i
      const size_t rhs = i;
      DomainScalar      *xj = x + rhs * xstride;
      const RangeScalar *yj = y + rhs * ystride;
      const Scalar  *rowvals;
      const Ordinal *rowinds;
      Scalar dval;
      size_t nE;
      // 
      // put y into x and solve system in-situ
      // this is copy-safe, in the scenario that x and y point to the same location.
      //
      for (size_t row=0; row < numRows; ++row) {
        xj[row] = yj[row];
      }
      // 
      if (upper && unitDiag) {
        // upper + unit
        for (size_t row=0; row < numRows-1; ++row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          for (size_t j=0; j < nE; ++j) {
            xj[rowinds[j]] -= (rowvals[j] * xj[row]);
          }
        }
      }
      else if (upper && !unitDiag) {
        // upper + non-unit; diag is first element in row
        for (size_t row=0; row < numRows-1; ++row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          dval = rowvals[0];
          xj[row] /= dval;
          for (size_t j=1; j < nE; ++j) {
            xj[rowinds[j]] -= (rowvals[j] * xj[row]);
          }
        }
        rowvals = vals_beg[numRows-1];
        dval = rowvals[0];
        xj[numRows-1] /= dval;
      }
      else if (!upper && unitDiag) {
        // lower + unit
        for (size_t row=numRows-1; row > 0; --row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          for (size_t j=0; j < nE; ++j) {
            xj[rowinds[j]] -= (rowvals[j] * xj[row]);
          }
        }
      }
      else if (!upper && !unitDiag) {
        // lower + non-unit; diag is last element in row
        for (size_t row=numRows-1; row > 0; --row) {
          nE = numEntries[row];
          rowvals = vals_beg[row];
          rowinds = inds_beg[row];
          dval = rowvals[nE-1];
          xj[row] /= dval;
          for (size_t j=0; j < nE-1; ++j) {
            xj[rowinds[j]] -= (rowvals[j] * xj[row]);
          }
        }
        rowvals = vals_beg[0];
        dval = rowvals[0];
        xj[0] /= dval;
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
    // 2D: array of pointers
    Teuchos::ArrayRCP<const Ordinal *> pbuf_inds2D_;
    Teuchos::ArrayRCP<const Scalar  *> pbuf_vals2D_;
    Teuchos::ArrayRCP<size_t>          pbuf_numEntries_;

    size_t numRows_;
    bool indsInit_, valsInit_, isPacked_, isEmpty_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseSolve<Scalar,Ordinal,Node>::DefaultSparseSolve(const Teuchos::RCP<Node> &node)
  : node_(node)
  , indsInit_(false)
  , valsInit_(false)
  , isPacked_(false)
  , isEmpty_(false) {
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultSparseSolve<Scalar,Ordinal,Node>::~DefaultSparseSolve() {
  }

  template<class Scalar, class Ordinal, class Node>
  template <class GRAPH>
  Teuchos::DataAccess DefaultSparseSolve<Scalar,Ordinal,Node>::initializeStructure(GRAPH &graph, Teuchos::DataAccess cv) {
    // not implemented for general sparse graphs
    TEST_FOR_EXCEPTION(true, std::exception, 
        Teuchos::typeName(*this) << "::initializeStructure(): method is not implemented for graph of type " << Teuchos::typeName(graph));
  }

  template<class Scalar, class Ordinal, class Node>
  template <class MATRIX>
  Teuchos::DataAccess DefaultSparseSolve<Scalar,Ordinal,Node>::initializeValues(MATRIX &matrix, Teuchos::DataAccess cv) {
    // not implemented for general sparse matrices
    TEST_FOR_EXCEPTION(true, std::exception, 
        Teuchos::typeName(*this) << "::initializeValues(): method is not implemented for matrix of type " << Teuchos::typeName(matrix));
  }

  template <class Scalar, class Ordinal, class Node>
  Teuchos::DataAccess DefaultSparseSolve<Scalar,Ordinal,Node>::initializeStructure(CrsGraph<Ordinal,Node> &graph, Teuchos::DataAccess cv) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(cv != Teuchos::View, std::runtime_error,
        Teuchos::typeName(*this) << "::initializeStructure(): requires View access.");
    TEST_FOR_EXCEPTION(indsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    if (graph.isPacked()) {
      isEmpty_ = false;
      isPacked_ = true;
      pbuf_inds1D_    = graph.getPackedIndices();
      pbuf_offsets1D_ = graph.getPackedOffsets();
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
  Teuchos::DataAccess DefaultSparseSolve<Scalar,Ordinal,Node>::initializeValues(CrsMatrix<Scalar,Node> &matrix, Teuchos::DataAccess cv) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(cv != Teuchos::View, std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): requires View access.");
    TEST_FOR_EXCEPTION(valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): values already initialized.");
    TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isPacked_ != matrix.isPacked(), std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (isEmpty_ || matrix.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (matrix.isPacked()) {
      isEmpty_ = false;
      pbuf_vals1D_ = matrix.getPackedValues();
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
    isEmpty_ = false;
  }

  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultSparseSolve<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                      MultiVector<RangeScalar,Node> &X) const {
    typedef DefaultSparseSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1D;
    typedef DefaultSparseSolveOp2<Scalar,Ordinal,DomainScalar,RangeScalar>  Op2D;
    typedef DefaultSparseTransposeSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  TOp1D;
    typedef DefaultSparseTransposeSolveOp2<Scalar,Ordinal,DomainScalar,RangeScalar>  TOp2D;
    TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): this solve was not fully initialized.");
    TEST_FOR_EXCEPTION(X.getNumCols() != Y.getNumCols(), std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left hand side and right hand side multivectors have differing numbers of vectors.");
    TEST_FOR_EXCEPTION(X.getNumRows() < numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left-hand-side multivector does not have enough rows. Likely cause is that the column map was not provided to the Tpetra::CrsMatrix in the case of an implicit unit diagonal.");

    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEST_FOR_EXCEPTION(diag != Teuchos::UNIT_DIAG, std::runtime_error,
          Teuchos::typeName(*this) << "::solve(): solve of empty matrix only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else if (isPacked_ == true) {
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
        TOp1D wdp;
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
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op2D wdp;
        rbh.begin();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(pbuf_inds2D_);
        wdp.vals_beg   = rbh.template addConstBuffer<const Scalar *>(pbuf_vals2D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op2D>(0,numRHS,wdp);
      }
      else {
        TOp2D wdp;
        rbh.begin();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(pbuf_numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(pbuf_inds2D_);
        wdp.vals_beg   = rbh.template addConstBuffer<const Scalar *>(pbuf_vals2D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp2D>(0,numRHS,wdp);
      }
    }
    return;
  }

} // namespace Kokkos

#endif /* KOKKOS_BASESPARSESOLVE_H */
