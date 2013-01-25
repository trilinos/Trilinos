// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_CRSMATRIX_DEF_HPP
#define TPETRA_CRSMATRIX_DEF_HPP

// TODO: row-wise insertion of entries in globalAssemble() may be more efficient
// TODO: consider maintaining sorted entries at all times and leaning heavily on STL set_intersect/set_union methods for all insert/replace/suminto

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_as.hpp>

#include "Tpetra_CrsMatrixMultiplyOp.hpp" // must include for implicit instantiation to work
#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_CrsMatrix_decl.hpp"
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
//! Comparison operator for Tpetra::CrsIJV objects, used by Tpetra::CrsMatrix
template <class Ordinal, class Scalar>
bool std::operator<(const Tpetra::CrsIJV<Ordinal,Scalar> &ijv1, const Tpetra::CrsIJV<Ordinal,Scalar> &ijv2) {
  return ijv1.i < ijv2.i;
}
#endif

namespace Tpetra {
  //
  // Users must never rely on anything in the Details namespace.
  //
  namespace Details {
    /// \struct AbsMax
    /// \brief Functor that implements Tpetra::ABSMAX CombineMode.
    ///
    /// \tparam Scalar Same as the Scalar template parameter of CrsMatrix.
    template<class Scalar>
    struct AbsMax {
      Scalar operator() (const Scalar& x, const Scalar& y) {
        typedef Teuchos::ScalarTraits<Scalar> STS;

        return std::max (STS::magnitude (x), STS::magnitude (y));
      }
    };
  } // namespace Details

  template <class Ordinal, class Scalar>
  CrsIJV<Ordinal,Scalar>::CrsIJV() {}

  template <class Ordinal, class Scalar>
  CrsIJV<Ordinal,Scalar>::CrsIJV(Ordinal row, Ordinal col, const Scalar &val) {
    i = row;
    j = col;
    v = val;
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrix (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap,
             size_t maxNumEntriesPerRow,
             ProfileType pftype,
             const RCP<Teuchos::ParameterList>& params) :
    DistObject<char, LocalOrdinal, GlobalOrdinal, Node> (rowMap),
    insertGlobalValuesWarnedEfficiency_ (false),
    insertLocalValuesWarnedEfficiency_ (false)
  {
    try {
      myGraph_ = rcp (new Graph (rowMap, maxNumEntriesPerRow, pftype, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        typeName(*this) << "::CrsMatrix(): caught exception while allocating "
        "CrsGraph object: " << std::endl << e.what ());
    }
    staticGraph_ = myGraph_;
    // It is okay to create this now; this will prevent us from having
    // to check for it on every call to apply().  We will use a
    // nonowning RCP to wrap *this; this is safe as long as we do not
    // share sameScalarMultiplyOp_ with anyone.  Sharing it externally
    // would allow the object to persist past the CrsMatrix's
    // destruction.
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar> (rcp (this,false).getConst ());
    resumeFill(params);
    checkInternalState();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrix (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap,
             const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
             ProfileType pftype,
             const RCP<Teuchos::ParameterList>& params) :
    DistObject<char, LocalOrdinal, GlobalOrdinal, Node> (rowMap),
    insertGlobalValuesWarnedEfficiency_ (false),
    insertLocalValuesWarnedEfficiency_ (false)
  {
    try {
      myGraph_ = rcp (new Graph (rowMap, NumEntriesPerRowToAlloc, pftype, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: "
          << std::endl << e.what() << std::endl);
    }
    staticGraph_ = myGraph_;
    // See comment in constructor implementation above.
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar> (rcp (this,false).getConst ());
    resumeFill(params);
    checkInternalState();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrix (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
             const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
             size_t maxNumEntriesPerRow,
             ProfileType pftype,
             const RCP<Teuchos::ParameterList>& params) :
    DistObject<char, LocalOrdinal, GlobalOrdinal, Node> (rowMap),
    insertGlobalValuesWarnedEfficiency_ (false),
    insertLocalValuesWarnedEfficiency_ (false)
  {
    try {
      myGraph_ = rcp (new Graph (rowMap, colMap, maxNumEntriesPerRow, pftype, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        typeName(*this) << "::CrsMatrix(): caught exception while allocating "
        "CrsGraph object: " << std::endl << e.what ());
    }
    staticGraph_ = myGraph_;
    // See comment in constructor implementation above.
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar> (rcp (this,false).getConst ());
    resumeFill(params);
    checkInternalState();
  }

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrix (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
             const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
             const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
             ProfileType pftype,
             const RCP<Teuchos::ParameterList>& params) :
    DistObject<char, LocalOrdinal, GlobalOrdinal, Node> (rowMap),
    insertGlobalValuesWarnedEfficiency_ (false),
    insertLocalValuesWarnedEfficiency_ (false)
  {
    try {
      myGraph_ = rcp (new Graph (rowMap, colMap, NumEntriesPerRowToAlloc,
                                 pftype, params));
    }
    catch (std::exception &e) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        typeName(*this) << "::CrsMatrix(): caught exception while allocating "
        "CrsGraph object: " << std::endl << e.what ());
    }
    staticGraph_ = myGraph_;
    // See comment in constructor implementation above.
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar> (rcp (this,false).getConst ());
    resumeFill(params);
    checkInternalState();
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrix (const RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &graph,
             const RCP<Teuchos::ParameterList>& params)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node> (graph->getRowMap ())
  , staticGraph_ (graph)
  , insertGlobalValuesWarnedEfficiency_ (false)
  , insertLocalValuesWarnedEfficiency_ (false)
  {
    const std::string tfecfFuncName("CrsMatrix(graph)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(staticGraph_.is_null (),
      std::runtime_error, ": When calling the CrsMatrix constructor that "
      "accepts a static graph, the pointer to the graph must not be null.");
    // We prohibit the case where the graph is not yet filled.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! staticGraph_->isFillComplete (),
      std::runtime_error, ": The specified graph is not fill-complete. You "
      "must invoke fillComplete() on the graph before using it to construct a "
      "CrsMatrix.  Note that calling resumeFill() makes the graph not fill-"
      "complete, even if you had previously called fillComplete().  In that "
      "case, you must call fillComplete() on the graph again.");
    // See comment in constructor implementation above.
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );
    // the graph has entries, and the matrix should have entries as well, set to zero. no need or point in lazy allocating in this case.
    // first argument LocalIndices is ignored; the graph is already allocated (local or global, we don't care here)
    allocateValues (LocalIndices, GraphAlreadyAllocated);
    resumeFill(params);
    checkInternalState();
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  ~CrsMatrix() {}

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  const RCP<const Teuchos::Comm<int> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getComm() const {
    return getCrsGraph ()->getComm ();
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  RCP<Node>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNode() const {
    return staticGraph_->getNode ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ProfileType CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getProfileType() const {
    return getCrsGraph()->getProfileType();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isFillComplete() const {
    return fillComplete_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isFillActive() const {
    return !fillComplete_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isStorageOptimized() const {
    return getCrsGraph()->isStorageOptimized();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isLocallyIndexed() const {
    return getCrsGraph()->isLocallyIndexed();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isGloballyIndexed() const {
    return getCrsGraph()->isGloballyIndexed();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasColMap() const {
    return getCrsGraph()->hasColMap();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalNumEntries() const {
    return getCrsGraph()->getGlobalNumEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeNumEntries() const {
    return getCrsGraph()->getNodeNumEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalNumRows() const {
    return getCrsGraph()->getGlobalNumRows();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalNumCols() const {
    return getCrsGraph()->getGlobalNumCols();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeNumRows() const {
    return getCrsGraph()->getNodeNumRows();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeNumCols() const {
    return getCrsGraph()->getNodeNumCols();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalNumDiags() const {
    return getCrsGraph()->getGlobalNumDiags();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeNumDiags() const {
    return getCrsGraph()->getNodeNumDiags();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
    return getCrsGraph()->getNumEntriesInGlobalRow(globalRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    return getCrsGraph()->getNumEntriesInLocalRow(localRow);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalMaxNumRowEntries() const {
    return getCrsGraph()->getGlobalMaxNumRowEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNodeMaxNumRowEntries() const {
    return getCrsGraph()->getNodeMaxNumRowEntries();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getIndexBase() const {
    return getRowMap()->getIndexBase();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRowMap() const {
    return getCrsGraph()->getRowMap();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getColMap() const {
    return getCrsGraph()->getColMap();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getDomainMap() const {
    return getCrsGraph()->getDomainMap();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRangeMap() const {
    return getCrsGraph()->getRangeMap();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGraph() const {
    if (staticGraph_ != null) return staticGraph_;
    return myGraph_;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getCrsGraph() const {
    if (staticGraph_ != null) return staticGraph_;
    return myGraph_;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isLowerTriangular() const {
    return getCrsGraph()->isLowerTriangular();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isUpperTriangular() const {
    return getCrsGraph()->isUpperTriangular();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isStaticGraph() const {
    return (myGraph_ == null);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasTransposeApply() const {
    return true;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::supportsRowViews() const {
    return true;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                    Internal utility methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  allocateValues (ELocalGlobal lg, GraphAllocationStatus gas)
  {
#ifdef HAVE_TPETRA_DEBUG
    // If the graph indices are already allocated, then gas should be
    // GraphAlreadyAllocated.  Otherwise, gas should be
    // GraphNotYetAllocated.
    if ((gas == GraphAlreadyAllocated) != staticGraph_->indicesAreAllocated()) {
      const std::string err1 ("allocateValues: The caller has asserted that "
                              "the graph is ");
      const std::string err2 ("already allocated, but the static graph says "
                              "that its indices are ");
      const std::string err3 ("already allocated.  Please report this bug to "
                              "the Tpetra developers.");
      TEUCHOS_TEST_FOR_EXCEPTION(gas == GraphAlreadyAllocated && ! staticGraph_->indicesAreAllocated(),
        std::logic_error, err1 << err2 << "not " << err3);
      TEUCHOS_TEST_FOR_EXCEPTION(gas != GraphAlreadyAllocated && staticGraph_->indicesAreAllocated(),
        std::logic_error, err1 << "not " << err2 << err3);
    }

    // If the graph is unallocated, then it had better be a
    // matrix-owned graph.  ("Matrix-owned graph" means that the
    // matrix gets to define the graph structure.  If the CrsMatrix
    // constructor that takes an RCP<const CrsGraph> was used, then
    // the matrix does _not_ own the graph.)
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! staticGraph_->indicesAreAllocated() && myGraph_.is_null(),
      std::logic_error,
      "allocateValues: The static graph says that its indices are not "
      "allocated, but the graph is not owned by the matrix.  Please report "
      "this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    if (gas == GraphNotYetAllocated) {
      myGraph_->allocateIndices (lg);
    }
    // ask graph to allocate our values, with the same structure
    // this will allocate values2D_ one way or the other
    if (getProfileType() == StaticProfile) {
      // "Static profile" means that the number of matrix entries in
      // each row was fixed at the time the CrsMatrix constructor was
      // called.  This lets us use 1-D storage for the matrix's
      // values.  ("1-D storage" means the same as that used by the
      // three arrays in the classic compressed sparse row format.)
      values1D_ = staticGraph_->template allocateValues1D<Scalar>();
    }
    else {
      // "Dynamic profile" means the number of matrix entries in each
      // row is not fixed and may expand.  Thus, we store the matrix's
      // values in "2-D storage," meaning an array of arrays.  The
      // outer array has as many inner arrays as there are rows in the
      // matrix, and each inner array stores the values in that row.
      values2D_ = staticGraph_->template allocateValues2D<Scalar>();
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  fillLocalGraphAndMatrix (const RCP<ParameterList> &params)
  {
    typedef LocalOrdinal LO;

    const size_t numRows = getNodeNumRows();
    // This method's goal is to fill in the three arrays (compressed
    // sparse row format) that define the sparse graph's and matrix's
    // structure, and the sparse matrix's values.
    ArrayRCP<LO>     inds;
    ArrayRCP<size_t> ptrs;
    ArrayRCP<Scalar> vals;

    // Get references to the data in myGraph_, so we can modify them
    // as well.  Note that we only call fillLocalGraphAndMatrix() if
    // the matrix owns the graph, which means myGraph_ is not null.
    ArrayRCP<LO>            &lclInds1D_     = myGraph_->lclInds1D_;
    ArrayRCP<ArrayRCP<LO> > &lclInds2D_     = myGraph_->lclInds2D_;
    ArrayRCP<size_t>        &rowPtrs_       = myGraph_->rowPtrs_;
    ArrayRCP<size_t>        &numRowEntries_ = myGraph_->numRowEntries_;
    size_t & nodeNumEntries_   = myGraph_->nodeNumEntries_;
    size_t & nodeNumAllocated_ = myGraph_->nodeNumAllocated_;

    if (getProfileType() == DynamicProfile) {
      // DynamicProfile means that the matrix's column indices and
      // values are currently stored in a 2-D "unpacked" format, in
      // the arrays-of-arrays lclInds2D_ (for column indices) and
      // values2D_ (for values).  We allocate 1-D storage and then
      // copy from 2-D storage in lclInds2D_ resp. values2D_ into 1-D
      // storage in inds resp. vals.
      ptrs = sparse_ops_type::allocRowPtrs (getRowMap ()->getNode (), numRowEntries_ ());
      inds = sparse_ops_type::template allocStorage<LO> (getRowMap ()->getNode (), ptrs ());
      vals = sparse_ops_type::template allocStorage<Scalar> (getRowMap ()->getNode (), ptrs ());

      // numRowEntries_ tells the number of valid entries
      // in each row (as opposed to the allocated size)
      for (size_t row=0; row < numRows; ++row) {
        const size_t numentrs = numRowEntries_[row];
        std::copy (lclInds2D_[row].begin(),
                   lclInds2D_[row].begin() + numentrs,
                   inds + ptrs[row]);
        std::copy (values2D_[row].begin(),
                   values2D_[row].begin() + numentrs,
                   vals + ptrs[row]);
      }
    }
    else if (getProfileType() == StaticProfile) {
      // StaticProfile means that the matrix's column indices and
      // values are currently stored in a 1-D format.  However, this
      // format is "unpacked"; it doesn't necessarily have the same
      // row offsets as indicated by the ptrs array returned by
      // allocRowPtrs.  This could happen, for example, if the user
      // specified StaticProfile in the constructor and fixed the
      // number of matrix entries in each row, but didn't fill all
      // those entries.
      if (nodeNumEntries_ != nodeNumAllocated_) {
        // We have to pack the 1-D storage, since the user didn't fill
        // up all requested storage.  We compute the row offsets
        // (ptrs) from numRowEntries_, which has the true number of
        // inserted entries in each row (vs. the number that we
        // requested when constructing the matrix, which is used by
        // the unpacked row offsets array rowPtrs_).
        ptrs = sparse_ops_type::allocRowPtrs (getRowMap ()->getNode (), numRowEntries_ ());
        inds = sparse_ops_type::template allocStorage<LO> (getRowMap ()->getNode (), ptrs ());
        vals = sparse_ops_type::template allocStorage<Scalar> (getRowMap ()->getNode (), ptrs ());
        for (size_t row=0; row < numRows; ++row) {
          // rowPtrs_ contains the unpacked row offsets, so use it to
          // copy data out of unpacked 1-D storage.
          const size_t numentrs = numRowEntries_[row];
          std::copy (lclInds1D_.begin() + rowPtrs_[row],
                     lclInds1D_.begin() + rowPtrs_[row] + numentrs,
                     inds + ptrs[row]);
          std::copy (values1D_.begin() + rowPtrs_[row],
                     values1D_.begin() + rowPtrs_[row] + numentrs,
                     vals + ptrs[row]);
        }
      }
      else {
        // The user filled up all requested storage, so we don't have
        // to pack.
        ptrs = rowPtrs_;
        inds = lclInds1D_;
        vals = values1D_;
      }
    }

    // May we ditch the old allocations for the packed (and otherwise
    // "optimized") allocations, later in this routine?  Optimize
    // storage if the graph is not static, or if the graph already has
    // optimized storage.
    const bool default_OptimizeStorage =
      ! isStaticGraph () || staticGraph_->isStorageOptimized ();
    const bool requestOptimizedStorage =
      (params != null && params->get ("Optimize Storage", default_OptimizeStorage))
      ||
      (params == null && default_OptimizeStorage);

    // The graph has optimized storage when indices are allocated,
    // numRowEntries_ is null, and there are more than zero rows on
    // this process.  It's impossible for the graph to have dynamic
    // profile (getProfileType() == DynamicProfile) and be optimized
    // (isStorageOptimized()).
    if (requestOptimizedStorage) {
      // Free the old, unpacked, unoptimized allocations.
      // Change the graph from dynamic to static allocaiton profile
      //
      // delete old data
      lclInds2D_ = null;
      numRowEntries_ = null;
      values2D_ = null;
      // keep the new, packed, optimized allocations
      lclInds1D_ = inds;
      rowPtrs_   = ptrs;
      values1D_  = vals;
      // we're packed: number allocated is the same as the number of valid entries
      nodeNumAllocated_ = nodeNumEntries_;
      // we've switched to a static allocation profile (if we weren't already)
      myGraph_->pftype_ = StaticProfile;
    }

    RCP<ParameterList> lclparams;
    if (params == null) {
      lclparams = parameterList ();
    }
    else {
      lclparams = sublist (params, "Local Graph");
    }

    // Make the local graph, using the ptrs and inds arrays we build
    // above.  The local graph should be null, but we delete it first
    // so that any memory can be freed before we allocate the new one.
    myGraph_->lclGraph_ = null;
    myGraph_->lclGraph_ =
      rcp (new local_graph_type (getRowMap ()->getNodeNumElements (),
                                 getColMap ()->getNodeNumElements (),
                                 getRowMap ()->getNode (), lclparams));
    myGraph_->lclGraph_->setStructure (ptrs, inds);

    // Now the graph has ptrs and inds, so we don't need to keep them here.
    ptrs = null;
    inds = null;

    // Make the local matrix, using the local graph and vals array.
    if (params == null) {
      lclparams = parameterList ();
    }
    else {
      lclparams = sublist (params, "Local Matrix");
    }
    // The local matrix should be null, but we delete it first so that
    // any memory can be freed before we allocate the new one.
    lclMatrix_ = null;
    lclMatrix_ = rcp (new local_matrix_type (staticGraph_->getLocalGraph (), lclparams));
    lclMatrix_->setValues (vals);
    // Now the matrix has vals, so we don't need to keep it here.
    vals = null;

    // Finalize the local graph and matrix together.
    if (params == null) {
      lclparams = parameterList ();
    }
    else {
      lclparams = sublist (params, "Local Sparse Ops");
    }
    // Figure out if the matrix has a unit diagonal, and whether it is
    // upper or lower triangular (or neither).
    const Teuchos::EDiag diag = getNodeNumDiags() < getNodeNumRows() ?
      Teuchos::UNIT_DIAG :
      Teuchos::NON_UNIT_DIAG;
    Teuchos::EUplo uplo = Teuchos::UNDEF_TRI;
    if (isUpperTriangular ()) {
      uplo = Teuchos::UPPER_TRI;
    }
    else if (isLowerTriangular ()) {
      uplo = Teuchos::LOWER_TRI;
    }
    sparse_ops_type::finalizeGraphAndMatrix (uplo, diag,
                                             *myGraph_->getLocalGraphNonConst (),
                                             *lclMatrix_,
                                             lclparams);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  fillLocalMatrix (const RCP<ParameterList> &params)
  {
    const size_t numRows = getNodeNumRows();

    // The goals of this routine are first, to allocate and fill
    // packed 1-D storage (see below for an explanation) in the vals
    // array, and second, to give vals to the local matrix and
    // finalize the local matrix.  We only need ptrs, the packed 1-D
    // row offsets, within the scope of this routine, since we're only
    // filling the local matrix here (use fillLocalGraphAndMatrix() to
    // fill both the graph and the matrix at the same time).
    ArrayRCP<size_t> ptrs;
    ArrayRCP<Scalar> vals;

    // get data from staticGraph_
    ArrayRCP<LocalOrdinal>            lclInds1D     = staticGraph_->lclInds1D_;
    ArrayRCP<ArrayRCP<LocalOrdinal> > lclInds2D     = staticGraph_->lclInds2D_;
    ArrayRCP<size_t>                  rowPtrs       = staticGraph_->rowPtrs_;
    ArrayRCP<size_t>                  numRowEntries = staticGraph_->numRowEntries_;
    size_t nodeNumEntries   = staticGraph_->nodeNumEntries_;
    size_t nodeNumAllocated = staticGraph_->nodeNumAllocated_;

    // May we ditch the old allocations for the packed (and otherwise
    // "optimized") allocations, later in this routine?  Request
    // optimized storage by default.
    bool requestOptimizedStorage = true;
    const bool default_OptimizeStorage =
      ! isStaticGraph () || staticGraph_->isStorageOptimized ();
    if (params != null && ! params->get ("Optimize Storage", default_OptimizeStorage)) {
      requestOptimizedStorage = false;
    }
    // If we're not allowed to change a static graph, then we can't
    // change the storage of the matrix, either.  This means that if
    // the graph's storage isn't already optimized, we can't optimize
    // the matrix's storage either.  Check and give warning, as
    // appropriate.
    if (staticGraph_->isStorageOptimized() == false && requestOptimizedStorage)
    {
      TPETRA_ABUSE_WARNING(true, std::runtime_error,
        "::fillLocalMatrix(): You requested optimized storage by setting the"
        "\"Optimize Storage\" flag to \"true\" in the parameter list, or by virtue"
        "of default behavior. However, the associated CrsGraph was filled separately"
        "and requested not to optimize storage. Therefore, the CrsMatrix cannot"
        "optimize storage.");
      requestOptimizedStorage = false;
    }

    if (getProfileType() == DynamicProfile) {
      // DynamicProfile means that the matrix's values are currently
      // stored in a 2-D "unpacked" format, in the array-of-arrays
      // values2D_.  We allocate 1-D storage and then copy from 2-D
      // storage in values2D_ into 1-D storage in vals.  Since we're
      // only allocating the local matrix here, not the local graph,
      // we don't need to keep the row offsets array ptrs, but we do
      // need it here temporarily in order to convert to 1-D storage.
      // (The allocStorage() function needs it.)  We'll free ptrs
      // later in this method.
      ptrs = sparse_ops_type::allocRowPtrs (getRowMap ()->getNode (), numRowEntries ());
      vals = sparse_ops_type::template allocStorage<Scalar> (getRowMap ()->getNode (), ptrs ());
      // TODO (mfh 05 Dec 2012) We should really parallelize this copy
      // operation.  This is not currently required in the
      // sparse_ops_type interface.  Some implementations of that
      // interface (such as AltSparseOps) do provide a copyStorage()
      // method, but this does not currently cover the case of copying
      // from 2-D to 1-D storage.
      for (size_t row=0; row < numRows; ++row) {
        const size_t numentrs = numRowEntries[row];
        std::copy (values2D_[row].begin(),
                   values2D_[row].begin() + numentrs,
                   vals+ptrs[row]);
      }
    }
    else if (getProfileType() == StaticProfile) {
      // StaticProfile means that the matrix's values are currently
      // stored in a 1-D format.  However, this format is "unpacked";
      // it doesn't necessarily have the same row offsets as indicated
      // by the ptrs array returned by allocRowPtrs.  This could
      // happen, for example, if the user specified StaticProfile in
      // the constructor and fixed the number of matrix entries in
      // each row, but didn't fill all those entries.
      //
      // As above, we don't need to keep the "packed" row offsets
      // array ptrs here, but we do need it here temporarily, so we
      // have to allocate it.  We'll free ptrs later in this method.
      //
      // Note that this routine checks whether storage has already
      // been packed.  This is a common case for solution of nonlinear
      // PDEs using the finite element method, as long as the
      // structure of the sparse matrix does not change between linear
      // solves.
      if (nodeNumEntries != nodeNumAllocated) {
        // We have to pack the 1-D storage, since the user didn't fill
        // up all requested storage.
        ptrs = sparse_ops_type::allocRowPtrs (getRowMap ()->getNode (), numRowEntries ());
        vals = sparse_ops_type::template allocStorage<Scalar> (getRowMap ()->getNode (), ptrs ());
        // TODO (mfh 05 Dec 2012) We should really parallelize this
        // copy operation.  This is not currently required in the
        // sparse_ops_type interface.  Some implementations of that
        // interface (such as AltSparseOps) do provide a copyStorage()
        // method, but I have to check whether it requires that the
        // input have the same packed offsets as the output.
        for (size_t row=0; row < numRows; ++row) {
          const size_t numentrs = numRowEntries[row];
          std::copy (values1D_.begin() + rowPtrs[row],
                     values1D_.begin() + rowPtrs[row]+numentrs,
                     vals + ptrs[row]);
        }
      }
      else {
        // The user filled up all requested storage, so we don't have
        // to pack.
        vals = values1D_;
      }
    }
    // We're done with the packed row offsets array now.
    ptrs = null;

    // May we ditch the old allocations for the packed one?
    if (requestOptimizedStorage) {
      // The user requested optimized storage, so we can dump the
      // unpacked 2-D and 1-D storage, and keep the packed storage.
      values2D_ = null;
      values1D_ = vals;
    }

    // build the matrix, hand over the values
    RCP<ParameterList> lclparams;
    if (params == null) {
      lclparams = parameterList ();
    }
    else {
      lclparams = sublist (params, "Local Matrix");
    }

    // The local matrix should be null at this point Just in case it isn't
    // (future-proofing), delete it first in order to free memory before we allocate a new
    // one.  Otherwise, we risk storing two matrices temporarily, since the destructor
    // of the old matrix won't be called until the new matrix's constructor finishes.
    lclMatrix_ = null;
    lclMatrix_ = rcp (new local_matrix_type (staticGraph_->getLocalGraph (), lclparams));
    lclMatrix_->setValues (vals);
    vals = null;

    // Finalize the local matrix.
    if (params == null) {
      lclparams = parameterList ();
    }
    else {
      lclparams = sublist (params, "Local Sparse Ops");
    }

    // mfh 05 Dec 2012: This is the place where the matrix's data
    // might get reorganized into a different data structure, to match
    // the data structure of the graph.  This is not necessarily the
    // optimized final data structure (as used by apply() for sparse
    // matrix-vector multiply).  That happens at the following line in fillComplete():
    //
    // lclMatOps_->setGraphAndMatrix (staticGraph_->getLocalGraph (), lclMatrix_);
    //
    // The requirement to allow the graph to be stored separately from
    // the matrix is important for some applications.  It may save
    // memory if multiple matrices share the same structure, and it
    // allows the graph (and therefore the storage layout of the
    // matrix's values) to be precomputed.
    sparse_ops_type::finalizeMatrix (*staticGraph_->getLocalGraph (), *lclMatrix_, lclparams);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                  User-visible class methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  insertLocalValues (LocalOrdinal localRow,
                     const ArrayView<const LocalOrdinal> &indices,
                     const ArrayView<const Scalar>       &values)
  {
    const std::string tfecfFuncName("insertLocalValues()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillActive (), std::runtime_error,
      " requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isStaticGraph (),  std::runtime_error,
      " cannot insert indices with static graph; use replaceLocalValues() instead.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(myGraph_->isGloballyIndexed(),
      std::runtime_error, ": graph indices are global; use insertGlobalValues().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! hasColMap (), std::runtime_error,
      " cannot insert local indices without a column map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(),
      std::runtime_error, ": values.size() must equal indices.size().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! getRowMap()->isNodeLocalElement(localRow), std::runtime_error,
      ": Local row index " << localRow << " does not belong to this process.");
    if (! myGraph_->indicesAreAllocated ()) {
      allocateValues (LocalIndices, GraphNotYetAllocated);
    }
    // use column map to filter the entries:
    Array<LocalOrdinal> f_inds(indices);
    Array<Scalar>       f_vals(values);
    typename Graph::SLocalGlobalNCViews inds_ncview;
    inds_ncview.linds = f_inds();
    const size_t numFilteredEntries =
      myGraph_->template filterIndicesAndValues<LocalIndices, Scalar> (inds_ncview,
                                                                       f_vals ());
    if (numFilteredEntries > 0) {
      RowInfo rowInfo = myGraph_->getRowInfo(localRow);
      const size_t curNumEntries = rowInfo.numEntries;
      const size_t newNumEntries = curNumEntries + numFilteredEntries;
      if (newNumEntries > rowInfo.allocSize) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(getProfileType() == StaticProfile,
          std::runtime_error, ": new indices exceed statically allocated graph "
          "structure.");

        // Only print an efficiency warning once per CrsMatrix
        // instance, per method name (insertLocalValues() or
        // insertGlobalValues()).
        if (! insertLocalValuesWarnedEfficiency_) {
          TPETRA_EFFICIENCY_WARNING(true, std::runtime_error,
            "::insertLocalValues():" << std::endl
            << "Pre-allocated space has been exceeded, requiring new allocation.  "
            "This is allowed but not efficient in terms of run time.  "
            "To improve efficiency, suggest a larger number of entries per row in the constructor.  "
            "You may either specify a maximum number of entries for all the rows, or a per-row maximum.  "
            "This CrsMatrix instance will not print further messages of this kind, in order not to clutter output.");
          insertLocalValuesWarnedEfficiency_ = true;
        }
        // Make space for the new matrix entries.
        rowInfo = myGraph_->template updateAllocAndValues<LocalIndices, Scalar> (rowInfo, newNumEntries, values2D_[localRow]);
      }
      typename Graph::SLocalGlobalViews inds_view;
      inds_view.linds = f_inds(0,numFilteredEntries);
      myGraph_->template insertIndicesAndValues<LocalIndices, LocalIndices> (rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), f_vals.begin());
#ifdef HAVE_TPETRA_DEBUG
      {
        const size_t chkNewNumEntries = myGraph_->getNumEntriesInLocalRow (localRow);
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(chkNewNumEntries != newNumEntries,
          std::logic_error, ": Internal logic error. Please contact Tpetra team.");
      }
#endif
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isLocallyIndexed(), std::logic_error,
      ": At end of insertLocalValues(), this CrsMatrix is not locally indexed.  "
      "Please report this bug to the Tpetra developers.");
#endif
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  insertGlobalValues (GlobalOrdinal globalRow,
                      const ArrayView<const GlobalOrdinal> &indices,
                      const ArrayView<const Scalar>        &values)
  {
    const std::string tfecfFuncName("insertGlobalValues()");

    // mfh 14 Dec 2012: Defer test for static graph until we know that
    // globalRow is in the row Map.  If it's not in the row Map, it
    // doesn't matter whether or not the graph is static; the data
    // just go into nonlocals_, for later use by globalAssemble().
    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    //   isStaticGraph(), std::runtime_error,
    //   ": matrix was constructed with static graph. Cannot insert new entries.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      values.size() != indices.size(), std::runtime_error,
      ": values.size() must equal indices.size().  values.size() = "
      << values.size() << ", but indices.size() = " << indices.size() << ".");

    const LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);

    if (lrow != LOT::invalid()) { // globalRow is in our row Map.
      // If the matrix has a static graph, this process is now allowed
      // to insert into rows it owns.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        this->isStaticGraph(), std::runtime_error,
        ": The CrsMatrix was constructed with a static graph.  In that case, "
        "it's forbidded to insert new entries into rows owned by the calling process.");
      if (! myGraph_->indicesAreAllocated ()) {
        allocateValues (GlobalIndices, GraphNotYetAllocated);
      }
      typename Graph::SLocalGlobalViews inds_view;
      ArrayView<const Scalar> vals_view;

      // We have to declare these Arrays here rather than in the
      // hasColMap() if branch, so that views to them will remain
      // valid for the whole scope.
      Array<GlobalOrdinal> filtered_indices;
      Array<Scalar>        filtered_values;
      if (hasColMap()) { // We have a column Map.
        //
        // Use column map to filter the indices and corresponding
        // values, so that we only insert entries into columns we own.
        typename Graph::SLocalGlobalNCViews inds_ncview;
        ArrayView<Scalar> vals_ncview;
        filtered_indices.assign(indices.begin(), indices.end());
        filtered_values.assign(values.begin(), values.end());
        inds_ncview.ginds = filtered_indices();
        const size_t numFilteredEntries = myGraph_->template filterIndicesAndValues<GlobalIndices,Scalar>(inds_ncview,filtered_values());
        inds_view.ginds = filtered_indices(0,numFilteredEntries);
        vals_view       = filtered_values(0,numFilteredEntries);
      }
      else { // we don't have a column Map.
        inds_view.ginds = indices;
        vals_view       = values;
      }
      const size_t numFilteredEntries = vals_view.size();
      // add the new indices and values
      if (numFilteredEntries > 0) {
        RowInfo rowInfo = myGraph_->getRowInfo(lrow);
        const size_t curNumEntries = rowInfo.numEntries;
        const size_t newNumEntries = curNumEntries + numFilteredEntries;
        if (newNumEntries > rowInfo.allocSize) {
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            getProfileType() == StaticProfile, std::runtime_error,
            ": new indices exceed statically allocated graph structure.");

          // Only print an efficiency warning once per CrsMatrix
          // instance, per method name (insertLocalValues() or
          // insertGlobalValues()).
          if (! insertGlobalValuesWarnedEfficiency_) {
            TPETRA_EFFICIENCY_WARNING(true, std::runtime_error,
              "::insertGlobalValues():" << std::endl
              << "Pre-allocated space has been exceeded, requiring new allocation.  "
              "This is allowed but not efficient in terms of run time.  "
              "To improve efficiency, suggest a larger number of entries per row in the constructor.  "
              "You may either specify a maximum number of entries for all the rows, or a per-row maximum.  "
              "This CrsMatrix instance will not print further messages of this kind, in order not to clutter output.");
            insertGlobalValuesWarnedEfficiency_ = true;
          }
          // Update allocation only as much as necessary
          rowInfo = myGraph_->template updateAllocAndValues<GlobalIndices,Scalar>(rowInfo, newNumEntries, values2D_[lrow]);
        }
        if (isGloballyIndexed()) {
          // <GlobalIndices, GlobalIndices> template parameters
          // involve getGlobalViewNonConst() and direct copying, which
          // should be reasonably fast.
          myGraph_->template insertIndicesAndValues<GlobalIndices,GlobalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), vals_view.begin());
        }
        else {
          // <GlobalIndices, LocalIndices> template parameters involve
          // calling the Map's getLocalElement() once per entry to
          // insert.  This may be slow.
          myGraph_->template insertIndicesAndValues<GlobalIndices,LocalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), vals_view.begin());
        }
#ifdef HAVE_TPETRA_DEBUG
        {
          const size_t chkNewNumEntries = myGraph_->getNumEntriesInLocalRow(lrow);
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(chkNewNumEntries != newNumEntries,
            std::logic_error, ": There should be a total of " << newNumEntries
            << " entries in the row, but the graph now reports " << chkNewNumEntries
            << " entries.  Please report this bug to the Tpetra developers.");
        }
#endif // HAVE_TPETRA_DEBUG
      }
    }
    else { // The calling process doesn't own the given row, so add
           // the new data to the list of nonlocals.
      typename ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
      typename ArrayView<const Scalar       >::iterator val =  values.begin();
      nonlocals_[globalRow].reserve( nonlocals_[globalRow].size() + indices.size() );
      for (; val != values.end(); ++val, ++ind) {
        nonlocals_[globalRow].push_back(std::make_pair(*ind, *val));
      }
    }
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  replaceLocalValues (LocalOrdinal localRow,
                      const ArrayView<const LocalOrdinal> &indices,
                      const ArrayView<const Scalar> &values)
  {
    // secondArg is a binary function that returns its second
    // argument.  This replaces entries in the given row with their
    // corresponding entry of values.
    typedef secondArg<Scalar, Scalar> f_type;
    this->template transformLocalValues<f_type> (localRow, indices,
                                                 values, f_type ());
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  replaceGlobalValues (GlobalOrdinal globalRow,
                       const ArrayView<const GlobalOrdinal> &indices,
                       const ArrayView<const Scalar>        &values)
  {
    // secondArg is a binary function that returns its second
    // argument.  This replaces entries in the given row with their
    // corresponding entry of values.
    typedef secondArg<Scalar, Scalar> f_type;
    this->template transformGlobalValues<f_type> (globalRow, indices,
                                                  values, f_type ());
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  sumIntoGlobalValues (GlobalOrdinal globalRow,
                       const ArrayView<const GlobalOrdinal> &indices,
                       const ArrayView<const Scalar>        &values)

  {
    try {
      typedef std::plus<Scalar> f_type;
      this->template transformGlobalValues<f_type> (globalRow, indices,
                                                    values, f_type ());
    }
    catch (Details::InvalidGlobalRowIndex<GlobalOrdinal>& e) {
      // For nonlocal data, use insertGlobalValues().  Since globalRow
      // is not in the row Map, this will put the data in nonlocals_.
      // globalAssemble() (called during fillComplete()) will exchange
      // that data and sum it in using sumIntoGlobalValues().
      this->insertGlobalValues (globalRow, indices, values);
    }
  }


  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node,
            class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  sumIntoLocalValues (LocalOrdinal localRow,
                      const ArrayView<const LocalOrdinal>  &indices,
                      const ArrayView<const Scalar>        &values)
  {
    typedef std::plus<Scalar> f_type;
    this->template transformLocalValues<f_type> (localRow, indices,
                                                 values, f_type ());
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayView<const Scalar>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  getView (RowInfo rowinfo) const
  {
    ArrayView<const Scalar> view;
    if (values1D_ != null && rowinfo.allocSize > 0) {
      view = values1D_(rowinfo.offset1D,rowinfo.allocSize);
    }
    else if (values2D_ != null) {
      view = values2D_[rowinfo.localRow]();
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayView<Scalar> CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getViewNonConst(RowInfo rowinfo)
  {
    ArrayView<Scalar> view;
    if (values1D_ != null && rowinfo.allocSize > 0) {
      view = values1D_(rowinfo.offset1D,rowinfo.allocSize);
    }
    else if (values2D_ != null) {
      view = values2D_[rowinfo.localRow]();
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalRowCopy(
                                LocalOrdinal localRow,
                                const ArrayView<LocalOrdinal> &indices,
                                const ArrayView<Scalar>       &values,
                                size_t &numEntries) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      isGloballyIndexed() && ! hasColMap(),
      std::runtime_error,
      "Tpetra::CrsMatrix::getLocalRowCopy: The matrix is globally indexed and "
      "does not have a column Map yet.  That means we don't have local indices "
      "for columns yet, so it doesn't make sense to call this method.  If the "
      "matrix doesn't have a column Map yet, you should call fillComplete on "
      "it first.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! getRowMap ()->isNodeLocalElement (localRow), 
      std::runtime_error,
      "Tpetra::CrsMatrix::getLocalRowCopy: The calling process " 
      << this->getComm ()->getRank () << " does not own the given local row " 
      << localRow << ".");
    const RowInfo rowinfo = staticGraph_->getRowInfo(localRow);
    numEntries = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION(
      static_cast<size_t>(indices.size()) < numEntries || 
      static_cast<size_t>(values.size()) < numEntries,
      std::runtime_error, 
      "Tpetra::CrsMatrix::getLocalRowCopy: The given row " << localRow 
      << " has " << numEntries << " entries.  One or both of the given array "
      "views are not long enough to store that many entries.  indices can "
      "store " << indices.size() << " entries and values can store " 
      << values.size() << " entries.");
    if (staticGraph_->isLocallyIndexed()) {
      ArrayView<const LocalOrdinal> indrowview = staticGraph_->getLocalView(rowinfo);
      ArrayView<const Scalar>       valrowview = getView(rowinfo);
      std::copy( indrowview.begin(), indrowview.begin() + numEntries, indices.begin() );
      std::copy( valrowview.begin(), valrowview.begin() + numEntries,  values.begin() );
    }
    else if (staticGraph_->isGloballyIndexed()) {
      ArrayView<const GlobalOrdinal> indrowview = staticGraph_->getGlobalView(rowinfo);
      ArrayView<const Scalar>        valrowview = getView(rowinfo);
      std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      for (size_t j=0; j < numEntries; ++j) {
        indices[j] = getColMap()->getLocalElement(indrowview[j]);
      }
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above if indices are allocated
      TEUCHOS_TEST_FOR_EXCEPTION(
        staticGraph_->indicesAreAllocated(), 
	std::logic_error, "Tpetra::CrsMatrix::getLocalRowCopy: Control flow "
	"should never reach here.  Please report this bug to the Tpetra "
	"developers.");
#endif // HAVE_TPETRA_DEBUG
      numEntries = 0;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalRowCopy(
                                GlobalOrdinal globalRow,
                                const ArrayView<GlobalOrdinal> &indices,
                                const ArrayView<Scalar>        &values,
                                size_t &numEntries) const
  {
    // Only locally owned rows can be queried, otherwise complain
    const std::string tfecfFuncName("getGlobalRowCopy()");
    const LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(lrow == LOT::invalid(), std::runtime_error, ": globalRow does not belong to this node.");
    const RowInfo rowinfo = staticGraph_->getRowInfo(lrow);
    numEntries = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(static_cast<size_t>(indices.size()) < numEntries || static_cast<size_t>(values.size()) < numEntries,
        std::runtime_error, ": size of indices,values must be sufficient to store the specified row.");
    if (staticGraph_->isGloballyIndexed()) {
      ArrayView<const GlobalOrdinal> indrowview = staticGraph_->getGlobalView(rowinfo);
      ArrayView<const Scalar>        valrowview = getView(rowinfo);
      std::copy( indrowview.begin(), indrowview.begin() + numEntries, indices.begin() );
      std::copy( valrowview.begin(), valrowview.begin() + numEntries,  values.begin() );
    }
    else if (staticGraph_->isLocallyIndexed()) {
      ArrayView<const LocalOrdinal> indrowview = staticGraph_->getLocalView(rowinfo);
      ArrayView<const Scalar>       valrowview = getView(rowinfo);
      std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      for (size_t j=0; j < numEntries; ++j) {
        indices[j] = getColMap()->getGlobalElement(indrowview[j]);
      }
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above if indices are allocated
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( staticGraph_->indicesAreAllocated() == true, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
#endif
      numEntries = 0;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalRowView(
                                LocalOrdinal localRow,
                                ArrayView<const LocalOrdinal> &indices,
                                ArrayView<const Scalar>       &values) const
  {
    const std::string tfecfFuncName("getLocalRowView()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed() == true, std::runtime_error, ": local indices cannot be provided.");
    indices = null;
    values  = null;
    if (getRowMap()->isNodeLocalElement(localRow) == true) {
      const RowInfo rowinfo = staticGraph_->getRowInfo(localRow);
      if (rowinfo.numEntries > 0) {
        indices = staticGraph_->getLocalView(rowinfo);
        indices = indices(0,rowinfo.numEntries);
        values  = getView(rowinfo);
        values  = values(0,rowinfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)indices.size() != getNumEntriesInLocalRow(localRow) || indices.size() != values.size(),
        std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalRowView(
                                GlobalOrdinal globalRow,
                                ArrayView<const GlobalOrdinal> &indices,
                                ArrayView<const Scalar>        &values) const
  {
    const std::string tfecfFuncName("getGlobalRowView()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() == true, std::runtime_error, ": global indices cannot be provided.");
    indices = null;
    values  = null;
    const LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    if (lrow != OrdinalTraits<LocalOrdinal>::invalid()) {
      const RowInfo rowinfo = staticGraph_->getRowInfo(lrow);
      if (rowinfo.numEntries > 0) {
        indices = staticGraph_->getGlobalView(rowinfo);
        indices = indices(0,rowinfo.numEntries);
        values  = getView(rowinfo);
        values  = values(0,rowinfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)indices.size() != getNumEntriesInGlobalRow(globalRow) || indices.size() != values.size(),
        std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::scale(const Scalar &alpha)
  {
    const std::string tfecfFuncName("scale()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false, std::runtime_error, " requires that fill is active.");
    // scale all values in the matrix
    // it is easiest to scale all allocated values, instead of scaling only the ones with valid entries
    // however, if there are no valid entries, we can short-circuit
    // furthermore, if the values aren't allocated, we can short-circuit (unallocated values are zero, scaling to zero)
    const size_t     nlrs = staticGraph_->getNodeNumRows(),
                 numAlloc = staticGraph_->getNodeAllocationSize(),
               numEntries = staticGraph_->getNodeNumEntries();
    if (staticGraph_->indicesAreAllocated() == false || numAlloc == 0 || numEntries == 0) {
      // do nothing
    }
    else {
      if (staticGraph_->getProfileType() == StaticProfile) {
        typename ArrayRCP<Scalar>::iterator it;
        for (it = values1D_.begin(); it != values1D_.end(); ++it) {
          (*it) *= alpha;
        }
      }
      else if (staticGraph_->getProfileType() == DynamicProfile) {
        typename ArrayRCP<Scalar>::iterator it;
        for (size_t row=0; row < nlrs; ++row) {
          if (values2D_[row] != null) {
            for (it = values2D_[row].begin(); it != values2D_[row].end(); ++it) {
              (*it) *= alpha;
            }
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setAllToScalar(const Scalar &alpha)
  {
    const std::string tfecfFuncName("setAllToScalar()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false, std::runtime_error, " requires that fill is active.");
    // replace all values in the matrix
    // it is easiest to replace all allocated values, instead of replacing only the ones with valid entries
    // however, if there are no valid entries, we can short-circuit
    // furthermore, if the values aren't allocated, we can short-circuit (no entry have been inserted so far)
    const size_t     nlrs = staticGraph_->getNodeNumRows(),
                 numAlloc = staticGraph_->getNodeAllocationSize(),
               numEntries = staticGraph_->getNodeNumEntries();
    if (staticGraph_->indicesAreAllocated() == false || numAlloc == 0 || numEntries == 0) {
      // do nothing
    }
    else {
      if (staticGraph_->getProfileType() == StaticProfile) {
        std::fill( values1D_.begin(), values1D_.end(), alpha );
      }
      else if (staticGraph_->getProfileType() == DynamicProfile) {
        for (size_t row=0; row < nlrs; ++row) {
          std::fill( values2D_[row].begin(), values2D_[row].end(), alpha );
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &dvec) const
  {
    const std::string tfecfFuncName("getLocalDiagCopy()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isFillComplete() == false, std::runtime_error, " until fillComplete() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(dvec.getMap()->isSameAs(*getRowMap()) == false, std::runtime_error, ": dvec must have the same map as the CrsMatrix.");
    const size_t STINV = OrdinalTraits<size_t>::invalid();
#ifdef HAVE_TPETRA_DEBUG
    size_t numDiagFound = 0;
#endif
    const size_t nlrs = getNodeNumRows();
    ArrayRCP<Scalar> vecView = dvec.get1dViewNonConst();
    RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> > colMap = getColMap();
    for (size_t r=0; r < nlrs; ++r) {
      vecView[r] = ScalarTraits<Scalar>::zero();
      GlobalOrdinal rgid = getRowMap()->getGlobalElement(r);
      if (colMap->isNodeGlobalElement(rgid)) {
        LocalOrdinal rlid = colMap->getLocalElement(rgid);
        RowInfo rowinfo = staticGraph_->getRowInfo(r);
        if (rowinfo.numEntries > 0) {
          const size_t j = staticGraph_->findLocalIndex(rowinfo, rlid);
          ArrayView<const Scalar> view = this->getView(rowinfo);
          if (j != STINV) {
            vecView[r] = view[j];
#ifdef HAVE_TPETRA_DEBUG
            ++numDiagFound;
#endif
          }
        }
      }
    }
    vecView = null;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(numDiagFound != getNodeNumDiags(), std::logic_error, ": logic error. Please contact Tpetra team.");
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::leftScale(
    const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x)
  {
    const std::string tfecfFuncName("leftScale()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!isFillComplete(), std::runtime_error, ": matrix must be fill complete.");
    RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xp = null;
    if(getRangeMap()->isSameAs(*(x.getMap()))){
      // Take from Epetra: If we have a non-trivial exporter, we must
      // import elements that are permuted or are on other processors.
      // (We will use the exporter to perform the import ("reverse
      // mode").)
      if(getCrsGraph()->getExporter() != null){
        RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tempVec
          = rcp(new Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(getRowMap()));
        tempVec->doImport(x, *(getCrsGraph()->getExporter()), INSERT);
        xp = tempVec;
      }
      else{
        xp = rcpFromRef(x);
      }
    }
    else if(getRowMap()->isSameAs(*(x.getMap()))){
      xp = rcpFromRef(x);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::invalid_argument, ": The "
        "input scaling vector x's Map must be the same as either the row Map or "
        "the range Map of the CrsMatrix.");
    }
    ArrayRCP<const Scalar> vectorVals = xp->getData(0);
    ArrayView<Scalar> rowValues = null;
    for(LocalOrdinal i = OrdinalTraits<LocalOrdinal>::zero(); Teuchos::as<size_t>(i) < getNodeNumRows(); ++i){
      const RowInfo rowinfo = staticGraph_->getRowInfo(i);
      rowValues = getViewNonConst(rowinfo);
      Scalar scaleValue = vectorVals[i];
      for(LocalOrdinal j=OrdinalTraits<LocalOrdinal>::zero(); Teuchos::as<size_t>(j)<rowinfo.numEntries; ++j){
        rowValues[j] *= scaleValue;
      }
      rowValues = null;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::rightScale(
    const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x)
  {
    const std::string tfecfFuncName("rightScale()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!isFillComplete(), std::runtime_error, ": matrix must be fill complete.");
    RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xp = null;
    if(getDomainMap()->isSameAs(*(x.getMap()))){
      // Take from Epetra:
      // If we have a non-trivial exporter, we must import elements that are
      // permuted or are on other processors.  (We will use the exporter to
      // perform the import.)
      if(getCrsGraph()->getImporter() != null){
        RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tempVec
          = rcp(new Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(getColMap()));
        tempVec->doImport(x, *(getCrsGraph()->getImporter()), INSERT);
        xp = tempVec;
      }
      else{
        xp = rcpFromRef(x);
      }
    }
    else if(getRowMap()->isSameAs(*(x.getMap()))){
      xp = rcpFromRef(x);
    }
    else{
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::runtime_error, ": The vector x must be the same as either the row map or the range map");
    }

    ArrayRCP<const Scalar> vectorVals = xp->getData(0);
    ArrayView<Scalar> rowValues = null;
    for(
      LocalOrdinal i = OrdinalTraits<LocalOrdinal>::zero();
      Teuchos::as<size_t>(i) < getNodeNumRows();
      ++i)
    {
      const RowInfo rowinfo = staticGraph_->getRowInfo(i);
      rowValues = getViewNonConst(rowinfo);
      ArrayView<const LocalOrdinal> colInices;
      getCrsGraph()->getLocalRowView(i, colInices);
      for(
        LocalOrdinal j = OrdinalTraits<LocalOrdinal>::zero();
        Teuchos::as<size_t>(j) < rowinfo.numEntries;
        ++j
      )
      {
        rowValues[j] *= vectorVals[colInices[j]];
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  typename ScalarTraits<Scalar>::magnitudeType
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getFrobeniusNorm() const
  {
    // TODO: push localFrobNorm() down to the LocalMatOps class
    //
    // check the cache first
    Magnitude frobNorm = frobNorm_;
    if (frobNorm == -ScalarTraits<Magnitude>::one()) {
      Magnitude mySum = ScalarTraits<Magnitude>::zero();
      if (getNodeNumEntries() > 0) {
        if (isStorageOptimized()) {
          // can do this in one pass through A
          typename ArrayRCP<const Scalar>::iterator valit, valend;
          valit = values1D_.begin();
          valend = valit + getNodeNumEntries();
          while (valit != valend) {
            const Scalar val = *valit++;
            mySum += ST::magnitude( ST::conjugate(val) * val );
          }
        }
        else if (getProfileType() == StaticProfile)
        {
          // must hit each row individually
          const size_t numRows = getNodeNumRows();
          for (size_t r=0; r != numRows; ++r)
          {
            typename ArrayRCP<const Scalar>::iterator valit, valend;
            RowInfo rowInfo = myGraph_->getRowInfo(r);
            valit = values1D_.begin() + rowInfo.offset1D;
            valend = valit + rowInfo.numEntries;
            while (valit != valend) {
              const Scalar val = *valit++;
              mySum += ST::magnitude( ST::conjugate(val) * val );
            }
          }
        }
        else if (getProfileType() == DynamicProfile)
        {
          // must hit each row individually
          const size_t numRows = getNodeNumRows();
          for (size_t r=0; r != numRows; ++r)
          {
            typename ArrayRCP<const Scalar>::iterator valit, valend;
            RowInfo rowInfo = myGraph_->getRowInfo(r);
            valit = values2D_[r].begin();
            valend = valit + rowInfo.numEntries;
            while (valit != valend) {
              const Scalar val = *valit++;
              mySum += ST::magnitude( ST::conjugate(val) * val );
            }
          }
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, typeName(*this) << "::getFrobeniusNorm(): Internal logic error. Please contact Tpetra team.");
        }
      }
      Magnitude totalSum;
      Teuchos::reduceAll(*(getComm()), Teuchos::REDUCE_SUM, mySum, outArg(totalSum));
      frobNorm = ScalarTraits<Magnitude>::squareroot(totalSum);
    }
    if (isFillComplete()) {
      // cache the result
      frobNorm_ = frobNorm;
    }
    return frobNorm;
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::globalAssemble()
  {
    using Teuchos::SerialDenseMatrix;
    using std::pair;
    using std::make_pair;
    // Iterator over nonlocals_, the nonlocal data stored by previous
    // calls to insertGlobalValues() for nonowned rows.
    typedef typename std::map<GlobalOrdinal,Array<pair<GlobalOrdinal,Scalar> > >::const_iterator NLITER;
    typedef typename Array<pair<GlobalOrdinal,Scalar> >::const_iterator NLRITER;

    const int numImages = getComm()->getSize();
    const int myImageID = getComm()->getRank();
    const std::string tfecfFuncName("globalAssemble()");
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *getRowMap()->getComm() );
#endif // HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false, std::runtime_error, " requires that fill is active.");

    // Determine (via a global all-reduce) if any processes have
    // global entries to share.
    size_t MyNonlocals = nonlocals_.size(),
           MaxGlobalNonlocals;
    Teuchos::reduceAll<int,size_t>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,
      outArg(MaxGlobalNonlocals));
    if (MaxGlobalNonlocals == 0) return;  // no entries to share

    // FIXME (mfh 14 Dec 2012) The code below reimplements an Export
    // operation.  It would be better just to use an Export.  See
    // Comment #34 in discussion of Bug 5782.

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images: globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int,GlobalOrdinal> > IdsAndRows;
    std::map<GlobalOrdinal,int> NLR2Id;
    SerialDenseMatrix<int,char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all non-local rows
      // we want a list of the rows for which we have data
      Array<GlobalOrdinal> NLRs;
      std::set<GlobalOrdinal> setOfRows;
      for (NLITER iter = nonlocals_.begin(); iter != nonlocals_.end(); ++iter)
      {
        setOfRows.insert(iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize(setOfRows.size());
      std::copy(setOfRows.begin(), setOfRows.end(), NLRs.begin());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds(NLRs.size());
      {
        LookupStatus stat = getRowMap()->getRemoteIndexList(NLRs(),NLRIds());
        char lclerror = ( stat == IDNotPresent ? 1 : 0 );
        char gblerror;
        Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,lclerror,outArg(gblerror));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(gblerror, std::runtime_error, ": non-local entries correspond to invalid rows.");
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages,0);
      typename Array<GlobalOrdinal>::const_iterator nlr;
      typename Array<int>::const_iterator id;
      for (nlr = NLRs.begin(), id = NLRIds.begin();
           nlr != NLRs.end(); ++nlr, ++id)
      {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        IdsAndRows.push_back(make_pair(*id,*nlr));
      }
      for (int j=0; j<numImages; ++j)
      {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      Teuchos::gatherAll(*getComm(),numImages,localNeighbors.getRawPtr(),numImages*numImages,globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images.
      // globalNeighbors(i,j) != 0 means that j sends to i/that i receives from j
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    //////////////////////////////////////////////////////////////////////////////////////

    // loop over all columns to know from which images I can expect to receive something
    for (int j=0; j<numImages; ++j)
    {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    size_t numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<CrsIJV<GlobalOrdinal,Scalar> > IJVSendBuffer;
    Array<size_t> sendSizes(sendIDs.size(), 0);
    size_t numSends = 0;
    for (typename Array<pair<int,GlobalOrdinal> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow)
    {
      int            id = IdAndRow->first;
      GlobalOrdinal row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(sendIDs[numSends] != id, std::logic_error, ": internal logic error. Contact Tpetra team.");
      }
      // copy data for row into contiguous storage
      for (NLRITER jv = nonlocals_[row].begin(); jv != nonlocals_[row].end(); ++jv)
      {
        IJVSendBuffer.push_back( CrsIJV<GlobalOrdinal,Scalar>(row,jv->first,jv->second) );
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size() > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<typename Array<int>::size_type>(numSends) != sendIDs.size(),
        std::logic_error, ": internal logic error. Contact Tpetra team.");

    //
    // don't need this data anymore
    // clear it before we start allocating a bunch of new memory
    nonlocals_.clear();

    //////////////////////////////////////////////////////////////////////////////////////
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    //////////////////////////////////////////////////////////////////////////////////////
    // perform non-blocking sends: send sizes to our recipients
    Array<RCP<Teuchos::CommRequest<int> > > sendRequests;
    for (size_t s=0; s < numSends ; ++s) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back( Teuchos::isend<int,size_t>(*getComm(),rcpFromRef(sendSizes[s]),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<RCP<Teuchos::CommRequest<int> > > recvRequests;
    Array<size_t> recvSizes(numRecvs);
    for (size_t r=0; r < numRecvs; ++r) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      recvRequests.push_back( Teuchos::ireceive(*getComm(),rcp(&recvSizes[r],false),recvIDs[r]) );
    }
    // wait on all
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*getComm(),sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*getComm(),recvRequests());
    }
    Teuchos::barrier(*getComm());
    sendRequests.clear();
    recvRequests.clear();

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW SEND/RECEIVE ALL ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // from the size info, build the ArrayViews into IJVSendBuffer
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > sendBuffers(numSends,null);
    {
      size_t cur = 0;
      for (size_t s=0; s<numSends; ++s) {
        sendBuffers[s] = IJVSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (size_t s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<GlobalOrdinal,Scalar> > tmparcp = arcp(sendBuffers[s].getRawPtr(),0,sendBuffers[s].size(),false);
      sendRequests.push_back( Teuchos::isend<int,CrsIJV<GlobalOrdinal,Scalar> >(*getComm(),tmparcp,sendIDs[s]) );
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    size_t totalRecvSize = std::accumulate(recvSizes.begin(),recvSizes.end(),0);
    Array<CrsIJV<GlobalOrdinal,Scalar> > IJVRecvBuffer(totalRecvSize);
    // from the size info, build the ArrayViews into IJVRecvBuffer
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > recvBuffers(numRecvs,null);
    {
      size_t cur = 0;
      for (size_t r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJVRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (size_t r=0; r < numRecvs ; ++r)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<GlobalOrdinal,Scalar> > tmparcp = arcp(recvBuffers[r].getRawPtr(),0,recvBuffers[r].size(),false);
      recvRequests.push_back( Teuchos::ireceive(*getComm(),tmparcp,recvIDs[r]) );
    }
    // perform waits
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*getComm(),sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*getComm(),recvRequests());
    }
    Teuchos::barrier(*getComm());
    sendRequests.clear();
    recvRequests.clear();


    ////////////////////////////////////////////////////////////////////////////////////
    // NOW PROCESS THE RECEIVED ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // TODO: instead of adding one entry at a time, add one row at a time.
    //       this requires resorting; they arrived sorted by sending node, so that entries could be non-contiguous if we received
    //       multiple entries for a particular row from different processors.
    //       it also requires restoring the data, which may make it not worth the trouble.

    if (this->isStaticGraph ()) {
      for (typename Array<CrsIJV<GlobalOrdinal,Scalar> >::const_iterator ijv = IJVRecvBuffer.begin(); ijv != IJVRecvBuffer.end(); ++ijv) {
        sumIntoGlobalValues (ijv->i, tuple (ijv->j), tuple (ijv->v));
      }
    }
    else { // Dynamic graph; can use insertGlobalValues ()
      for (typename Array<CrsIJV<GlobalOrdinal,Scalar> >::const_iterator ijv = IJVRecvBuffer.begin(); ijv != IJVRecvBuffer.end(); ++ijv) {
        try {
          insertGlobalValues(ijv->i, tuple(ijv->j), tuple(ijv->v));
        }
        catch (std::runtime_error &e) {
          std::ostringstream outmsg;
          outmsg << e.what() << std::endl
                 << "caught in globalAssemble() in " << __FILE__ << ":" << __LINE__
                 << std::endl ;
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, outmsg.str());
        }
      }
    }

    // WHEW! THAT WAS TIRING!
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node,
            class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  resumeFill (const RCP<ParameterList> &params)
  {
    const std::string tfecfFuncName("resumeFill");
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *getRowMap()->getComm() );
#endif // HAVE_TPETRA_DEBUG

    // We use -1 to indicate that the Frobenius norm need to be recomputed.
    frobNorm_ = -Teuchos::ScalarTraits<Magnitude>::one();
    if (! isStaticGraph()) { // Don't resume fill of a nonowned graph.
      myGraph_->resumeFill (params);
    }
    clearGlobalConstants();
    lclMatrix_ = null;
    lclMatOps_ = null;
    fillComplete_ = false;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive() || isFillComplete(), std::logic_error,
      "::resumeFill(): At end of method, either fill is not active or fill is "
      "complete.  This violates stated post-conditions.  Please report this bug "
      "to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::computeGlobalConstants() {
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::clearGlobalConstants() {
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(const RCP<ParameterList> &params) {
    fillComplete(getRowMap(),getRowMap(),params);
  }



  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  fillComplete (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap,
                const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap,
                const RCP<ParameterList> &params)
  {
    const std::string tfecfFuncName("fillComplete()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive() || isFillComplete(),
      std::runtime_error, ": Matrix fill state must be active (isFillActive() "
      "must be true) before calling fillComplete().");
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *getRowMap()->getComm() );
#endif // HAVE_TPETRA_DEBUG

    // allocate if unallocated
    if (! getCrsGraph()->indicesAreAllocated()) {
      // allocate global, in case we do not have a column map
      allocateValues( GlobalIndices, GraphNotYetAllocated );
    }
    // Global assemble, if we need to (we certainly don't need to if
    // there's only one process).  This call only costs a single
    // all-reduce if we don't need global assembly.
    if (getComm()->getSize() > 1) {
      // this calls insertGlobalValues(), one entry at a time.
      globalAssemble ();
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(nonlocals_.size() > 0,
        std::runtime_error, ": cannot have nonlocal entries on a serial run.  "
        "An invalid entry (i.e., with row index not in the row Map) must have "
        "been submitted to the CrsMatrix.");
    }

    if (isStaticGraph ()) {
      const bool domainMapsMatch = staticGraph_->getDomainMap() == domainMap;
      const bool rangeMapsMatch = staticGraph_->getRangeMap() == rangeMap;
      // FIXME (mfh 19 Mar 2012) Why can't we allow the Maps to be
      // different objects, but semantically the same (in the sense of
      // isSameAs())?
      // (cgb 24 May 2012) We can/should. We can fix now or wait for a user to complain.
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! domainMapsMatch, std::runtime_error,
        ": The CrsMatrix's domain Map does not match the graph's domain Map.  "
        "The graph cannot be changed because it was given to the CrsMatrix "
        "constructor as const.");

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! rangeMapsMatch, std::runtime_error,
        ": The CrsMatrix's range Map does not match the graph's range Map.  "
        "The graph cannot be changed because it was given to the CrsMatrix "
        "constructor as const.");
    }
    else {
      // Set the graph's domain and range Maps.
      // This may clear the Import/Export objects.
      myGraph_->setDomainRangeMaps(domainMap, rangeMap);
      // Make the graph's column Map, if necessary.
      if (! myGraph_->hasColMap()) {
        myGraph_->makeColMap();
      }
      // make indices local
      if (myGraph_->isGloballyIndexed()) {
        myGraph_->makeIndicesLocal();
      }
      if (! myGraph_->isSorted()) sortEntries();
      if (! myGraph_->isMerged()) mergeRedundantEntries();
      myGraph_->makeImportExport(); // Make Import and Export objects
      myGraph_->computeGlobalConstants();
      myGraph_->fillComplete_ = true;
      myGraph_->checkInternalState();
    }
    computeGlobalConstants();
    // fill local objects; will fill and finalize local graph if appropriate
    if (myGraph_ != null) {
      // The matrix owns the graph, so fill the local graph at the
      // same time as the local matrix.
      fillLocalGraphAndMatrix(params);
    }
    else {
      // The matrix does _not_ own the graph, and the graph's
      // structure is already fixed, so just fill the local matrix.
      fillLocalMatrix(params);
    }
    //
    // Set up the local sparse kernels.
    //
    lclMatOps_ = rcp (new sparse_ops_type (getNode ()));
    // This is where we take the local graph and matrix, and turn them
    // into (possibly optimized) sparse kernels.
    lclMatOps_->setGraphAndMatrix (staticGraph_->getLocalGraph (), lclMatrix_);

    // Once we've initialized the sparse kernels, we're done with the
    // local objects.  We may now release them and their memory, since
    // they will persist in the local sparse ops if necessary.  We
    // keep the local graph if the parameters tell us to do so.
    lclMatrix_ = null;
    if (myGraph_ != null) {
      bool preserveLocalGraph = false;
      if (params != null) {
        preserveLocalGraph = params->get ("Preserve Local Graph", false);
      }
      if (! preserveLocalGraph) {
        myGraph_->lclGraph_ = null;
      }
    }
    // Now we're fill complete!
    fillComplete_ = true;

    // Sanity checks at the end.
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isFillActive(), std::logic_error,
      ": We're at the end of fillComplete(), but isFillActive() is true.  "
      "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! isFillComplete(), std::logic_error,
      ": We're at the end of fillComplete(), but isFillActive() is true.  "
      "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sortEntries()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        typeName(*this) << "::sortEntries(): cannot sort with static graph.");
    if (myGraph_->isSorted() == false) {
      for (size_t row=0; row < getNodeNumRows(); ++row) {
        RowInfo rowInfo = myGraph_->getRowInfo(row);
        myGraph_->template sortRowIndicesAndValues<Scalar>(rowInfo,this->getViewNonConst(rowInfo));
      }
      // we just sorted every row
      myGraph_->setSorted(true);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::mergeRedundantEntries()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        typeName(*this) << "::mergeRedundantEntries(): cannot merge with static graph.");
    if (myGraph_->isMerged() == false) {
      for (size_t row=0; row < getNodeNumRows(); ++row) {
        RowInfo rowInfo = myGraph_->getRowInfo(row);
        myGraph_->template mergeRowIndicesAndValues<typename ArrayRCP<Scalar>::iterator>(rowInfo,this->getViewNonConst(rowInfo).begin(), std::plus<Scalar>());
      }
      // we just merged every row
      myGraph_->setMerged(true);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::apply(
                                        const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                                        MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                                        Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
    TEUCHOS_TEST_FOR_EXCEPTION( isFillComplete() == false, std::runtime_error,
        typeName(*this) << "::apply(): cannot call apply() until fillComplete() has been called.");
    sameScalarMultiplyOp_->apply(X,Y,mode,alpha,beta);
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  gaussSeidel (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
               const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
               const Scalar& dampingFactor,
               const ESweepDirection direction,
               const int numSweeps) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( isFillComplete() == false, std::runtime_error,
      typeName(*this) << "::gaussSeidel: cannot call this method until "
      "fillComplete() has been called.");
    sameScalarMultiplyOp_->gaussSeidel (B, X, D, dampingFactor, direction, numSweeps);
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  gaussSeidelCopy (MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                   const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                   const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
                   const Scalar& dampingFactor,
                   const ESweepDirection direction,
                   const int numSweeps) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( isFillComplete() == false, std::runtime_error,
      typeName(*this) << "::gaussSeidelCopy: cannot call this method until "
      "fillComplete() has been called.");
    sameScalarMultiplyOp_->gaussSeidelCopy (X, B, D, dampingFactor,
                                            direction, numSweeps);
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class DomainScalar, class RangeScalar>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::localMultiply(
                                        const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                                              MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                                              Teuchos::ETransp mode, RangeScalar alpha, RangeScalar beta) const
  {
    using Teuchos::NO_TRANS;
    const std::string tfecfFuncName("localMultiply()");
    typedef ScalarTraits<RangeScalar> RST;
    const Kokkos::MultiVector<DomainScalar,Node> *lclX = &X.getLocalMV();
    Kokkos::MultiVector<RangeScalar,Node>        *lclY = &Y.getLocalMVNonConst();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(mode == NO_TRANS && X.getMap() != getColMap() && *X.getMap() != *getColMap(), std::runtime_error, " X is not distributed according to the appropriate map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(mode != NO_TRANS && X.getMap() != getRowMap() && *X.getMap() != *getRowMap(), std::runtime_error, " X is not distributed according to the appropriate map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(mode == NO_TRANS && Y.getMap() != getRowMap() && *Y.getMap() != *getRowMap(), std::runtime_error, " Y is not distributed according to the appropriate map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(mode != NO_TRANS && Y.getMap() != getColMap() && *Y.getMap() != *getColMap(), std::runtime_error, " Y is not distributed according to the appropriate map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!isFillComplete(),                                              std::runtime_error, " until fillComplete() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(X.getNumVectors() != Y.getNumVectors(),                         std::runtime_error, ": X and Y must have the same number of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(X.isConstantStride() == false || Y.isConstantStride() == false, std::runtime_error, ": X and Y must be constant stride.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(lclX==lclY,                                                     std::runtime_error, ": X and Y cannot share data.");
#endif
    //
    // Call the matvec
    if (beta == RST::zero()) {
      // Y = alpha*op(M)*X with overwrite semantics
      lclMatOps_->template multiply<DomainScalar,RangeScalar>(mode, alpha, *lclX, *lclY);
    }
    else {
      // Y = alpha*op(M) + beta*Y
      lclMatOps_->template multiply<DomainScalar,RangeScalar>(mode, alpha, *lclX, beta, *lclY);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class DomainScalar, class RangeScalar>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  localGaussSeidel (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                    MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
                    const RangeScalar& dampingFactor,
                    const Kokkos::ESweepDirection direction) const
  {
    Kokkos::MultiVector<DomainScalar,Node>& x = X.getLocalMVNonConst ();
    const Kokkos::MultiVector<RangeScalar,Node>& b = B.getLocalMV ();
    const Kokkos::MultiVector<RangeScalar,Node>& d = D.getLocalMV ();

    lclMatOps_->template gaussSeidel<DomainScalar, RangeScalar> (b, x, d, dampingFactor, direction);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class DomainScalar, class RangeScalar>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::localSolve(
                                    const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>  &Y,
                                          MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                                          Teuchos::ETransp mode) const
  {
    using Teuchos::NO_TRANS;
    const std::string tfecfFuncName("localSolve()");
    const Kokkos::MultiVector<RangeScalar,Node> *lclY = &Y.getLocalMV();
    Kokkos::MultiVector<DomainScalar,Node>      *lclX = &X.getLocalMVNonConst();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!isFillComplete(),                                              std::runtime_error, " until fillComplete() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(X.getNumVectors() != Y.getNumVectors(),                         std::runtime_error, ": X and Y must have the same number of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(X.isConstantStride() == false || Y.isConstantStride() == false, std::runtime_error, ": X and Y must be constant stride.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isUpperTriangular() == false && isLowerTriangular() == false,   std::runtime_error, ": can only solve() triangular matrices.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS,      std::logic_error, " does not currently support transposed solve for complex scalar types.");
#endif
    //
    // Call the solve
    if (mode == Teuchos::NO_TRANS) {
      lclMatOps_->template solve<DomainScalar,RangeScalar>(Teuchos::NO_TRANS, *lclY, *lclX);
    }
    else {
      lclMatOps_->template solve<DomainScalar,RangeScalar>(Teuchos::CONJ_TRANS, *lclY, *lclX);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class T>
  RCP<CrsMatrix<T,LocalOrdinal,GlobalOrdinal,Node> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::convert() const
  {
    const std::string tfecfFuncName("convert()");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isFillComplete() == false, std::runtime_error, ": fill must be complete.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(getCrsGraph()->getLocalGraph() == null, std::runtime_error,
        ": local graph data was deleted during fillComplete().\n"
        "To allow convert(), set the following to fillComplete():\n"
        "   \"Preserve Local Graph\" == true ");
    RCP<CrsMatrix<T,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > newmat;
    newmat = rcp(new CrsMatrix<T,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>(getCrsGraph()));
    const Map<LocalOrdinal,GlobalOrdinal,Node> &rowMap = *getRowMap();
    Array<T> newvals;
    for (LocalOrdinal li=rowMap.getMinLocalIndex(); li <= rowMap.getMaxLocalIndex(); ++li)
    {
      ArrayView<const LocalOrdinal> rowinds;
      ArrayView<const Scalar>       rowvals;
      this->getLocalRowView(li,rowinds,rowvals);
      if (rowvals.size() > 0) {
        newvals.resize(rowvals.size());
        std::transform( rowvals.begin(), rowvals.end(), newvals.begin(), Teuchos::asFunc<T>() );
        newmat->replaceLocalValues(li, rowinds, newvals());
      }
    }
    newmat->fillComplete(this->getDomainMap(), this->getRangeMap());
    return newmat;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::checkInternalState() const
  {
#ifdef HAVE_TPETRA_DEBUG
    const std::string tfecfFuncName("checkInternalState()");
    const std::string err(": Likely internal logic error. Please contact Tpetra team.");
    RCP<Node> node = getNode();
    // check the internal state of this data structure
    // this is called by numerous state-changing methods, in a debug build, to ensure that the object
    // always remains in a valid state

    // we must have a static graph
    //
    // a dynamic graph, depending on which constructor was used.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( staticGraph_ == null,                                             std::logic_error, err);
    // i only ever have a local matrix for a small moment in time
    TEUCHOS_TEST_FOR_EXCEPTION( lclMatrix_ != null,                                                          std::logic_error, err );
    // if active, i have no local sparse ops
    TEUCHOS_TEST_FOR_EXCEPTION( isFillActive() && lclMatOps_ != null,                                        std::logic_error, err );
    // if filled, i have a local sparse ops
    TEUCHOS_TEST_FOR_EXCEPTION( isFillComplete() && lclMatOps_ == null,                                      std::logic_error, err );
    // myGraph == null means that the matrix has a static graph.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( myGraph_ != null && myGraph_ != staticGraph_,                     std::logic_error, err);
    // if matrix is fill complete, then graph must be fill complete
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( fillComplete_ == true && staticGraph_->isFillComplete() == false, std::logic_error, err);
    // if matrix is storage optimized, it should have a 1D allocation
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isStorageOptimized() == true && values2D_ != null,                std::logic_error, err);
    // if matrix/graph are static profile, then 2D allocation should not be present
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getProfileType() == StaticProfile  && values2D_ != null,          std::logic_error, err);
    // if matrix/graph are dynamic profile, then 1D allocation should not be present
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getProfileType() == DynamicProfile && values1D_ != null,          std::logic_error, err);
    // if values are allocated and they are non-zero in number, then one of the allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( staticGraph_->indicesAreAllocated()
                        && staticGraph_->getNodeAllocationSize() > 0 && staticGraph_->getNodeNumRows() > 0
                        && values2D_ == null && values1D_ == null,                                   std::logic_error, err);
    // we cannot have both a 1D and 2D allocation
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( values1D_ != null && values2D_ != null,                           std::logic_error, err);
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::description() const {
    std::ostringstream oss;
    oss << DistObject<char, LocalOrdinal,GlobalOrdinal,Node>::description();
    if (isFillComplete()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::as;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;

    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) {
      vl = VERB_LOW;
    }
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank();
    const int numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t> (width, as<size_t> (11)) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl;
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << getGlobalNumDiags() << std::endl;
        out << "Global max number of row entries = " << getGlobalMaxNumRowEntries() << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        getRowMap()->describe(out,vl);
        //
        if (getColMap() != null) {
          if (getColMap() == getRowMap()) {
            if (myImageID == 0) out << "\nColumn map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nColumn map: " << std::endl;
            getColMap()->describe(out,vl);
          }
        }
        if (getDomainMap() != null) {
          if (getDomainMap() == getRowMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else if (getDomainMap() == getColMap()) {
            if (myImageID == 0) out << "\nDomain map is col map.";
          }
          else {
            if (myImageID == 0) out << "\nDomain map: " << std::endl;
            getDomainMap()->describe(out,vl);
          }
        }
        if (getRangeMap() != null) {
          if (getRangeMap() == getDomainMap()) {
            if (myImageID == 0) out << "\nRange map is domain map." << std::endl;
          }
          else if (getRangeMap() == getRowMap()) {
            if (myImageID == 0) out << "\nRange map is row map." << std::endl;
          }
          else {
            if (myImageID == 0) out << "\nRange map: " << std::endl;
            getRangeMap()->describe(out,vl);
          }
        }
        if (myImageID == 0) out << std::endl;
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl;
            if (staticGraph_->indicesAreAllocated() == false) {
              out << "Node not allocated" << std::endl;
            }
            else {
              out << "Node number of allocated entries = " << staticGraph_->getNodeAllocationSize() << std::endl;
            }
            out << "Node number of entries = " << getNodeNumEntries() << std::endl;
            if (isFillComplete()) {
              out << "Node number of diagonals = " << getNodeNumDiags() << std::endl;
            }
            out << "Node max number of entries = " << getNodeMaxNumRowEntries() << std::endl;
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row"
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << std::setw(width) << "(Index,Value)";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              const size_t nE = getNumEntriesInLocalRow(r);
              GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
              out << std::setw(width) << myImageID
                  << std::setw(width) << gid
                  << std::setw(width) << nE;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  ArrayView<const GlobalOrdinal> rowinds;
                  ArrayView<const Scalar> rowvals;
                  getGlobalRowView(gid,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << rowinds[j]
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const LocalOrdinal> rowinds;
                  ArrayView<const Scalar> rowvals;
                  getLocalRowView(r,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << getColMap()->getGlobalElement(rowinds[j])
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  bool
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  checkSizes (const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source)
  {
    // It's not clear what kind of compatibility checks on sizes can be performed here.
    // Epetra_CrsGraph doesn't check any sizes for compatibility.

    // right now, we'll only support import/exporting between CrsMatrix<Scalar>
    // if the source dist object isn't CrsMatrix or some offspring, flag this operation as incompatible.
    try {
      typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> this_type;
      const this_type& A = dynamic_cast<const this_type&> (source);
      (void) A;
    }
    catch (...) {
      // If the input isn't even a CrsMatrix, then certainly the sizes
      // don't match.
      return false;
    }
    return true;
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  copyAndPermute (const DistObject<char, LocalOrdinal,GlobalOrdinal,Node> & source,
                  size_t numSameIDs,
                  const ArrayView<const LocalOrdinal> &permuteToLIDs,
                  const ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    using Teuchos::as;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Node NT;

    // Method name string for TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC.
    const std::string tfecfFuncName("copyAndPermute()");

    // This dynamic cast should succeed, because we've already tested
    // it in checkSizes().
    typedef CrsMatrix<Scalar, LO, GO, NT, LocalMatOps> this_type;
    const this_type& sourceMatrix = dynamic_cast<const this_type&> (source);

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(permuteToLIDs.size() != permuteFromLIDs.size(),
      std::invalid_argument, "permuteToLIDs.size() = " << permuteToLIDs.size()
      << "!= permuteFromLIDs.size() = " << permuteFromLIDs.size() << ".");

    const bool sourceIsLocallyIndexed = sourceMatrix.isLocallyIndexed();
    //
    // Copy the first numSame row from source to target (this matrix).
    // This involves copying rows corresponding to LIDs [0, numSame-1].
    //
    Array<GO> rowInds;
    Array<Scalar> rowVals;
    LO sourceLID = 0;
    for (size_t i = 0; i < numSameIDs; ++i, ++sourceLID) {
      // Global ID for the current row index in the source matrix.
      // The first numSameIDs GIDs in the two input lists are the
      // same, so sourceGID == targetGID in this case.
      const GO sourceGID = sourceMatrix.getMap()->getGlobalElement (sourceLID);
      const GO targetGID = sourceGID;

      // Input views for the combineGlobalValues() call below.
      ArrayView<const GO> rowIndsConstView;
      ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {
        const size_t rowLength = sourceMatrix.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > as<size_t> (rowInds.size())) {
          rowInds.resize (rowLength);
          rowVals.resize (rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        ArrayView<GO> rowIndsView = rowInds.view (0, rowLength);
        ArrayView<Scalar> rowValsView = rowVals.view (0, rowLength);

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        sourceMatrix.getGlobalRowCopy (sourceGID, rowIndsView, rowValsView, checkRowLength);

#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(rowLength != checkRowLength,
          std::logic_error, ": For global row index " << sourceGID << ", the source"
          " matrix's getNumEntriesInGlobalRow() method returns a row length of "
          << rowLength << ", but the getGlobalRowCopy() method reports that "
          "the row length is " << checkRowLength << ".  Please report this bug "
          "to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

        rowIndsConstView = rowIndsView.view (0, rowLength);
        rowValsConstView = rowValsView.view (0, rowLength);
      }
      else { // source matrix is globally indexed.
        sourceMatrix.getGlobalRowView (sourceGID, rowIndsConstView, rowValsConstView);
      }

      // Combine the data into the target matrix.
      if (isStaticGraph()) {
        // Applying a permutation to a matrix with a static graph
        // means REPLACE-ing entries.
        combineGlobalValues (targetGID, rowIndsConstView, rowValsConstView, REPLACE);
      }
      else {
        // Applying a permutation to a matrix with a dynamic graph
        // means INSERT-ing entries.  This has the same effect as
        // ADD, if the target graph already has an entry there.
        combineGlobalValues (targetGID, rowIndsConstView, rowValsConstView, INSERT);
      }
    } // For each of the consecutive source and target IDs that are the same

    //
    // Permute the remaining rows.
    //
    for (size_t p = 0; p < (size_t)permuteToLIDs.size(); ++p) {
      const GO sourceGID = sourceMatrix.getMap()->getGlobalElement (permuteFromLIDs[p]);
      const GO targetGID = this->getMap()->getGlobalElement (permuteToLIDs[p]);

      // Input views for the combineGlobalValues() call below.
      ArrayView<const GO> rowIndsConstView;
      ArrayView<const Scalar> rowValsConstView;

      if (sourceIsLocallyIndexed) {
        const size_t rowLength = sourceMatrix.getNumEntriesInGlobalRow (sourceGID);
        if (rowLength > as<size_t> (rowInds.size())) {
          rowInds.resize (rowLength);
          rowVals.resize (rowLength);
        }
        // Resizing invalidates an Array's views, so we must make new
        // ones, even if rowLength hasn't changed.
        ArrayView<GO> rowIndsView = rowInds.view (0, rowLength);
        ArrayView<Scalar> rowValsView = rowVals.view (0, rowLength);

        // The source matrix is locally indexed, so we have to get a
        // copy.  Really it's the GIDs that have to be copied (because
        // they have to be converted from LIDs).
        size_t checkRowLength = 0;
        sourceMatrix.getGlobalRowCopy (sourceGID, rowIndsView, rowValsView, checkRowLength);

#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(rowLength != checkRowLength,
          std::logic_error, ": For the source matrix's global row index "
          << sourceGID << ", the source matrix's getNumEntriesInGlobalRow() method "
          "returns a row length of " << rowLength << ", but the "
          "getGlobalRowCopy() method reports that the row length is "
          << checkRowLength << ".  Please report this bug to the Tpetra "
          "developers.");
#endif // HAVE_TPETRA_DEBUG

        rowIndsConstView = rowIndsView.view (0, rowLength);
        rowValsConstView = rowValsView.view (0, rowLength);
      }
      else {
        sourceMatrix.getGlobalRowView (sourceGID, rowIndsConstView, rowValsConstView);
      }

      // Combine the data into the target matrix.
      if (isStaticGraph()) {
        combineGlobalValues (targetGID, rowIndsConstView, rowValsConstView, REPLACE);
      }
      else {
        combineGlobalValues (targetGID, rowIndsConstView, rowValsConstView, INSERT);
      }
    } // For each ID to permute
  }


  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  packAndPrepare (const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                  const ArrayView<const LocalOrdinal>& exportLIDs,
                  Array<char>& exports,
                  const ArrayView<size_t>& numPacketsPerLID,
                  size_t& constantNumPackets,
                  Distributor &distor)
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const LO>::size_type size_type;
    const std::string tfecfFuncName ("packAndPrepare()");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(exportLIDs.size() != numPacketsPerLID.size(),
      std::invalid_argument, "exportLIDs.size() = " << exportLIDs.size()
      << "!= numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

    // This dynamic cast should succeed.  The packAndPrepare() method
    // is invoked by DistObject::doTransfer(), which calls
    // checkSizes() and copyAndPermute() first.  At least one of the
    // latter two methods must have successfully done the cast in
    // order for doTransfer() to get to this point.
    typedef CrsMatrix<Scalar, LO, GO, Node, LocalMatOps> this_type;
    const this_type& src_mat = dynamic_cast<const this_type&> (source);
    const bool src_is_locally_indexed = src_mat.isLocallyIndexed();
    constantNumPackets = 0;

    // first, set the contents of numPacketsPerLID, and accumulate a total-num-packets:
    // grab the max row size, while we're at it. may need it below.
    // Subtle: numPacketsPerLID is for byte-packets, so it needs to be multiplied
    const size_t SizeOfOrdValPair = sizeof(GO) + sizeof(Scalar);
    size_t totalNumEntries = 0;
    size_t maxExpRowLength = 0;
    for (size_type i = 0; i < exportLIDs.size(); ++i) {
      const GO expGID = src_mat.getMap()->getGlobalElement(exportLIDs[i]);
      const size_t row_length = src_mat.getNumEntriesInGlobalRow(expGID);
      numPacketsPerLID[i] = row_length * SizeOfOrdValPair;
      totalNumEntries += row_length;
      maxExpRowLength = (row_length > maxExpRowLength) ? row_length : maxExpRowLength;
    }

    // Pack export data by interleaving rows' indices and values in
    // the following way:
    //
    // [inds_row0 vals_row0 inds_row1 vals_row1 ... inds_rowN vals_rowN]
    if (totalNumEntries > 0) {
      // exports is an array of char (bytes). it needs room for all of the indices and values
      const size_t totalNumBytes = totalNumEntries * SizeOfOrdValPair;
      exports.resize(totalNumBytes);

      ArrayView<char> avIndsC, avValsC;
      ArrayView<GO>   avInds;
      ArrayView<Scalar>   avVals;

      // now loop again and pack rows of indices into exports:
      // if global indices exist in the source, then we can use view semantics
      // otherwise, we are forced to use copy semantics (for the indices; for simplicity, we'll use them for values as well)
      size_t curOffsetInBytes = 0;
      if (src_is_locally_indexed) {
        Array<GO> row_inds (maxExpRowLength);
        Array<Scalar> row_vals (maxExpRowLength);
        for (size_type i = 0; i < exportLIDs.size(); ++i) {
          // Get a copy of the current row's data.  We get a copy and
          // not a view, because the indices are stored as local
          // indices, not as global indices.
          //
          // TODO (mfh 14 Mar 2012) It might save some copying to add
          // a method that gets a view of the values but a copy of the
          // global indices.
          const GO GID = src_mat.getMap()->getGlobalElement(exportLIDs[i]);
          size_t rowSize;
          src_mat.getGlobalRowCopy(GID, row_inds(), row_vals(), rowSize);

          // Get views of the spots in the exports array in which to
          // put the indices resp. values.  The type cast makes the
          // views look like GO resp. Scalar, when the array they are
          // viewing is really an array of char.
          //
          // TODO (mfh 14 Mar 2012): Why do we need the reinterpret
          // cast?  Why can't we just store pairs?  Is it because
          // there are no Comm functions for sending and receiving
          // pairs?  How hard can that be to implement?
          avIndsC = exports(curOffsetInBytes, rowSize*sizeof(GO));
          avValsC = exports(curOffsetInBytes+rowSize*sizeof(GO), rowSize*sizeof(Scalar));
          avInds = av_reinterpret_cast<GO> (avIndsC);
          avVals = av_reinterpret_cast<Scalar> (avValsC);
          // Copy the source matrix's row data into the views of the
          // exports array for indices resp. values.
          std::copy (row_inds.begin(), row_inds.begin()+rowSize, avInds.begin());
          std::copy (row_vals.begin(), row_vals.begin()+rowSize, avVals.begin());
          // Keep track of how many bytes we packed.
          curOffsetInBytes += SizeOfOrdValPair * rowSize;
        }
      }
      else { // the source matrix's indices are stored as GIDs, not LIDs.
        ArrayView<const GO> row_inds;
        ArrayView<const Scalar> row_vals;
        for (size_type i = 0; i < exportLIDs.size(); ++i) {
          // Get a view of the current row's data.  We don't need to
          // get a copy, since the source matrix's indices are stored
          // as GIDs.
          const GO GID = src_mat.getMap()->getGlobalElement(exportLIDs[i]);
          src_mat.getGlobalRowView(GID, row_inds, row_vals);
          const size_t rowSize = static_cast<size_t> (row_inds.size());
          // Get views of the spots in the exports array in which to
          // put the indices resp. values.  See notes and FIXME above.
          avIndsC = exports(curOffsetInBytes, rowSize*sizeof(GO));
          avValsC = exports(curOffsetInBytes+rowSize*sizeof(GO), rowSize*sizeof(Scalar));
          avInds = av_reinterpret_cast<GO> (avIndsC);
          avVals = av_reinterpret_cast<Scalar> (avValsC);
          // Copy the source matrix's row data into the views of the
          // exports array for indices resp. values.
          std::copy (row_inds.begin(), row_inds.end(), avInds.begin());
          std::copy (row_vals.begin(), row_vals.end(), avVals.begin());
          curOffsetInBytes += SizeOfOrdValPair * rowSize;
        }
      }
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(curOffsetInBytes != totalNumBytes,
        std::logic_error, "packAndPrepare: At end of method, the final offset "
        "bytes curOffsetInBytes=" << curOffsetInBytes << " != total number of "
        "bytes totalNumBytes=" << totalNumBytes << ".  Please report this bug "
        "to the Tpetra developers.");
#endif //  HAVE_TPETRA_DEBUG
    }
  }

  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  combineGlobalValues (const GlobalOrdinal globalRowIndex,
                       const ArrayView<const GlobalOrdinal> columnIndices,
                       const ArrayView<const Scalar> values,
                       const Tpetra::CombineMode combineMode)
  {
    if (isStaticGraph()) {
      // INSERT doesn't make sense for a static graph, since you
      // aren't allowed to change the structure of the graph.
      // However, all the other combine modes work.
      if (combineMode == ADD) {
        sumIntoGlobalValues (globalRowIndex, columnIndices(), values());
      }
      else if (combineMode == REPLACE) {
        replaceGlobalValues (globalRowIndex, columnIndices(), values());
      }
      else if (combineMode == ABSMAX) {
        using Details::AbsMax;
        AbsMax<Scalar> f;
        this->template transformGlobalValues<AbsMax<Scalar> > (globalRowIndex,
                                                               columnIndices(),
                                                               values(), f);
      }
      else if (combineMode == INSERT) {
        TEUCHOS_TEST_FOR_EXCEPTION(isStaticGraph() && combineMode == INSERT,
          std::invalid_argument, "combineGlobalValues: INSERT combine mode "
          "is not allowed if the matrix has a static graph (i.e., was "
          "constructed with the CrsMatrix constructor that takes a const "
          "CrsGraph pointer).");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "combineGlobalValues: Should never get here!  Please report this bug"
          "to the Tpetra developers.");
      }
    }
    else { // The matrix has a dynamic graph.
      if (combineMode == ADD || combineMode == INSERT) {
        // For a dynamic graph, all incoming column indices are
        // inserted into the target graph.  Duplicate indices will
        // have their values summed.  In this context, ADD and INSERT
        // are equivalent.  We need to call insertGlobalValues()
        // anyway if the column indices don't yet exist in this row,
        // so we just call insertGlobalValues() for both cases.
        insertGlobalValues (globalRowIndex, columnIndices(), values());
      }
      // FIXME (mfh 14 Mar 2012):
      //
      // Implementing ABSMAX or REPLACE for a dynamic graph would
      // require modifying assembly to attach a possibly different
      // combine mode to each inserted (i, j, A_ij) entry.  For
      // example, consider two different Export operations to the same
      // target CrsMatrix, the first with ABSMAX combine mode and the
      // second with REPLACE.  This isn't a common use case, so we
      // won't mess with it for now.
      else if (combineMode == ABSMAX) {
        TEUCHOS_TEST_FOR_EXCEPTION(! isStaticGraph() && combineMode == ABSMAX,
          std::logic_error, "combineGlobalValues: ABSMAX combine mode when "
          "the matrix has a dynamic graph is not yet implemented.");
      }
      else if (combineMode == REPLACE) {
        TEUCHOS_TEST_FOR_EXCEPTION(! isStaticGraph() && combineMode == REPLACE,
          std::logic_error, "combineGlobalValues: REPLACE combine mode when "
          "the matrix has a dynamic graph is not yet implemented.");
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "combineGlobalValues: Should never get here!  Please report this bug"
          "to the Tpetra developers.");
      }
    }
  }



  template<class Scalar,
           class LocalOrdinal,
           class GlobalOrdinal,
           class Node,
           class LocalMatOps>
  void
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  unpackAndCombine (const ArrayView<const LocalOrdinal> &importLIDs,
                    const ArrayView<const char> &imports,
                    const ArrayView<size_t> &numPacketsPerLID,
                    size_t constantNumPackets,
                    Distributor & /* distor */,
                    CombineMode combineMode)
  {
    const std::string tfecfFuncName("unpackAndCombine()");
    const CombineMode validModes[4] = {ADD, REPLACE, ABSMAX, INSERT};
    const char* validModeNames[4] = {"ADD", "REPLACE", "ABSMAX", "INSERT"};
    const int numValidModes = 4;

    if (std::find (validModes, validModes+numValidModes, combineMode) ==
        validModes+numValidModes) {
      std::ostringstream os;
      os << "unpackAndCombine: Invalid combine mode.  Valid modes are {";
      for (int k = 0; k < numValidModes; ++k) {
        os << validModeNames[k];
        if (k < numValidModes - 1) {
          os << ", ";
        }
      }
      os << "}.";
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::invalid_argument, os.str());
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(importLIDs.size() != numPacketsPerLID.size(),
      std::invalid_argument, "importLIDs.size() = " << importLIDs.size()
      << "!= numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

    // FIXME (mfh 05 Dec 2012) Here are all the assumptions encoded in
    // the following line of code:
    //
    // 1. The data (index,value) for each element are packed tightly,
    //    with no extra space in between.
    //
    // 2. sizeof(Scalar) says how much data were used to represent a
    //    Scalar in its packed form.
    //
    // 3. All processes and all instances of Scalar use the same
    //    amount of data to represent a Scalar.  (GlobalOrdinal is
    //    typically a built-in integer type, so this is generally true
    //    for GlobalOrdinal.)
    //
    const size_t SizeOfOrdValPair = sizeof(GlobalOrdinal)+sizeof(Scalar);
    const size_t totalNumBytes = imports.size(); // * sizeof(char), which is one.
    const size_t totalNumEntries = totalNumBytes / SizeOfOrdValPair;

    if (totalNumEntries > 0) {
      // data packed as follows:
      // [inds_row0 vals_row0 inds_row1 vals_row1 ...]
      ArrayView<const char> avIndsC, avValsC;
      ArrayView<const GlobalOrdinal> avInds;
      ArrayView<const Scalar>        avVals;

      size_t curOffsetInBytes = 0;
      typedef typename ArrayView<const LocalOrdinal>::size_type size_type;
      for (size_type i = 0; i < importLIDs.size(); ++i) {
        // get row info
        const LocalOrdinal LID = importLIDs[i];
        const GlobalOrdinal myGID = this->getMap()->getGlobalElement(LID);
        const size_t rowSize = numPacketsPerLID[i] / SizeOfOrdValPair;
        // Needs to be in here in case of zero length rows.  If not,
        // the lines following the if statement error out if the row
        // length is zero. KLN 13/06/2011
        //
        // mfh 05 Dec 2012: The problem to which Kurtis refers in the
        // above comment may no longer be an issue, since
        // ArrayView::view() (which implements ArrayView::operator())
        // now allows views of length zero.
        if (rowSize == 0) {
          continue;
        }
        // Get views of the import (incoming data) buffers.  Again,
        // this code assumes that sizeof(Scalar) is the number of
        // bytes used by each Scalar.  It also assumes that
        // Teuchos::Comm has correctly deserialized Scalar in place in
        // avValsC.
        avIndsC = imports(curOffsetInBytes, rowSize * sizeof(GlobalOrdinal));
        avValsC = imports(curOffsetInBytes + rowSize * sizeof(GlobalOrdinal),
                          rowSize * sizeof(Scalar));
        avInds = av_reinterpret_cast<const GlobalOrdinal> (avIndsC);
        avVals = av_reinterpret_cast<const Scalar       > (avValsC);

        combineGlobalValues (myGID, avInds(), avVals(), combineMode);
        curOffsetInBytes += rowSize * SizeOfOrdValPair;
      }
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(curOffsetInBytes != totalNumBytes,
        std::logic_error, "After unpacking and combining all the imports, the "
        "final offset in bytes curOffsetInBytes=" << curOffsetInBytes << " != "
        "total number of bytes totalNumBytes=" << totalNumBytes << ".  Please "
        "report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                         Deprecated methods                              //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsMatrix< SCALAR , LO , GO , NODE >; \
  template RCP< CrsMatrix< SCALAR , LO , GO , NODE > >   \
                CrsMatrix< SCALAR , LO , GO , NODE >::convert< SCALAR > () const;

#define TPETRA_CRSMATRIX_CONVERT_INSTANT(SO,SI,LO,GO,NODE) \
  \
  template RCP< CrsMatrix< SO , LO , GO , NODE > >   \
                CrsMatrix< SI , LO , GO , NODE >::convert< SO > () const;

#define TPETRA_CRSMATRIX_IMPORT_AND_FILL_COMPLETE_INSTANT(SCALAR, LO, GO, NODE) \
  template<>                                                                        \
  RCP<CrsMatrix<SCALAR, LO, GO, NODE> >                                \
  importAndFillCompleteCrsMatrix (const RCP<const CrsMatrix<SCALAR, LO, GO, NODE> >& sourceMatrix, \
                                  const Import<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type>& importer, \
                                  const RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& domainMap, \
                                  const RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& rangeMap,  \
                                                               const RCP<Teuchos::ParameterList>& params);

#define TPETRA_CRSMATRIX_EXPORT_AND_FILL_COMPLETE_INSTANT(SCALAR, LO, GO, NODE) \
  template<>                                                                        \
  RCP<CrsMatrix<SCALAR, LO, GO, NODE> >                                \
  exportAndFillCompleteCrsMatrix (const RCP<const CrsMatrix<SCALAR, LO, GO, NODE> >& sourceMatrix, \
                                  const Export<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type>& exporter, \
                                  const RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& domainMap, \
                                  const RCP<const Map<CrsMatrix<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
                                                               CrsMatrix<SCALAR, LO, GO, NODE>::node_type> >& rangeMap,  \
                                                               const RCP<Teuchos::ParameterList>& params);

#endif
