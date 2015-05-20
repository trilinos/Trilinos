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

#ifndef TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DEF_HPP
#define TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DEF_HPP

/// \file Tpetra_Experimental_BlockCrsMatrix_def.hpp
/// \brief Definition of Tpetra::Experimental::BlockCrsMatrix

#include "Tpetra_Experimental_BlockCrsMatrix_decl.hpp"
#include <Tpetra_Details_PackTraits.hpp>

namespace Tpetra {
namespace Experimental {

  template<class Scalar, class LO, class GO, class Node>
  std::ostream&
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  markLocalErrorAndGetStream ()
  {
    * (this->localError_) = true;
    if ((*errs_).is_null ()) {
      *errs_ = Teuchos::rcp (new std::ostringstream ());
    }
    return **errs_;
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix () :
    dist_object_type (Teuchos::rcp (new map_type ())), // nonnull, so DistObject doesn't throw
    graph_ (Teuchos::rcp (new map_type ()), 0), // FIXME (mfh 16 May 2014) no empty ctor yet
    blockSize_ (static_cast<LO> (0)),
#if defined(HAVE_TPETRACLASSIC_SERIAL) || defined(HAVE_TPETRACLASSIC_TBB) || defined(HAVE_TPETRACLASSIC_THREADPOOL) || defined(HAVE_TPETRACLASSIC_OPENMP)
    ptr_ (NULL),
#endif
    ind_ (NULL),
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    columnPadding_ (0), // no padding by default
    rowMajor_ (true), // row major blocks by default
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()), // ptr to a null ptr
    computedDiagonalGraph_(false)
  {
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    rowMeshMap_ (* (graph.getRowMap ())),
    blockSize_ (blockSize),
#if defined(HAVE_TPETRACLASSIC_SERIAL) || defined(HAVE_TPETRACLASSIC_TBB) || defined(HAVE_TPETRACLASSIC_THREADPOOL) || defined(HAVE_TPETRACLASSIC_OPENMP)
    ptr_ (NULL), // to be initialized below
#endif
    ind_ (NULL), // to be initialized below
    val_ (NULL), // to be initialized below
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    columnPadding_ (0), // no padding by default
    rowMajor_ (true), // row major blocks by default
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()), // ptr to a null ptr
    computedDiagonalGraph_(false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

    graphRCP_ = Teuchos::rcpFromRef(graph_);

    // Trick to test whether LO is nonpositive, without a compiler
    // warning in case LO is unsigned (which is generally a bad idea
    // anyway).  I don't promise that the trick works, but it
    // generally does with gcc at least, in my experience.
    const bool blockSizeIsNonpositive = (blockSize + 1 <= 1);
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSizeIsNonpositive, std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input blockSize = " << blockSize <<
      " <= 0.  The block size must be positive.");

    domainPointMap_ = BMV::makePointMap (* (graph.getDomainMap ()), blockSize);
    rangePointMap_ = BMV::makePointMap (* (graph.getRangeMap ()), blockSize);

#if defined(HAVE_TPETRACLASSIC_SERIAL) || defined(HAVE_TPETRACLASSIC_TBB) || defined(HAVE_TPETRACLASSIC_THREADPOOL) || defined(HAVE_TPETRACLASSIC_OPENMP)
    ptr_ = graph.getNodeRowPtrs ().getRawPtr ();
#else
    {
      typedef typename crs_graph_type::local_graph_type::row_map_type row_map_type;
      typedef typename row_map_type::HostMirror::non_const_type nc_host_row_map_type;

      row_map_type ptr_d = graph.getLocalGraph ().row_map;
      // FIXME (mfh 23 Mar 2015) Once we write a Kokkos kernel for the
      // mat-vec, we won't need a host version of this.
      nc_host_row_map_type ptr_h_nc = Kokkos::create_mirror_view (ptr_d);
      Kokkos::deep_copy (ptr_h_nc, ptr_d);
      ptr_ = ptr_h_nc;
    }
#endif

    ind_ = graph.getNodePackedIndices ().getRawPtr ();
    valView_.resize (graph.getNodeNumEntries () * offsetPerBlock ());
    val_ = valView_.getRawPtr ();
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const map_type& domainPointMap,
                  const map_type& rangePointMap,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    rowMeshMap_ (* (graph.getRowMap ())),
    domainPointMap_ (domainPointMap),
    rangePointMap_ (rangePointMap),
    blockSize_ (blockSize),
#if defined(HAVE_TPETRACLASSIC_SERIAL) || defined(HAVE_TPETRACLASSIC_TBB) || defined(HAVE_TPETRACLASSIC_THREADPOOL) || defined(HAVE_TPETRACLASSIC_OPENMP)
    ptr_ (NULL), // to be initialized below
#endif
    ind_ (NULL), // to be initialized below
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    columnPadding_ (0), // no padding by default
    rowMajor_ (true), // row major blocks by default
    localError_ (new bool (false)),
    errs_ (new Teuchos::RCP<std::ostringstream> ()), // ptr to a null ptr
    computedDiagonalGraph_(false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

    graphRCP_ = Teuchos::rcpFromRef(graph_);

    // Trick to test whether LO is nonpositive, without a compiler
    // warning in case LO is unsigned (which is generally a bad idea
    // anyway).  I don't promise that the trick works, but it
    // generally does with gcc at least, in my experience.
    const bool blockSizeIsNonpositive = (blockSize + 1 <= 1);
    TEUCHOS_TEST_FOR_EXCEPTION(
      blockSizeIsNonpositive, std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input blockSize = " << blockSize <<
      " <= 0.  The block size must be positive.");

#if defined(HAVE_TPETRACLASSIC_SERIAL) || defined(HAVE_TPETRACLASSIC_TBB) || defined(HAVE_TPETRACLASSIC_THREADPOOL) || defined(HAVE_TPETRACLASSIC_OPENMP)
    ptr_ = graph.getNodeRowPtrs ().getRawPtr ();
#else
    {
      typedef typename crs_graph_type::local_graph_type::row_map_type row_map_type;
      typedef typename row_map_type::HostMirror::non_const_type nc_host_row_map_type;

      row_map_type ptr_d = graph.getLocalGraph ().row_map;
      // FIXME (mfh 23 Mar 2015) Once we write a Kokkos kernel for the
      // mat-vec, we won't need a host version of this.
      nc_host_row_map_type ptr_h_nc = Kokkos::create_mirror_view (ptr_d);
      Kokkos::deep_copy (ptr_h_nc, ptr_d);
      ptr_ = ptr_h_nc;
    }
#endif
    ind_ = graph.getNodePackedIndices ().getRawPtr ();
    valView_.resize (graph.getNodeNumEntries () * offsetPerBlock ());
    val_ = valView_.getRawPtr ();
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getDomainMap () const
  { // Copy constructor of map_type does a shallow copy.
    // We're only returning by RCP for backwards compatibility.
    return Teuchos::rcp (new map_type (domainPointMap_));
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getRangeMap () const
  { // Copy constructor of map_type does a shallow copy.
    // We're only returning by RCP for backwards compatibility.
    return Teuchos::rcp (new map_type (rangePointMap_));
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getRowMap () const
  {
    return graph_.getRowMap();
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const typename BlockCrsMatrix<Scalar, LO, GO, Node>::map_type>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getColMap () const
  {
    return graph_.getColMap();
  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumRows() const
  {
    return graph_.getGlobalNumRows();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeNumRows() const
  {
    return graph_.getNodeNumRows();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeMaxNumRowEntries() const
  {
    return graph_.getNodeMaxNumRowEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  apply (const mv_type& X,
         mv_type& Y,
         Teuchos::ETransp mode,
         Scalar alpha,
         Scalar beta) const
  {
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    TEUCHOS_TEST_FOR_EXCEPTION(
      mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS && mode != Teuchos::CONJ_TRANS,
      std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::apply: "
      "Invalid 'mode' argument.  Valid values are Teuchos::NO_TRANS, "
      "Teuchos::TRANS, and Teuchos::CONJ_TRANS.");

    BMV X_view;
    BMV Y_view;
    const LO blockSize = getBlockSize ();
    try {
      X_view = BMV (X, * (graph_.getDomainMap ()), blockSize);
      Y_view = BMV (Y, * (graph_.getRangeMap ()), blockSize);
    }
    catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: Either the input MultiVector X or the output MultiVector Y "
        "cannot be viewed as a BlockMultiVector, given this BlockCrsMatrix's "
        "graph.  BlockMultiVector's constructor threw the following exception: "
        << e.what ());
    }

    try {
      // mfh 16 May 2014: Casting 'this' to nonconst is icky, but it's
      // either that or mark fields of this class as 'mutable'.  The
      // problem is that applyBlock wants to do lazy initialization of
      // temporary block multivectors.
      const_cast<this_type*> (this)->applyBlock (X_view, Y_view, mode, alpha, beta);
    } catch (std::invalid_argument& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: The implementation method applyBlock complained about having "
        "an invalid argument.  It reported the following: " << e.what ());
    } catch (std::logic_error& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: The implementation method applyBlock complained about a "
        "possible bug in its implementation.  It reported the following: "
        << e.what ());
    } catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: The implementation method applyBlock threw an exception which "
        "is neither std::invalid_argument nor std::logic_error, but is a "
        "subclass of std::exception.  It reported the following: "
        << e.what ());
    } catch (...) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
        "apply: The implementation method applyBlock threw an exception which "
        "is not an instance of a subclass of std::exception.  This probably "
        "indicates a bug.  Please report this to the Tpetra developers.");
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlock (const BlockMultiVector<Scalar, LO, GO, Node>& X,
              BlockMultiVector<Scalar, LO, GO, Node>& Y,
              Teuchos::ETransp mode,
              const Scalar alpha,
              const Scalar beta)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getBlockSize () != Y.getBlockSize (), std::invalid_argument,
      "Tpetra::Experimental::BlockCrsMatrix::applyBlock: "
      "X and Y have different block sizes.  X.getBlockSize() = "
      << X.getBlockSize () << " != Y.getBlockSize() = "
      << Y.getBlockSize () << ".");

    if (mode == Teuchos::NO_TRANS) {
      applyBlockNoTrans (X, Y, alpha, beta);
    } else if (mode == Teuchos::TRANS || mode == Teuchos::CONJ_TRANS) {
      applyBlockTrans (X, Y, mode, alpha, beta);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "Tpetra::Experimental::BlockCrsMatrix::"
        "applyBlock: Invalid 'mode' argument.  Valid values are "
        "Teuchos::NO_TRANS, Teuchos::TRANS, and Teuchos::CONJ_TRANS.");
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  setAllToScalar (const Scalar& alpha)
  {
    const LO numLocalMeshRows = static_cast<LO> (rowMeshMap_.getNodeNumElements ());
    for (LO lclRow = 0; lclRow < numLocalMeshRows; ++lclRow) {
      const size_t meshBeg = ptr_[lclRow];
      const size_t meshEnd = ptr_[lclRow+1];
      for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
        little_block_type A_cur = getNonConstLocalBlockFromAbsOffset (absBlkOff);
        A_cur.fill (static_cast<impl_scalar_type> (alpha));
      }
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  replaceLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptr_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t hint = 0; // Guess for the relative offset into the current row
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid column indices in colInds

    for (LO k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset =
        findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        A_old.assign (A_new);
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const
  {

    const map_type& rowMap = * (graph_.getRowMap());
    const map_type& colMap = * (graph_.getColMap ());

    const size_t myNumRows = rowMeshMap_.getNodeNumElements();
    if (static_cast<size_t> (offsets.size ()) != myNumRows) {
      offsets.resize (static_cast<size_t> (myNumRows));
    }

#ifdef HAVE_TPETRA_DEBUG
    bool allRowMapDiagEntriesInColMap = true;
    bool allDiagEntriesFound = true;
#endif // HAVE_TPETRA_DEBUG

    for (size_t r = 0; r < myNumRows; ++r) {
      const GO rgid = rowMap.getGlobalElement (r);
      const LO rlid = colMap.getLocalElement (rgid);

#ifdef HAVE_TPETRA_DEBUG
      if (rlid == Teuchos::OrdinalTraits<LO>::invalid ()) {
        allRowMapDiagEntriesInColMap = false;
      }
#endif // HAVE_TPETRA_DEBUG

      if (rlid != Teuchos::OrdinalTraits<LO>::invalid ()) {
        RowInfo rowinfo = graph_.getRowInfo (r);
        if (rowinfo.numEntries > 0) {
          offsets[r] = graph_.findLocalIndex (rowinfo, rlid);
        }
        else {
          offsets[r] = Teuchos::OrdinalTraits<size_t>::invalid ();
#ifdef HAVE_TPETRA_DEBUG
          allDiagEntriesFound = false;
#endif // HAVE_TPETRA_DEBUG
        }
      }
    }

#ifdef HAVE_TPETRA_DEBUG
    using Teuchos::reduceAll;
    using std::endl;
    const char tfecfFuncName[] = "getLocalDiagOffsets";

    const bool localSuccess =
      allRowMapDiagEntriesInColMap && allDiagEntriesFound;
    int localResults[3];
    localResults[0] = allRowMapDiagEntriesInColMap ? 1 : 0;
    localResults[1] = allDiagEntriesFound ? 1 : 0;
    // min-all-reduce will compute least rank of all the processes
    // that didn't succeed.
    localResults[2] =
      ! localSuccess ? getComm ()->getRank () : getComm ()->getSize ();
    int globalResults[3];
    globalResults[0] = 0;
    globalResults[1] = 0;
    globalResults[2] = 0;
    reduceAll<int, int> (* (getComm ()), Teuchos::REDUCE_MIN,
                         3, localResults, globalResults);
    if (globalResults[0] == 0 || globalResults[1] == 0) {
      std::ostringstream os; // build error message
      const bool both =
        globalResults[0] == 0 && globalResults[1] == 0;
      os << ": At least one process (including Process " << globalResults[2]
         << ") had the following issue" << (both ? "s" : "") << ":" << endl;
      if (globalResults[0] == 0) {
        os << "  - The column Map does not contain at least one diagonal entry "
          "of the matrix." << endl;
      }
      if (globalResults[1] == 0) {
        os << "  - There is a row on that / those process(es) that does not "
          "contain a diagonal entry." << endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::runtime_error, os.str());
    }
#endif // HAVE_TPETRA_DEBUG
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  computeDiagonalGraph ()
  {
    using Teuchos::rcp;

    if (computedDiagonalGraph_) {
      // FIXME (mfh 12 Aug 2014) Consider storing the "diagonal graph"
      // separately from the matrix.  It should really go in the
      // preconditioner, not here.  We could do this by adding a
      // method that accepts a nonconst diagonal graph, and updates
      // it.  btw it would probably be a better idea to use a
      // BlockMultiVector to store the diagonal, not a graph.
      return;
    }

    const size_t maxDiagEntPerRow = 1;
    // NOTE (mfh 12 Aug 2014) We could also pass in the column Map
    // here.  However, we still would have to do LID->GID lookups to
    // make sure that we are using the correct diagonal column
    // indices, so it probably wouldn't help much.
    diagonalGraph_ =
      rcp (new crs_graph_type (graph_.getRowMap (), maxDiagEntPerRow,
                               Tpetra::StaticProfile));
    const map_type& meshRowMap = * (graph_.getRowMap ());

    Teuchos::Array<GO> diagGblColInds (maxDiagEntPerRow);

    for (LO lclRowInd = meshRowMap.getMinLocalIndex ();
         lclRowInd <= meshRowMap.getMaxLocalIndex (); ++lclRowInd) {
      const GO gblRowInd = meshRowMap.getGlobalElement (lclRowInd);
      diagGblColInds[0] = gblRowInd;
      diagonalGraph_->insertGlobalIndices (gblRowInd, diagGblColInds ());
    }
    diagonalGraph_->fillComplete (graph_.getDomainMap (),
                                  graph_.getRangeMap ());
    computedDiagonalGraph_ = true;
  }

  template <class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<CrsGraph<LO, GO, Node> >
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getDiagonalGraph () const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! computedDiagonalGraph_, std::runtime_error, "Tpetra::Experimental::"
      "BlockCrsMatrix::getDiagonalGraph: You must call computeDiagonalGraph() "
      "before calling this method.");
    return diagonalGraph_;
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  localGaussSeidel (const BlockMultiVector<Scalar, LO, GO, Node>& B,
                    BlockMultiVector<Scalar, LO, GO, Node>& X,
                    BlockCrsMatrix<Scalar, LO, GO, Node> & factorizedDiagonal,
                    const int * factorizationPivots,
                    const Scalar omega,
                    const ESweepDirection direction) const
  {
    const LO numLocalMeshRows =
      static_cast<LO> (rowMeshMap_.getNodeNumElements ());
    const LO numVecs = static_cast<LO> (X.getNumVectors ());

    // If using (new) Kokkos, replace localMem with thread-local
    // memory.  Note that for larger block sizes, this will affect the
    // two-level parallelization.  Look to Stokhos for best practice
    // on making this fast for GPUs.
    const LO blockSize = getBlockSize ();
    Teuchos::Array<impl_scalar_type> localMem (blockSize);
    Teuchos::Array<impl_scalar_type> localMat (blockSize*blockSize);
    little_vec_type X_lcl (localMem.getRawPtr (), blockSize, 1);

    const LO * columnIndices;
    Scalar * Dmat;
    LO numIndices;

    // FIXME (mfh 12 Aug 2014) This probably won't work if LO is unsigned.
    LO rowBegin = 0, rowEnd = 0, rowStride = 0;
    if (direction == Forward) {
      rowBegin = 1;
      rowEnd = numLocalMeshRows+1;
      rowStride = 1;
    }
    else if (direction == Backward) {
      rowBegin = numLocalMeshRows;
      rowEnd = 0;
      rowStride = -1;
    }
    else if (direction == Symmetric) {
      this->localGaussSeidel (B, X, factorizedDiagonal, factorizationPivots, omega, Forward);
      this->localGaussSeidel (B, X, factorizedDiagonal, factorizationPivots, omega, Backward);
      return;
    }

    const Scalar one_minus_omega = Teuchos::ScalarTraits<Scalar>::one()-omega;
    const Scalar     minus_omega = -omega;

    if (numVecs == 1) {
      for (LO lclRow = rowBegin; lclRow != rowEnd; lclRow += rowStride) {
        const LO actlRow = lclRow - 1;

        little_vec_type B_cur = B.getLocalBlock (actlRow, 0);
        X_lcl.assign (B_cur);
        X_lcl.scale (omega);

        const size_t meshBeg = ptr_[actlRow];
        const size_t meshEnd = ptr_[actlRow+1];
        for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
          const LO meshCol = ind_[absBlkOff];
          const_little_block_type A_cur =
            getConstLocalBlockFromAbsOffset (absBlkOff);

          little_vec_type X_cur = X.getLocalBlock (meshCol, 0);

          // X_lcl += alpha*A_cur*X_cur
          const Scalar alpha = meshCol == actlRow ? one_minus_omega : minus_omega;
          X_lcl.matvecUpdate (alpha, A_cur, X_cur);
        } // for each entry in the current local row of the matrx

        factorizedDiagonal.getLocalRowView (actlRow, columnIndices,
                                            Dmat, numIndices);
        little_block_type D_lcl =
          getNonConstLocalBlockFromInput (reinterpret_cast<impl_scalar_type*> (Dmat), 0);

        D_lcl.solve (X_lcl, &factorizationPivots[actlRow*blockSize_]);
        little_vec_type X_update = X.getLocalBlock (actlRow, 0);
        X_update.assign(X_lcl);
      } // for each local row of the matrix
    }
    else {
      for (LO lclRow = rowBegin; lclRow != rowEnd; lclRow += rowStride) {
        for (LO j = 0; j < numVecs; ++j) {
          LO actlRow = lclRow-1;

          little_vec_type B_cur = B.getLocalBlock (actlRow, j);
          X_lcl.assign (B_cur);
          X_lcl.scale (omega);

          const size_t meshBeg = ptr_[actlRow];
          const size_t meshEnd = ptr_[actlRow+1];
          for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
            const LO meshCol = ind_[absBlkOff];
            const_little_block_type A_cur =
              getConstLocalBlockFromAbsOffset (absBlkOff);

            little_vec_type X_cur = X.getLocalBlock (meshCol, j);

            // X_lcl += alpha*A_cur*X_cur
            const Scalar alpha = meshCol == actlRow ? one_minus_omega : minus_omega;
            X_lcl.matvecUpdate (alpha, A_cur, X_cur);
          } // for each entry in the current local row of the matrx

          factorizedDiagonal.getLocalRowView (actlRow, columnIndices,
                                              Dmat, numIndices);
          little_block_type D_lcl =
            getNonConstLocalBlockFromInput (reinterpret_cast<impl_scalar_type*> (Dmat), 0);

          D_lcl.solve (X_lcl, &factorizationPivots[actlRow*blockSize_]);

          little_vec_type X_update = X.getLocalBlock (actlRow, j);
          X_update.assign(X_lcl);
        } // for each entry in the current local row of the matrix
      } // for each local row of the matrix
    }
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  gaussSeidelCopy (MultiVector<Scalar,LO,GO,Node> &X,
                   const MultiVector<Scalar,LO,GO,Node> &B,
                   const MultiVector<Scalar,LO,GO,Node> &D,
                   const Scalar& dampingFactor,
                   const ESweepDirection direction,
                   const int numSweeps,
                   const bool zeroInitialGuess) const
  {
    // FIXME (mfh 12 Aug 2014) This method has entirely the wrong
    // interface for block Gauss-Seidel.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "gaussSeidelCopy: Not implemented.");
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  reorderedGaussSeidelCopy (MultiVector<Scalar,LO,GO,Node>& X,
                            const MultiVector<Scalar,LO,GO,Node>& B,
                            const MultiVector<Scalar,LO,GO,Node>& D,
                            const ArrayView<LO>& rowIndices,
                            const Scalar& dampingFactor,
                            const ESweepDirection direction,
                            const int numSweeps,
                            const bool zeroInitialGuess) const
  {
    // FIXME (mfh 12 Aug 2014) This method has entirely the wrong
    // interface for block Gauss-Seidel.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "reorderedGaussSeidelCopy: Not implemented.");
  }

  template <class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar,LO,GO,Node>::
  getLocalDiagCopy (BlockCrsMatrix<Scalar,LO,GO,Node>& diag,
                    const Teuchos::ArrayView<const size_t>& offsets) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero ();

    const size_t myNumRows = rowMeshMap_.getNodeNumElements();
    const LO* columnIndices;
    Scalar* vals;
    LO numColumns;
    Teuchos::Array<LO> cols(1);

    // FIXME (mfh 12 Aug 2014) Should use a "little block" for this instead.
    Teuchos::Array<Scalar> zeroMat (blockSize_*blockSize_, ZERO);
    for (size_t i = 0; i < myNumRows; ++i) {
      cols[0] = i;
      if (offsets[i] == Teuchos::OrdinalTraits<size_t>::invalid ()) {
        diag.replaceLocalValues (i, cols.getRawPtr (), zeroMat.getRawPtr (), 1);
      }
      else {
        getLocalRowView (i, columnIndices, vals, numColumns);
        diag.replaceLocalValues (i, cols.getRawPtr(), &vals[offsets[i]*blockSize_*blockSize_], 1);
      }
    }
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  absMaxLocalValues (const LO localRowInd,
                     const LO colInds[],
                     const Scalar vals[],
                     const LO numColInds) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptr_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t hint = 0; // Guess for the relative offset into the current row
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid column indices in colInds

    for (LO k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset =
        findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        A_old.absmax (A_new);
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  sumIntoLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptr_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t hint = 0; // Guess for the relative offset into the current row
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid column indices in colInds

    for (LO k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset =
        findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        A_old.update (static_cast<impl_scalar_type> (STS::one ()), A_new);
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowView (const LO localRowInd,
                   const LO*& colInds,
                   Scalar*& vals,
                   LO& numInds) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      colInds = NULL;
      vals = NULL;
      numInds = 0;
      return Teuchos::OrdinalTraits<LO>::invalid ();
    }
    else {
      const size_t absBlockOffsetStart = ptr_[localRowInd];
      colInds = ind_ + absBlockOffsetStart;

      impl_scalar_type* const vOut = val_ + absBlockOffsetStart * offsetPerBlock ();
      vals = reinterpret_cast<Scalar*> (vOut);

      numInds = ptr_[localRowInd + 1] - absBlockOffsetStart;
      return 0; // indicates no error
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowCopy (LO LocalRow,
                   const Teuchos::ArrayView<LO>& Indices,
                   const Teuchos::ArrayView<Scalar>& Values,
                   size_t &NumEntries) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "getLocalRowCopy: Copying is not implemented. You should use a view.");
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowOffsets (const LO localRowInd,
                      ptrdiff_t offsets[],
                      const LO colInds[],
                      const LO numColInds) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We got no offsets, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }

    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t hint = 0; // Guess for the relative offset into the current row
    LO validCount = 0; // number of valid column indices in colInds

    for (LO k = 0; k < numColInds; ++k) {
      const size_t relBlockOffset =
        findRelOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (relBlockOffset != STINV) {
        offsets[k] = relBlockOffset;
        hint = relBlockOffset + 1;
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  replaceLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptr_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset = offsets[k];
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        A_old.assign (A_new);
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  absMaxLocalValuesByOffsets (const LO localRowInd,
                              const ptrdiff_t offsets[],
                              const Scalar vals[],
                              const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptr_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset = offsets[k];
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        A_old.absmax (A_new);
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  sumIntoLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const
  {
    if (! rowMeshMap_.isNodeLocalElement (localRowInd)) {
      // We modified no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }
    const impl_scalar_type* const vIn = reinterpret_cast<const impl_scalar_type*> (vals);

    const size_t absRowBlockOffset = ptr_[localRowInd];
    const size_t perBlockSize = static_cast<LO> (offsetPerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t pointOffset = 0; // Current offset into input values
    LO validCount = 0; // number of valid offsets

    for (LO k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
      const size_t relBlockOffset = offsets[k];
      if (relBlockOffset != STINV) {
        const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
        little_block_type A_old =
          getNonConstLocalBlockFromAbsOffset (absBlockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vIn, pointOffset);
        A_old.update (static_cast<impl_scalar_type> (STS::one ()), A_new);
        ++validCount;
      }
    }
    return validCount;
  }


  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNumEntriesInLocalRow (const LO localRowInd) const
  {
    const size_t numEntInGraph = graph_.getNumEntriesInLocalRow (localRowInd);
    if (numEntInGraph == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      return static_cast<LO> (0); // the calling process doesn't have that row
    } else {
      return static_cast<LO> (numEntInGraph);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlockTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                   BlockMultiVector<Scalar, LO, GO, Node>& Y,
                   const Teuchos::ETransp mode,
                   const Scalar alpha,
                   const Scalar beta)
  {
    (void) X;
    (void) Y;
    (void) mode;
    (void) alpha;
    (void) beta;

    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::apply: "
      "transpose and conjugate transpose modes are not yet implemented.");
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  applyBlockNoTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                     BlockMultiVector<Scalar, LO, GO, Node>& Y,
                     const Scalar alpha,
                     const Scalar beta)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Tpetra::Import<LO, GO, Node> import_type;
    typedef Tpetra::Export<LO, GO, Node> export_type;
    const Scalar zero = STS::zero ();
    const Scalar one = STS::one ();
    RCP<const import_type> import = graph_.getImporter ();
    // "export" is a reserved C++ keyword, so we can't use it.
    RCP<const export_type> theExport = graph_.getExporter ();

    // FIXME (mfh 20 May 2014) X.mv_ and Y.mv_ requires a friend
    // declaration, which is useful only for debugging.
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.mv_.getCopyOrView () != Teuchos::View, std::invalid_argument,
      "Tpetra::Experimental::BlockCrsMatrix::applyBlockNoTrans: "
      "The input BlockMultiVector X has deep copy semantics, "
      "not view semantics (as it should).");
    TEUCHOS_TEST_FOR_EXCEPTION(
      Y.mv_.getCopyOrView () != Teuchos::View, std::invalid_argument,
      "Tpetra::Experimental::BlockCrsMatrix::applyBlockNoTrans: "
      "The output BlockMultiVector Y has deep copy semantics, "
      "not view semantics (as it should).");

    if (alpha == zero) {
      if (beta == zero) {
        Y.putScalar (zero); // replace Inf or NaN (BLAS rules)
      } else if (beta != one) {
        Y.scale (beta);
      }
    } else { // alpha != 0
      const BMV* X_colMap = NULL;
      if (import.is_null ()) {
        try {
          X_colMap = &X;
        } catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
            "applyBlockNoTrans:" << std::endl << "Tpetra::MultiVector::"
            "operator= threw an exception: " << std::endl << e.what ());
        }
      } else {
        // X_colMap_ is a pointer to a pointer to BMV.  Ditto for
        // Y_rowMap_ below.  This lets us do lazy initialization
        // correctly with view semantics of BlockCrsMatrix.  All views
        // of this BlockCrsMatrix have the same outer pointer.  That
        // way, we can set the inner pointer in one view, and all
        // other views will see it.
        if ((*X_colMap_).is_null () ||
            (**X_colMap_).getNumVectors () != X.getNumVectors () ||
            (**X_colMap_).getBlockSize () != X.getBlockSize ()) {
          *X_colMap_ = rcp (new BMV (* (graph_.getColMap ()), getBlockSize (),
                                     static_cast<LO> (X.getNumVectors ())));
        }
        (**X_colMap_).doImport (X, *import, Tpetra::REPLACE);
        try {
          X_colMap = &(**X_colMap_);
        } catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
            "applyBlockNoTrans:" << std::endl << "Tpetra::MultiVector::"
            "operator= threw an exception: " << std::endl << e.what ());
        }
      }

      BMV* Y_rowMap = NULL;
      if (theExport.is_null ()) {
        Y_rowMap = &Y;
      } else if ((*Y_rowMap_).is_null () ||
                 (**Y_rowMap_).getNumVectors () != Y.getNumVectors () ||
                 (**Y_rowMap_).getBlockSize () != Y.getBlockSize ()) {
        *Y_rowMap_ = rcp (new BMV (* (graph_.getRowMap ()), getBlockSize (),
                                   static_cast<LO> (X.getNumVectors ())));
        try {
          Y_rowMap = &(**Y_rowMap_);
        } catch (std::exception& e) {
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
            "applyBlockNoTrans:" << std::endl << "Tpetra::MultiVector::"
            "operator= threw an exception: " << std::endl << e.what ());
        }
      }

      localApplyBlockNoTrans (*X_colMap, *Y_rowMap, alpha, beta);

      if (! theExport.is_null ()) {
        Y.doExport (*Y_rowMap, *theExport, Tpetra::REPLACE);
      }
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  localApplyBlockNoTrans (const BlockMultiVector<Scalar, LO, GO, Node>& X,
                          BlockMultiVector<Scalar, LO, GO, Node>& Y,
                          const Scalar alpha,
                          const Scalar beta)
  {
    // If using (new) Kokkos, prefer Kokkos::ArithTraits to
    // Teuchos::ScalarTraits.
    const Scalar zero = STS::zero ();
    const Scalar one = STS::one ();
    const LO numLocalMeshRows =
      static_cast<LO> (rowMeshMap_.getNodeNumElements ());
    const LO numVecs = static_cast<LO> (X.getNumVectors ());

    // If using (new) Kokkos, replace localMem with thread-local
    // memory.  Note that for larger block sizes, this will affect the
    // two-level parallelization.  Look to Stokhos for best practice
    // on making this fast for GPUs.
    const LO blockSize = getBlockSize ();
    Teuchos::Array<impl_scalar_type> localMem (blockSize);
    little_vec_type Y_lcl (localMem.getRawPtr (), blockSize, 1);

    if (numVecs == 1) {
      for (LO lclRow = 0; lclRow < numLocalMeshRows; ++lclRow) {
        little_vec_type Y_cur = Y.getLocalBlock (lclRow, 0);

        if (beta == zero) {
          Y_lcl.fill (zero);
        } else if (beta == one) {
          Y_lcl.assign (Y_cur);
        } else {
          Y_lcl.assign (Y_cur);
          Y_lcl.scale (beta);
        }

        const size_t meshBeg = ptr_[lclRow];
        const size_t meshEnd = ptr_[lclRow+1];
        for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
          const LO meshCol = ind_[absBlkOff];
          const_little_block_type A_cur =
            getConstLocalBlockFromAbsOffset (absBlkOff);
          little_vec_type X_cur = X.getLocalBlock (meshCol, 0);
          // Y_lcl += alpha*A_cur*X_cur
          Y_lcl.matvecUpdate (alpha, A_cur, X_cur);
        } // for each entry in the current local row of the matrx

        Y_cur.assign (Y_lcl);
      } // for each local row of the matrix
    }
    else {
      for (LO lclRow = 0; lclRow < numLocalMeshRows; ++lclRow) {
        for (LO j = 0; j < numVecs; ++j) {
          little_vec_type Y_cur = Y.getLocalBlock (lclRow, j);

          if (beta == zero) {
            Y_lcl.fill (zero);
          } else if (beta == one) {
            Y_lcl.assign (Y_cur);
          } else {
            Y_lcl.assign (Y_cur);
            Y_lcl.scale (beta);
          }

          const size_t meshBeg = ptr_[lclRow];
          const size_t meshEnd = ptr_[lclRow+1];
          for (size_t absBlkOff = meshBeg; absBlkOff < meshEnd; ++absBlkOff) {
            const LO meshCol = ind_[absBlkOff];
            const_little_block_type A_cur =
              getConstLocalBlockFromAbsOffset (absBlkOff);
            little_vec_type X_cur = X.getLocalBlock (meshCol, j);
            // Y_lcl += alpha*A_cur*X_cur
            Y_lcl.matvecUpdate (alpha, A_cur, X_cur);
          } // for each entry in the current local row of the matrix

          Y_cur.assign (Y_lcl);
        } // for each entry in the current row of Y
      } // for each local row of the matrix
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  findRelOffsetOfColumnIndex (const LO localRowIndex,
                              const LO colIndexToFind,
                              const size_t hint) const
  {
    const size_t absStartOffset = ptr_[localRowIndex];
    const size_t absEndOffset = ptr_[localRowIndex+1];
    const size_t numEntriesInRow = absEndOffset - absStartOffset;

    // If the hint was correct, then the hint is the offset to return.
    if (hint < numEntriesInRow && ind_[absStartOffset+hint] == colIndexToFind) {
      // Always return the offset relative to the current row.
      return hint;
    }

    // The hint was wrong, so we must search for the given column
    // index in the column indices for the given row.
    size_t relOffset = Teuchos::OrdinalTraits<size_t>::invalid ();

    // We require that the graph have sorted rows.  However, binary
    // search only pays if the current row is longer than a certain
    // amount.  We set this to 32, but you might want to tune this.
    const size_t maxNumEntriesForLinearSearch = 32;
    if (numEntriesInRow > maxNumEntriesForLinearSearch) {
      // Use binary search.  It would probably be better for us to
      // roll this loop by hand.  If we wrote it right, a smart
      // compiler could perhaps use conditional loads and avoid
      // branches (according to Jed Brown on May 2014).
      const LO* beg = ind_ + absStartOffset;
      const LO* end = ind_ + absEndOffset;
      std::pair<const LO*, const LO*> p =
        std::equal_range (beg, end, colIndexToFind);
      if (p.first != p.second) {
        // offset is relative to the current row
        relOffset = static_cast<size_t> (p.first - beg);
      }
    }
    else { // use linear search
      for (size_t k = 0; k < numEntriesInRow; ++k) {
        if (colIndexToFind == ind_[absStartOffset + k]) {
          relOffset = k; // offset is relative to the current row
          break;
        }
      }
    }

    return relOffset;
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  offsetPerBlock () const
  {
    const LO numRows = blockSize_;

    LO numCols = blockSize_;
    if (columnPadding_ > 0) { // Column padding == 0 means no padding.
      const LO numColsRoundedDown = (blockSize_ / columnPadding_) * columnPadding_;
      numCols = (numColsRoundedDown < numCols) ?
        (numColsRoundedDown + columnPadding_) :
        numColsRoundedDown;
    }
    return numRows * numCols;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::const_little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getConstLocalBlockFromInput (const impl_scalar_type* val,
                               const size_t pointOffset) const
  {
    if (rowMajor_) {
      const size_t rowStride = (columnPadding_ == 0) ?
        static_cast<size_t> (blockSize_) : static_cast<size_t> (columnPadding_);
      return const_little_block_type (val + pointOffset, blockSize_, rowStride, 1);
    } else {
      return const_little_block_type (val + pointOffset, blockSize_, 1, blockSize_);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNonConstLocalBlockFromInput (impl_scalar_type* val,
                                  const size_t pointOffset) const
  {
    if (rowMajor_) {
      const size_t rowStride = (columnPadding_ == 0) ?
        static_cast<size_t> (blockSize_) : static_cast<size_t> (columnPadding_);
      return little_block_type (val + pointOffset, blockSize_, rowStride, 1);
    } else {
      return little_block_type (val + pointOffset, blockSize_, 1, blockSize_);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::const_little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getConstLocalBlockFromAbsOffset (const size_t absBlockOffset) const
  {
    if (absBlockOffset >= ptr_[rowMeshMap_.getNodeNumElements ()]) {
      // An empty block signifies an error.  We don't expect to see
      // this error in correct code, but it's helpful for avoiding
      // memory corruption in case there is a bug.
      return const_little_block_type (NULL, 0, 0, 0);
    } else {
      const size_t absPointOffset = absBlockOffset * offsetPerBlock ();
      return getConstLocalBlockFromInput (val_, absPointOffset);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNonConstLocalBlockFromAbsOffset (const size_t absBlockOffset) const
  {
    if (absBlockOffset >= ptr_[rowMeshMap_.getNodeNumElements ()]) {
      // An empty block signifies an error.  We don't expect to see
      // this error in correct code, but it's helpful for avoiding
      // memory corruption in case there is a bug.
      return little_block_type (NULL, 0, 0, 0);
    } else {
      const size_t absPointOffset = absBlockOffset * offsetPerBlock ();
      return getNonConstLocalBlockFromInput (const_cast<impl_scalar_type*> (val_),
                                             absPointOffset);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalBlock (const LO localRowInd, const LO localColInd) const
  {
    const size_t absRowBlockOffset = ptr_[localRowInd];

    size_t hint = 0;
    const size_t relBlockOffset =
      findRelOffsetOfColumnIndex (localRowInd, localColInd, hint);

    if (relBlockOffset != Teuchos::OrdinalTraits<size_t>::invalid ()) {
      const size_t absBlockOffset = absRowBlockOffset + relBlockOffset;
      return getNonConstLocalBlockFromAbsOffset (absBlockOffset);
    }
    else {
      return little_block_type (NULL, 0, 0, 0);
    }
  }

  // template<class Scalar, class LO, class GO, class Node>
  // void
  // BlockCrsMatrix<Scalar, LO, GO, Node>::
  // clearLocalErrorStateAndStream ()
  // {
  //   typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
  //   * (const_cast<this_type*> (this)->localError_) = false;
  //   *errs_ = Teuchos::null;
  // }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  checkSizes (const Tpetra::SrcDistObject& source)
  {
    using std::endl;
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    const this_type* src = dynamic_cast<const this_type* > (&source);

    if (src == NULL) {
      std::ostream& err = markLocalErrorAndGetStream ();
      err << "checkSizes: The source object of the Import or Export "
        "must be a BlockCrsMatrix with the same template parameters as the "
        "target object." << endl;
    }
    else {
      // Use a string of ifs, not if-elseifs, because we want to know
      // all the errors.
      if (src->getBlockSize () != this->getBlockSize ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source and target objects of the Import or "
            << "Export must have the same block sizes.  The source's block "
            << "size = " << src->getBlockSize () << " != the target's block "
            << "size = " << this->getBlockSize () << "." << endl;
      }
      if (! src->graph_.isFillComplete ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source object of the Import or Export is "
          "not fill complete.  Both source and target objects must be fill "
          "complete." << endl;
      }
      if (! this->graph_.isFillComplete ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The target object of the Import or Export is "
          "not fill complete.  Both source and target objects must be fill "
          "complete." << endl;
      }
      if (src->graph_.getColMap ().is_null ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The source object of the Import or Export does "
          "not have a column Map.  Both source and target objects must have "
          "column Maps." << endl;
      }
      if (this->graph_.getColMap ().is_null ()) {
        std::ostream& err = markLocalErrorAndGetStream ();
        err << "checkSizes: The target object of the Import or Export does "
          "not have a column Map.  Both source and target objects must have "
          "column Maps." << endl;
      }
    }

    return ! (* (this->localError_));
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  copyAndPermute (const Tpetra::SrcDistObject& source,
                  size_t numSameIDs,
                  const Teuchos::ArrayView<const LO>& permuteToLIDs,
                  const Teuchos::ArrayView<const LO>& permuteFromLIDs)
  {
    using std::endl;
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    const bool debug = false;

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": copyAndPermute: "
         << "numSameIDs: " << numSameIDs
         << ", permuteToLIDs.size(): " << permuteToLIDs.size ()
         << ", permuteFromLIDs.size(): " << permuteFromLIDs.size ()
         << endl;
      std::cerr << os.str ();
    }

    // There's no communication in this method, so it's safe just to
    // return on error.

    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << "copyAndPermute: The target object of the Import or Export is "
        "already in an error state." << endl;
      return;
    }
    if (permuteToLIDs.size () != permuteFromLIDs.size ()) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << "copyAndPermute: permuteToLIDs.size() = " << permuteToLIDs.size ()
          << " != permuteFromLIDs.size() = " << permuteFromLIDs.size () << "."
          << endl;
      return;
    }

    const this_type* src = dynamic_cast<const this_type* > (&source);
    if (src == NULL) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << "copyAndPermute: The source object of the Import or Export is "
        "either not a BlockCrsMatrix, or does not have the right template "
        "parameters.  checkSizes() should have caught this.  "
        "Please report this bug to the Tpetra developers." << endl;
      return;
    }
    if (* (src->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << "copyAndPermute: The source object of the Import or Export is "
        "already in an error state." << endl;
      return;
    }

    bool lclErr = false;
#ifdef HAVE_TPETRA_DEBUG
    std::set<LO> invalidSrcCopyRows;
    std::set<LO> invalidDstCopyRows;
    std::set<LO> invalidDstCopyCols;
    std::set<LO> invalidDstPermuteCols;
    std::set<LO> invalidPermuteFromRows;
#endif // HAVE_TPETRA_DEBUG

    // Copy the initial sequence of rows that are the same.
    //
    // The two graphs might have different column Maps, so we need to
    // do this using global column indices.  This is purely local, so
    // we only need to check for local sameness of the two column
    // Maps.

#ifdef HAVE_TPETRA_DEBUG
    const map_type& srcRowMap = * (src->graph_.getRowMap ());
#endif // HAVE_TPETRA_DEBUG
    const map_type& dstRowMap = * (this->graph_.getRowMap ());
    const map_type& srcColMap = * (src->graph_.getColMap ());
    const map_type& dstColMap = * (this->graph_.getColMap ());
    const bool canUseLocalColumnIndices = srcColMap.locallySameAs (dstColMap);
    const size_t numPermute = static_cast<size_t> (permuteFromLIDs.size ());

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": copyAndPermute: "
         << "canUseLocalColumnIndices: "
         << (canUseLocalColumnIndices ? "true" : "false")
         << ", numPermute: " << numPermute
         << endl;
      std::cerr << os.str ();
    }

    if (canUseLocalColumnIndices) {
      // Copy local rows that are the "same" in both source and target.
      for (LO localRow = 0; localRow < static_cast<LO> (numSameIDs); ++localRow) {
#ifdef HAVE_TPETRA_DEBUG
        if (! srcRowMap.isNodeLocalElement (localRow)) {
          lclErr = true;
          invalidSrcCopyRows.insert (localRow);
          continue; // skip invalid rows
        }
#endif // HAVE_TPETRA_DEBUG

        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        // If this call fails, that means the mesh row local index is
        // invalid.  That means the Import or Export is invalid somehow.
        LO err = src->getLocalRowView (localRow, lclSrcCols, vals, numEntries);
        if (err != 0) {
          lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
          (void) invalidSrcCopyRows.insert (localRow);
#endif // HAVE_TPETRA_DEBUG
        }
        else {
          err = this->replaceLocalValues (localRow, lclSrcCols, vals, numEntries);
          if (err != numEntries) {
            lclErr = true;
            if (! dstRowMap.isNodeLocalElement (localRow)) {
#ifdef HAVE_TPETRA_DEBUG
              invalidDstCopyRows.insert (localRow);
#endif // HAVE_TPETRA_DEBUG
            }
            else {
              // Once there's an error, there's no sense in saving
              // time, so we check whether the column indices were
              // invalid.  However, only remember which ones were
              // invalid in a debug build, because that might take a
              // lot of space.
              for (LO k = 0; k < numEntries; ++k) {
                if (! dstColMap.isNodeLocalElement (lclSrcCols[k])) {
                  lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                  (void) invalidDstCopyCols.insert (lclSrcCols[k]);
#endif // HAVE_TPETRA_DEBUG
                }
              }
            }
          }
        }
      } // for each "same" local row

      // Copy the "permute" local rows.
      for (size_t k = 0; k < numPermute; ++k) {
        const LO srcLclRow = static_cast<LO> (permuteFromLIDs[k]);
        const LO dstLclRow = static_cast<LO> (permuteToLIDs[k]);

        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        LO err = src->getLocalRowView (srcLclRow, lclSrcCols, vals, numEntries);
        if (err != 0) {
          lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
          invalidPermuteFromRows.insert (srcLclRow);
#endif // HAVE_TPETRA_DEBUG
        }
        else {
          err = this->replaceLocalValues (dstLclRow, lclSrcCols, vals, numEntries);
          if (err != numEntries) {
            lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
            for (LO c = 0; c < numEntries; ++c) {
              if (! dstColMap.isNodeLocalElement (lclSrcCols[c])) {
                invalidDstPermuteCols.insert (lclSrcCols[c]);
              }
            }
#endif // HAVE_TPETRA_DEBUG
          }
        }
      }
    }
    else { // must convert column indices to global
      // Reserve space to store the destination matrix's local column indices.
      const size_t maxNumEnt = src->graph_.getNodeMaxNumRowEntries ();
      Teuchos::Array<LO> lclDstCols (maxNumEnt);

      // Copy local rows that are the "same" in both source and target.
      for (LO localRow = 0; localRow < static_cast<LO> (numSameIDs); ++localRow) {
        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        // If this call fails, that means the mesh row local index is
        // invalid.  That means the Import or Export is invalid somehow.
        LO err = 0;
        try {
          err = src->getLocalRowView (localRow, lclSrcCols, vals, numEntries);
        } catch (std::exception& e) {
          if (debug) {
            std::ostringstream os;
            const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
            os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
               << localRow << ", src->getLocalRowView() threw an exception: "
               << e.what ();
            std::cerr << os.str ();
          }
          throw e;
        }

        if (err != 0) {
          lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
          invalidSrcCopyRows.insert (localRow);
#endif // HAVE_TPETRA_DEBUG
        }
        else {
          if (static_cast<size_t> (numEntries) > static_cast<size_t> (lclDstCols.size ())) {
            lclErr = true;
            if (debug) {
              std::ostringstream os;
              const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
              os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
                 << localRow << ", numEntries = " << numEntries << " > maxNumEnt = "
                 << maxNumEnt << endl;
              std::cerr << os.str ();
            }
          }
          else {
            // Convert the source matrix's local column indices to the
            // destination matrix's local column indices.
            Teuchos::ArrayView<LO> lclDstColsView = lclDstCols.view (0, numEntries);
            for (LO j = 0; j < numEntries; ++j) {
              lclDstColsView[j] = dstColMap.getLocalElement (srcColMap.getGlobalElement (lclSrcCols[j]));
              if (lclDstColsView[j] == Teuchos::OrdinalTraits<LO>::invalid ()) {
                lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                invalidDstCopyCols.insert (lclDstColsView[j]);
#endif // HAVE_TPETRA_DEBUG
              }
            }
            try {
              err = this->replaceLocalValues (localRow, lclDstColsView.getRawPtr (), vals, numEntries);
            } catch (std::exception& e) {
              if (debug) {
                std::ostringstream os;
                const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
                os << "Proc " << myRank << ": copyAndPermute: At \"same\" localRow "
                   << localRow << ", this->replaceLocalValues() threw an exception: "
                   << e.what ();
                std::cerr << os.str ();
              }
              throw e;
            }
            if (err != numEntries) {
              lclErr = true;
              if (debug) {
                std::ostringstream os;
                const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
                os << "Proc " << myRank << ": copyAndPermute: At \"same\" "
                  "localRow " << localRow << ", this->replaceLocalValues "
                  "returned " << err << " instead of numEntries = "
                   << numEntries << endl;
                std::cerr << os.str ();
              }
            }
          }
        }
      }

      // Copy the "permute" local rows.
      for (size_t k = 0; k < numPermute; ++k) {
        const LO srcLclRow = static_cast<LO> (permuteFromLIDs[k]);
        const LO dstLclRow = static_cast<LO> (permuteToLIDs[k]);

        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        LO err = 0;
        try {
          err = src->getLocalRowView (srcLclRow, lclSrcCols, vals, numEntries);
        } catch (std::exception& e) {
          if (debug) {
            std::ostringstream os;
            const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
            os << "Proc " << myRank << ": copyAndPermute: At \"permute\" "
              "srcLclRow " << srcLclRow << " and dstLclRow " << dstLclRow
               << ", src->getLocalRowView() threw an exception: " << e.what ();
            std::cerr << os.str ();
          }
          throw e;
        }

        if (err != 0) {
          lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
          invalidPermuteFromRows.insert (srcLclRow);
#endif // HAVE_TPETRA_DEBUG
        }
        else {
          if (static_cast<size_t> (numEntries) > static_cast<size_t> (lclDstCols.size ())) {
            lclErr = true;
          }
          else {
            // Convert the source matrix's local column indices to the
            // destination matrix's local column indices.
            Teuchos::ArrayView<LO> lclDstColsView = lclDstCols.view (0, numEntries);
            for (LO j = 0; j < numEntries; ++j) {
              lclDstColsView[j] = dstColMap.getLocalElement (srcColMap.getGlobalElement (lclSrcCols[j]));
              if (lclDstColsView[j] == Teuchos::OrdinalTraits<LO>::invalid ()) {
                lclErr = true;
#ifdef HAVE_TPETRA_DEBUG
                invalidDstPermuteCols.insert (lclDstColsView[j]);
#endif // HAVE_TPETRA_DEBUG
              }
            }
            err = this->replaceLocalValues (dstLclRow, lclDstColsView.getRawPtr (), vals, numEntries);
            if (err != numEntries) {
              lclErr = true;
            }
          }
        }
      }
    }

    if (lclErr) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
#ifdef HAVE_TPETRA_DEBUG
      err << "copyAndPermute: The graph structure of the source of the "
        "Import or Export must be a subset of the graph structure of the "
        "target.  ";
      err << "invalidSrcCopyRows = [";
      for (typename std::set<LO>::const_iterator it = invalidSrcCopyRows.begin ();
           it != invalidSrcCopyRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidSrcCopyRows.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstCopyRows = [";
      for (typename std::set<LO>::const_iterator it = invalidDstCopyRows.begin ();
           it != invalidDstCopyRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstCopyRows.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstCopyCols = [";
      for (typename std::set<LO>::const_iterator it = invalidDstCopyCols.begin ();
           it != invalidDstCopyCols.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstCopyCols.end ()) {
          err << ",";
        }
      }
      err << "], invalidDstPermuteCols = [";
      for (typename std::set<LO>::const_iterator it = invalidDstPermuteCols.begin ();
           it != invalidDstPermuteCols.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidDstPermuteCols.end ()) {
          err << ",";
        }
      }
      err << "], invalidPermuteFromRows = [";
      for (typename std::set<LO>::const_iterator it = invalidPermuteFromRows.begin ();
           it != invalidPermuteFromRows.end (); ++it) {
        err << *it;
        typename std::set<LO>::const_iterator itp1 = it;
        itp1++;
        if (itp1 != invalidPermuteFromRows.end ()) {
          err << ",";
        }
      }
      err << "]" << std::endl;
#else
      err << "copyAndPermute: The graph structure of the source of the "
        "Import or Export must be a subset of the graph structure of the "
        "target." << std::endl;
#endif // HAVE_TPETRA_DEBUG
    }

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      const bool lclSuccess = ! (* (this->localError_));
      os << "*** Proc " << myRank << ": copyAndPermute "
         << (lclSuccess ? "succeeded" : "FAILED");
      if (lclSuccess) {
        os << endl;
      } else {
        os << ": error messages: " << this->errorMessages (); // comes w/ endl
      }
      std::cerr << os.str ();
    }
  }

  namespace { // (anonymous)

    /// \brief Return the (maximum) number of bytes required to pack a
    ///   block row's entries.
    ///
    /// \param numEnt [in] Number of block entries in the row.
    ///
    /// \param numBytesPerValue [in] Maximum number of bytes per
    ///   scalar (not block!) entry (value) of the row.
    ///
    /// \param blkSize [in] Block size of the block sparse matrix.
    ///
    /// If \c Scalar (the type of entries in the matrix) is a plain
    /// old data (POD) type like \c float or \c double, or a struct of
    /// POD (like <tt>std::complex<double></tt>), then the second
    /// argument is just <tt>sizeof(Scalar)</tt>.  If a \c Scalar
    /// instance has a size determined at run time (e.g., when calling
    /// its constructor), then the second argument is the result of
    /// <tt>PackTraits<Scalar>::packValueCount</tt>, called on a
    /// <tt>Scalar</tt> value with the correct run-time size.
    template<class LO, class GO, class D>
    size_t
    packRowCount (const size_t numEnt,
                  const size_t numBytesPerValue,
                  const size_t blkSize)
    {
      using Tpetra::Details::PackTraits;

      if (numEnt == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return 0;
      }
      else {
        // We store the number of entries as a local index (LO).
        LO numEntLO = 0; // packValueCount wants this.
        GO gid;
        const size_t numEntLen = PackTraits<LO, D>::packValueCount (numEntLO);
        const size_t gidsLen = numEnt * PackTraits<GO, D>::packValueCount (gid);
        const size_t valsLen = numEnt * numBytesPerValue * blkSize * blkSize;
        return numEntLen + gidsLen + valsLen;
      }
    }

    /// \brief Unpack and return the number of (block) entries in the
    ///   packed row.
    ///
    /// \param imports [in] All the packed data.
    /// \param offset [in] Index of \c imports at which the row starts.
    /// \param numBytes [in] Number of bytes in the packed row.
    /// \param numBytesPerValue [in] Maximum number of bytes per
    ///   scalar (not block!) entry (value) of the row.
    ///
    /// \return Number of (block) entries in the packed row.
    template<class ST, class LO, class GO, class D>
    size_t
    unpackRowCount (const typename Tpetra::Details::PackTraits<LO, D>::input_buffer_type& imports,
                    const size_t offset,
                    const size_t numBytes,
                    const size_t numBytesPerValue)
    {
      using Kokkos::subview;
      using Tpetra::Details::PackTraits;
      typedef typename PackTraits<LO, D>::input_buffer_type input_buffer_type;
      typedef typename input_buffer_type::size_type size_type;

      if (numBytes == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return static_cast<size_t> (0);
      }
      else {
        LO numEntLO = 0;
        const size_t theNumBytes = PackTraits<LO, D>::packValueCount (numEntLO);
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(
          theNumBytes > numBytes, std::logic_error, "unpackRowCount: "
          "theNumBytes = " << theNumBytes << " < numBytes = " << numBytes
          << ".");
#endif // HAVE_TPETRA_DEBUG
        const std::pair<size_type, size_type> rng (offset, offset + theNumBytes);
        input_buffer_type inBuf = subview (imports, rng); // imports (offset, theNumBytes);
        const size_t actualNumBytes = PackTraits<LO, D>::unpackValue (numEntLO, inBuf);
#ifdef HAVE_TPETRA_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(
          theNumBytes > numBytes, std::logic_error, "unpackRowCount: "
          "actualNumBytes = " << actualNumBytes << " < numBytes = " << numBytes
          << ".");
#else
        (void) actualNumBytes;
#endif // HAVE_TPETRA_DEBUG
        return static_cast<size_t> (numEntLO);
      }
    }

    /// \brief Pack the block row (stored in the input arrays).
    ///
    /// \return The number of bytes packed.
    template<class ST, class LO, class GO, class D>
    size_t
    packRow (const typename Tpetra::Details::PackTraits<LO, D>::output_buffer_type& exports,
             const size_t offset,
             const size_t numEnt,
             const typename Tpetra::Details::PackTraits<GO, D>::input_array_type& gidsIn,
             const typename Tpetra::Details::PackTraits<ST, D>::input_array_type& valsIn,
             const size_t numBytesPerValue,
             const size_t blockSize)
    {
      using Kokkos::subview;
      using Tpetra::Details::PackTraits;
      // NOTE (mfh 02 Feb 2015) This assumes that output_buffer_type is
      // the same, no matter what type we're packing.  It's a reasonable
      // assumption, given that we go through the trouble of PackTraits
      // just so that we can pack data of different types in the same
      // buffer.
      typedef typename PackTraits<LO, D>::output_buffer_type output_buffer_type;
      typedef typename output_buffer_type::size_type size_type;
      typedef typename std::pair<size_type, size_type> pair_type;

      if (numEnt == 0) {
        // Empty rows always take zero bytes, to ensure sparsity.
        return 0;
      }
      const size_t numScalarEnt = numEnt * blockSize * blockSize;

      const GO gid = 0; // packValueCount wants this
      const LO numEntLO = static_cast<size_t> (numEnt);

      const size_t numEntBeg = offset;
      const size_t numEntLen = PackTraits<LO, D>::packValueCount (numEntLO);
      const size_t gidsBeg = numEntBeg + numEntLen;
      const size_t gidsLen = numEnt * PackTraits<GO, D>::packValueCount (gid);
      const size_t valsBeg = gidsBeg + gidsLen;
      const size_t valsLen = numScalarEnt * numBytesPerValue;

      output_buffer_type numEntOut =
        subview (exports, pair_type (numEntBeg, numEntBeg + numEntLen));
      output_buffer_type gidsOut =
        subview (exports, pair_type (gidsBeg, gidsBeg + gidsLen));
      output_buffer_type valsOut =
        subview (exports, pair_type (valsBeg, valsBeg + valsLen));

      size_t numBytesOut = 0;
      numBytesOut += PackTraits<LO, D>::packValue (numEntOut, numEntLO);
      numBytesOut += PackTraits<GO, D>::packArray (gidsOut, gidsIn, numEnt);
      numBytesOut += PackTraits<ST, D>::packArray (valsOut, valsIn, numScalarEnt);

      const size_t expectedNumBytes = numEntLen + gidsLen + valsLen;
      TEUCHOS_TEST_FOR_EXCEPTION(
        numBytesOut != expectedNumBytes, std::logic_error, "unpackRow: "
        "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
        << expectedNumBytes << ".");
      return numBytesOut;
    }

    // Return the number of bytes actually read / used.
    template<class ST, class LO, class GO, class D>
    size_t
    unpackRow (const typename Tpetra::Details::PackTraits<GO, D>::output_array_type& gidsOut,
               const typename Tpetra::Details::PackTraits<ST, D>::output_array_type& valsOut,
               const typename Tpetra::Details::PackTraits<int, D>::input_buffer_type& imports,
               const size_t offset,
               const size_t numBytes,
               const size_t numEnt,
               const size_t numBytesPerValue,
               const size_t blockSize)
    {
      using Kokkos::subview;
      using Tpetra::Details::PackTraits;
      // NOTE (mfh 02 Feb 2015) This assumes that input_buffer_type is
      // the same, no matter what type we're unpacking.  It's a
      // reasonable assumption, given that we go through the trouble of
      // PackTraits just so that we can pack data of different types in
      // the same buffer.
      typedef typename PackTraits<LO, D>::input_buffer_type input_buffer_type;
      typedef typename input_buffer_type::size_type size_type;
      typedef typename std::pair<size_type, size_type> pair_type;

      if (numBytes == 0) {
        // Rows with zero bytes always have zero entries.
        return 0;
      }
      const size_t numScalarEnt = numEnt * blockSize * blockSize;
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (imports.dimension_0 ()) <= offset,
        std::logic_error, "unpackRow: imports.dimension_0() = "
        << imports.dimension_0 () << " <= offset = " << offset << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (imports.dimension_0 ()) < offset + numBytes,
        std::logic_error, "unpackRow: imports.dimension_0() = "
        << imports.dimension_0 () << " < offset + numBytes = "
        << (offset + numBytes) << ".");

      const GO gid = 0; // packValueCount wants this
      const LO lid = 0; // packValueCount wants this

      const size_t numEntBeg = offset;
      const size_t numEntLen = PackTraits<LO, D>::packValueCount (lid);
      const size_t gidsBeg = numEntBeg + numEntLen;
      const size_t gidsLen = numEnt * PackTraits<GO, D>::packValueCount (gid);
      const size_t valsBeg = gidsBeg + gidsLen;
      const size_t valsLen = numScalarEnt * numBytesPerValue;

      input_buffer_type numEntIn =
        subview (imports, pair_type (numEntBeg, numEntBeg + numEntLen));
      input_buffer_type gidsIn =
        subview (imports, pair_type (gidsBeg, gidsBeg + gidsLen));
      input_buffer_type valsIn =
        subview (imports, pair_type (valsBeg, valsBeg + valsLen));

      size_t numBytesOut = 0;
      LO numEntOut;
      numBytesOut += PackTraits<LO, D>::unpackValue (numEntOut, numEntIn);
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<size_t> (numEntOut) != numEnt, std::logic_error,
        "unpackRow: Expected number of entries " << numEnt
        << " != actual number of entries " << numEntOut << ".");

      numBytesOut += PackTraits<GO, D>::unpackArray (gidsOut, gidsIn, numEnt);
      numBytesOut += PackTraits<ST, D>::unpackArray (valsOut, valsIn, numScalarEnt);

      TEUCHOS_TEST_FOR_EXCEPTION(
        numBytesOut != numBytes, std::logic_error, "unpackRow: numBytesOut = "
        << numBytesOut << " != numBytes = " << numBytes << ".");
      const size_t expectedNumBytes = numEntLen + gidsLen + valsLen;
      TEUCHOS_TEST_FOR_EXCEPTION(
        numBytesOut != expectedNumBytes, std::logic_error, "unpackRow: "
        "numBytesOut = " << numBytesOut << " != expectedNumBytes = "
        << expectedNumBytes << ".");
      return numBytesOut;
    }
  } // namespace (anonymous)

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  packAndPrepare (const Tpetra::SrcDistObject& source,
                  const Teuchos::ArrayView<const LO>& exportLIDs,
                  Teuchos::Array<packet_type>& exports,
                  const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                  size_t& constantNumPackets,
                  Tpetra::Distributor& /* distor */)
  {
    using std::endl;
    using Tpetra::Details::PackTraits;
    using Kokkos::MemoryUnmanaged;
    using Kokkos::subview;
    using Kokkos::View;
    typedef typename Tpetra::MultiVector<Scalar, LO, GO, Node>::impl_scalar_type ST;
    typedef typename Node::execution_space execution_space;
    typedef typename View<int*, execution_space>::HostMirror::execution_space HES;
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
    const bool debug = false;

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": packAndPrepare: exportLIDs.size() = "
         << exportLIDs.size () << ", numPacketsPerLID.size() = "
         << numPacketsPerLID.size () << endl;
      std::cerr << os.str ();
    }

    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << "packAndPrepare: The target object of the Import or Export is "
        "already in an error state." << endl;
      return;
    }

    const this_type* src = dynamic_cast<const this_type* > (&source);
    // Should have checked for these cases in checkSizes().
    if (src == NULL) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << "packAndPrepare: The source (input) object of the Import or "
        "Export is either not a BlockCrsMatrix, or does not have the right "
        "template parameters.  checkSizes() should have caught this.  "
        "Please report this bug to the Tpetra developers." << endl;
      return;
    }

    const crs_graph_type& srcGraph = src->graph_;
    const size_t blockSize = static_cast<size_t> (src->getBlockSize ());
    const size_type numExportLIDs = exportLIDs.size ();

    if (numExportLIDs != numPacketsPerLID.size ()) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << "packAndPrepare: exportLIDs.size() = " << numExportLIDs
          << " != numPacketsPerLID.size() = " << numPacketsPerLID.size ()
          << "." << endl;
      return;
    }

    // Graphs and matrices are allowed to have a variable number of
    // entries per row.  We could test whether all rows have the same
    // number of entries, but DistObject can only use this
    // optimization if all rows on _all_ processes have the same
    // number of entries.  Rather than do the all-reduce necessary to
    // test for this unlikely case, we tell DistObject (by setting
    // constantNumPackets to zero) to assume that different rows may
    // have different numbers of entries.
    constantNumPackets = 0;

    // Compute the number of bytes ("packets") per row to pack.  While
    // we're at it, compute the total # of block entries to send, and
    // the max # of block entries in any of the rows we're sending.
    size_t totalNumBytes = 0;
    size_t totalNumEntries = 0;
    size_t maxRowLength = 0;
    for (size_type i = 0; i < numExportLIDs; ++i) {
      const LO lclRow = exportLIDs[i];
      size_t numEnt = srcGraph.getNumEntriesInLocalRow (lclRow);
      // If any given LIDs are invalid, the above might return either
      // zero or the invalid size_t value.  If the former, we have no
      // way to tell, but that's OK; it just means the calling process
      // won't pack anything (it has nothing to pack anyway).  If the
      // latter, we replace it with zero (that row is not owned by the
      // calling process, so it has no entries to pack).
      if (numEnt == Teuchos::OrdinalTraits<size_t>::invalid ()) {
        numEnt = 0;
      }
      const size_t numScalarEnt = numEnt * blockSize * blockSize;

      // The 'if' branch implicitly assumes that packRowCount() returns
      // zero if numEnt == 0.
      size_t numBytesPerValue = 0;
      if (numEnt > 0) {
        // Get a locally indexed view of the current row's data.  If
        // the current row has > 0 entries, we need an entry in order
        // to figure out the byte count of the packed row.  (We really
        // only need it if ST's size is determined at run time.)
        Scalar* valsRaw = NULL;
        const LO* lidsRaw = NULL;
        LO actualNumEnt = 0;
        const LO errCode =
          src->getLocalRowView (lclRow, lidsRaw, valsRaw, actualNumEnt);

        if (numEnt != static_cast<size_t> (actualNumEnt)) {
          std::ostream& err = this->markLocalErrorAndGetStream ();
          err << "packAndPrepare: Local row " << i << " claims to have " <<
            numEnt << "entry/ies, but the View returned by getLocalRowView() "
            "has " << actualNumEnt << " entry/ies.  This should never happen.  "
            "Please report this bug to the Tpetra developers." << endl;
          return;
        }
        if (errCode == Teuchos::OrdinalTraits<LO>::invalid ()) {
          std::ostream& err = this->markLocalErrorAndGetStream ();
          err << "packAndPrepare: Local row " << i << " is not in the row Map "
            "of the source object on the calling process." << endl;
          return;
        }

        const ST* valsRawST = const_cast<const ST*> (reinterpret_cast<ST*> (valsRaw));
        View<const ST*, HES, MemoryUnmanaged> vals (valsRawST, numScalarEnt);

        // NOTE (mfh 07 Feb 2015) Since we're using the host memory
        // space here for now, this doesn't assume UVM.  That may change
        // in the future, if we ever start packing on the device.
        numBytesPerValue = PackTraits<ST, HES>::packValueCount (vals(0));
      }

      const size_t numBytes =
        packRowCount<LO, GO, HES> (numEnt, numBytesPerValue, blockSize);
      numPacketsPerLID[i] = numBytes;
      totalNumBytes += numBytes;
      totalNumEntries += numEnt;
      maxRowLength = std::max (maxRowLength, numEnt);
    }

    if (debug) {
      const int myRank = graph_.getComm ()->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": packAndPrepare: totalNumBytes = "
         << totalNumBytes << endl;
      std::cerr << os.str ();
    }

    // We use a "struct of arrays" approach to packing each row's
    // entries.  All the column indices (as global indices) go first,
    // then all their owning process ranks, and then the values.
    exports.resize (totalNumBytes);
    if (totalNumEntries > 0) {
      View<char*, HES, MemoryUnmanaged> exportsK (exports.getRawPtr (),
                                                  totalNumBytes);

      // Current position (in bytes) in the 'exports' output array.
      size_t offset = 0;

      // For each block row of the matrix owned by the calling
      // process, pack that block row's column indices and values into
      // the exports array.

      // Source matrix's column Map.  We verified in checkSizes() that
      // the column Map exists (is not null).
      const map_type& srcColMap = * (srcGraph.getColMap ());

      // Temporary buffer for global column indices.
      View<GO*, HES> gblColInds;
      {
        GO gid = 0;
        gblColInds = PackTraits<GO, HES>::allocateArray (gid, maxRowLength, "gids");
      }

      // Pack the data for each row to send, into the 'exports' buffer.
      for (size_type i = 0; i < numExportLIDs; ++i) {
        const LO lclRowInd = exportLIDs[i];
        const LO* lclColIndsRaw;
        Scalar* valsRaw;
        LO numEntLO;
        // It's OK to ignore the return value, since if the calling
        // process doesn't own that local row, then the number of
        // entries in that row on the calling process is zero.
        (void) src->getLocalRowView (lclRowInd, lclColIndsRaw, valsRaw, numEntLO);
        const size_t numEnt = static_cast<size_t> (numEntLO);
        const size_t numScalarEnt = numEnt * blockSize * blockSize;
        View<const LO*, HES, MemoryUnmanaged> lclColInds (lclColIndsRaw, numEnt);
        const ST* valsRawST = const_cast<const ST*> (reinterpret_cast<ST*> (valsRaw));
        View<const ST*, HES, MemoryUnmanaged> vals (valsRawST, numScalarEnt);

        // NOTE (mfh 07 Feb 2015) Since we're using the host memory
        // space here for now, this doesn't assume UVM.  That may
        // change in the future, if we ever start packing on device.
        const size_t numBytesPerValue = numEnt == 0 ?
          static_cast<size_t> (0) :
          PackTraits<ST, HES>::packValueCount (vals(0));

        // Convert column indices from local to global.
        for (size_t j = 0; j < numEnt; ++j) {
          gblColInds(j) = srcColMap.getGlobalElement (lclColInds(j));
        }

        // Copy the row's data into the current spot in the exports array.
        const size_t numBytes =
          packRow<ST, LO, GO, HES> (exportsK, offset, numEnt, gblColInds,
                                    vals, numBytesPerValue, blockSize);
        // Keep track of how many bytes we packed.
        offset += numBytes;
      } // for each LID (of a row) to send

      if (offset != totalNumBytes) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << "packAndPreapre: At end of method, the final offset (in bytes) "
            << offset << " does not equal the total number of bytes packed "
            << totalNumBytes << ".  "
            << "Please report this bug to the Tpetra developers." << endl;
        return;
      }
    } // if totalNumEntries > 0

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      const bool lclSuccess = ! (* (this->localError_));
      os << "*** Proc " << myRank << ": packAndPrepare "
         << (lclSuccess ? "succeeded" : "FAILED")
         << " (totalNumEntries = " << totalNumEntries << ") ***" << endl;
      std::cerr << os.str ();
    }
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  unpackAndCombine (const Teuchos::ArrayView<const LO>& importLIDs,
                    const Teuchos::ArrayView<const packet_type>& imports,
                    const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                    size_t /* constantNumPackets */, // not worthwhile to use this
                    Tpetra::Distributor& /* distor */,
                    Tpetra::CombineMode CM)
  {
    using std::endl;
    using Tpetra::Details::PackTraits;
    using Kokkos::MemoryUnmanaged;
    using Kokkos::subview;
    using Kokkos::View;
    typedef typename Tpetra::MultiVector<Scalar, LO, GO, Node>::impl_scalar_type ST;
    typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
    typedef typename Node::execution_space execution_space;
    typedef typename View<int*, execution_space>::HostMirror::execution_space HES;
    typedef std::pair<typename View<int*, HES>::size_type,
      typename View<int*, HES>::size_type> pair_type;
    typedef View<GO*, HES, MemoryUnmanaged> gids_out_type;
    typedef View<LO*, HES, MemoryUnmanaged> lids_out_type;
    typedef View<ST*, HES, MemoryUnmanaged> vals_out_type;
    typedef typename PackTraits<GO, HES>::input_buffer_type input_buffer_type;
    const char prefix[] = "Tpetra::Experimental::BlockCrsMatrix::unpackAndCombine: ";
    const bool debug = false;

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": unpackAndCombine" << endl;
      std::cerr << os.str ();
    }

    // It should not cause deadlock to return on error in this method,
    // since this method does not communicate.

    if (* (this->localError_)) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "The target object of the Import or Export is "
        "already in an error state." << endl;
      return;
    }
    if (importLIDs.size () != numPacketsPerLID.size ()) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "importLIDs.size() = " << importLIDs.size () << " != "
        "numPacketsPerLID.size() = " << numPacketsPerLID.size () << "." << endl;
      return;
    }
    if (CM != ADD && CM != INSERT && CM != REPLACE && CM != ABSMAX && CM != ZERO) {
      std::ostream& err = this->markLocalErrorAndGetStream ();
      err << prefix << "Invalid CombineMode value " << CM << ".  Valid "
          << "values include ADD, INSERT, REPLACE, ABSMAX, and ZERO."
          << endl;
      return;
    }

    // Target matrix's column Map.  Use to convert the global column
    // indices in the receive buffer to local indices.  We verified in
    // checkSizes() that the column Map exists (is not null).
    const map_type& tgtColMap = * (this->graph_.getColMap ());

    const size_type numImportLIDs = importLIDs.size ();
    if (CM == ZERO || numImportLIDs == 0) {
      if (debug) {
        std::ostringstream os;
        const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
        os << "Proc " << myRank << ": unpackAndCombine: Nothing to do" << endl;
        std::cerr << os.str ();
      }
      return; // nothing to do; no need to combine entries
    }

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      os << "Proc " << myRank << ": unpackAndCombine: Getting ready" << endl;
      std::cerr << os.str ();
    }

    input_buffer_type importsK (imports.getRawPtr (), imports.size ());
    const size_t blockSize = this->getBlockSize ();
    const size_t maxRowNumEnt = graph_.getNodeMaxNumRowEntries ();
    const size_t maxRowNumScalarEnt = maxRowNumEnt * blockSize * blockSize;

    size_t numBytesPerValue;
    {
      // FIXME (mfh 17 Feb 2015) What do I do about Scalar types with
      // run-time size?  We already assume that all entries in both
      // the source and target matrices have the same size.  If the
      // calling process owns at least one entry in either matrix, we
      // can use that entry to set the size.  However, it is possible
      // that the calling process owns no entries.  In that case,
      // we're in trouble.  One way to fix this would be for each
      // row's data to contain the run-time size.  This is only
      // necessary if the size is not a compile-time constant.
      Scalar val;
      numBytesPerValue = PackTraits<ST, HES>::packValueCount (val);
    }

    // Temporary space to cache incoming global column indices and
    // values.  Column indices come in as global indices, in case the
    // source object's column Map differs from the target object's
    // (this's) column Map.
    View<GO*, HES> gblColInds;
    View<LO*, HES> lclColInds;
    View<ST*, HES> vals;
    {
      GO gid = 0;
      LO lid = 0;
      // FIXME (mfh 17 Feb 2015) What do I do about Scalar types with
      // run-time size?  We already assume that all entries in both
      // the source and target matrices have the same size.  If the
      // calling process owns at least one entry in either matrix, we
      // can use that entry to set the size.  However, it is possible
      // that the calling process owns no entries.  In that case,
      // we're in trouble.  One way to fix this would be for each
      // row's data to contain the run-time size.  This is only
      // necessary if the size is not a compile-time constant.
      Scalar val;
      gblColInds = PackTraits<GO, HES>::allocateArray (gid, maxRowNumEnt, "gids");
      lclColInds = PackTraits<LO, HES>::allocateArray (lid, maxRowNumEnt, "lids");
      vals = PackTraits<ST, HES>::allocateArray (val, maxRowNumScalarEnt, "vals");
    }

    size_t offset = 0;
    for (size_type i = 0; i < numImportLIDs; ++i) {
      const size_t numBytes = numPacketsPerLID[i];
      if (numBytes == 0) {
        // Empty buffer for that row means that the row is empty.
        continue;
      }
      const size_t numEnt =
        unpackRowCount<ST, LO, GO, HES> (importsK, offset, numBytes, numBytesPerValue);
      if (numEnt > maxRowNumEnt) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix << "At i = " << i << ", numEnt = " << numEnt
            << " > maxRowNumEnt = " << maxRowNumEnt << endl;
        continue;
      }

      const size_t numScalarEnt = numEnt * blockSize * blockSize;
      const LO lclRow = importLIDs[i];

      gids_out_type gidsOut = subview (gblColInds, pair_type (0, numEnt));
      vals_out_type valsOut = subview (vals, pair_type (0, numScalarEnt));

      const size_t numBytesOut =
        unpackRow<ST, LO, GO, HES> (gidsOut, valsOut, importsK, offset, numBytes,
                                    numEnt, numBytesPerValue, blockSize);
      if (numBytes != numBytesOut) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix << "At i = " << i << ", numBytes = " << numBytes
            << " != numBytesOut = " << numBytesOut << ".";
        continue;
      }

      // Convert incoming global indices to local indices.
      lids_out_type lidsOut = subview (lclColInds, pair_type (0, numEnt));
      for (size_t k = 0; k < numEnt; ++k) {
        lidsOut(k) = tgtColMap.getLocalElement (gidsOut(k));
#ifdef HAVE_TPETRA_DEBUG
        if (lidsOut(k) == Teuchos::OrdinalTraits<LO>::invalid ()) {
          std::ostream& err = this->markLocalErrorAndGetStream ();
          err << prefix << "At i = " << i << ", GID " << gidsOut(k)
              << " is not owned by the calling process.";
          continue;
        }
#endif // HAVE_TPETRA_DEBUG
      }

      // Combine the incoming data with the matrix's current data.
      LO numCombd = 0;
      const LO* const lidsRaw = const_cast<const LO*> (lidsOut.ptr_on_device ());
      const Scalar* const valsRaw =
        reinterpret_cast<const Scalar*> (const_cast<const ST*> (valsOut.ptr_on_device ()));
      if (CM == ADD) {
        numCombd = this->sumIntoLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
      } else if (CM == INSERT || CM == REPLACE) {
        numCombd = this->replaceLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
      } else if (CM == ABSMAX) {
        numCombd = this->absMaxLocalValues (lclRow, lidsRaw, valsRaw, numEnt);
      }
#ifdef HAVE_TPETRA_DEBUG
      if (static_cast<LO> (numEnt) != numCombd) {
        std::ostream& err = this->markLocalErrorAndGetStream ();
        err << prefix << "At i = " << i << ", numEnt = " << numEnt
            << " != numCombd = " << numCombd << ".";
        continue;
      }
#else
      (void) numCombd; // ignore, just for now
#endif // HAVE_TPETRA_DEBUG

      // Don't update offset until current LID has succeeded.
      offset += numBytes;
    } // for each import LID i

    if (debug) {
      std::ostringstream os;
      const int myRank = this->graph_.getRowMap ()->getComm ()->getRank ();
      const bool lclSuccess = ! (* (this->localError_));
      os << "*** Proc " << myRank << ": unpackAndCombine "
         << (lclSuccess ? "succeeded" : "FAILED")
         << " ***" << endl;
      std::cerr << os.str ();
    }
  }


  template<class Scalar, class LO, class GO, class Node>
  std::string
  BlockCrsMatrix<Scalar, LO, GO, Node>::description () const
  {
    using Teuchos::TypeNameTraits;
    std::ostringstream os;
    os << "\"Tpetra::BlockCrsMatrix\": { "
       << "Template parameters: { "
       << "Scalar: " << TypeNameTraits<Scalar>::name ()
       << "LO: " << TypeNameTraits<LO>::name ()
       << "GO: " << TypeNameTraits<GO>::name ()
       << "Node: " << TypeNameTraits<Node>::name ()
       << " }"
       << ", Label: \"" << this->getObjectLabel () << "\""
       << ", Global dimensions: ["
       << graph_.getDomainMap ()->getGlobalNumElements () << ", "
       << graph_.getRangeMap ()->getGlobalNumElements () << "]"
       << ", Block size: " << getBlockSize ()
       << " }";
    return os.str ();
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::CommRequest;
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::outArg;
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    // using Teuchos::VERB_MEDIUM;
    // using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    using Teuchos::RCP;
    using Teuchos::wait;
    using std::endl;

    // Set default verbosity if applicable.
    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    if (vl == VERB_NONE) {
      return; // print nothing
    }

    // describe() always starts with a tab before it prints anything.
    Teuchos::OSTab tab0 (out);

    out << "\"Tpetra::BlockCrsMatrix\":" << endl;
    Teuchos::OSTab tab1 (out);

    out << "Template parameters:" << endl;
    {
      Teuchos::OSTab tab2 (out);
      out << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
          << "LO: " << TypeNameTraits<LO>::name () << endl
          << "GO: " << TypeNameTraits<GO>::name () << endl
          << "Node: " << TypeNameTraits<Node>::name () << endl;
    }
    out << "Label: \"" << this->getObjectLabel () << "\"" << endl
        << "Global dimensions: ["
        << graph_.getDomainMap ()->getGlobalNumElements () << ", "
        << graph_.getRangeMap ()->getGlobalNumElements () << "]" << endl;

    const LO blockSize = getBlockSize ();
    out << "Block size: " << blockSize << endl;

    if (vl >= VERB_EXTREME) {
      const Teuchos::Comm<int>& comm = * (graph_.getMap ()->getComm ());
      const int myRank = comm.getRank ();
      const int numProcs = comm.getSize ();

      // Print the calling process' data to the given output stream.
      RCP<std::ostringstream> lclOutStrPtr (new std::ostringstream ());
      RCP<FancyOStream> osPtr = getFancyOStream (lclOutStrPtr);
      FancyOStream& os = *osPtr;
      os << "Process " << myRank << ":" << endl;
      Teuchos::OSTab tab2 (os);

      const map_type& meshRowMap = * (graph_.getRowMap ());
      const map_type& meshColMap = * (graph_.getColMap ());
      for (LO meshLclRow = meshRowMap.getMinLocalIndex ();
           meshLclRow <= meshRowMap.getMaxLocalIndex ();
           ++meshLclRow) {
        const GO meshGblRow = meshRowMap.getGlobalElement (meshLclRow);
        os << "Row " << meshGblRow << ": {";

        const LO* lclColInds = NULL;
        Scalar* vals = NULL;
        LO numInds = 0;
        this->getLocalRowView (meshLclRow, lclColInds, vals, numInds);

        for (LO k = 0; k < numInds; ++k) {
          const GO gblCol = meshColMap.getGlobalElement (lclColInds[k]);

          os << "Col " << gblCol << ": [";
          for (LO i = 0; i < blockSize; ++i) {
            for (LO j = 0; j < blockSize; ++j) {
              os << vals[blockSize*blockSize*k + i*blockSize + j];
              if (j + 1 < blockSize) {
                os << ", ";
              }
            }
            if (i + 1 < blockSize) {
              os << "; ";
            }
          }
          os << "]";
          if (k + 1 < numInds) {
            os << ", ";
          }
        }
        os << "}" << endl;
      }

      // Print data on Process 0.  This will automatically respect the
      // current indentation level.
      if (myRank == 0) {
        out << lclOutStrPtr->str ();
        lclOutStrPtr = Teuchos::null; // clear it to save space
      }

      const int sizeTag = 1337;
      const int dataTag = 1338;

      ArrayRCP<char> recvDataBuf; // only used on Process 0

      // Send string sizes and data from each process in turn to
      // Process 0, and print on that process.
      for (int p = 1; p < numProcs; ++p) {
        if (myRank == 0) {
          // Receive the incoming string's length.
          ArrayRCP<size_t> recvSize (1);
          recvSize[0] = 0;
          RCP<CommRequest<int> > recvSizeReq =
            ireceive<int, size_t> (recvSize, p, sizeTag, comm);
          wait<int> (comm, outArg (recvSizeReq));
          const size_t numCharsToRecv = recvSize[0];

          // Allocate space for the string to receive.  Reuse receive
          // buffer space if possible.  We can do this because in the
          // current implementation, we only have one receive in
          // flight at a time.  Leave space for the '\0' at the end,
          // in case the sender doesn't send it.
          if (static_cast<size_t>(recvDataBuf.size()) < numCharsToRecv + 1) {
            recvDataBuf.resize (numCharsToRecv + 1);
          }
          ArrayRCP<char> recvData = recvDataBuf.persistingView (0, numCharsToRecv);
          // Post the receive of the actual string data.
          RCP<CommRequest<int> > recvDataReq =
            ireceive<int, char> (recvData, p, dataTag, comm);
          wait<int> (comm, outArg (recvDataReq));

          // Print the received data.  This will respect the current
          // indentation level.  Make sure that the string is
          // null-terminated.
          recvDataBuf[numCharsToRecv] = '\0';
          out << recvDataBuf.getRawPtr ();
        }
        else if (myRank == p) { // if I am not Process 0, and my rank is p
          // This deep-copies the string at most twice, depending on
          // whether std::string reference counts internally (it
          // generally does, so this won't deep-copy at all).
          const std::string stringToSend = lclOutStrPtr->str ();
          lclOutStrPtr = Teuchos::null; // clear original to save space

          // Send the string's length to Process 0.
          const size_t numCharsToSend = stringToSend.size ();
          ArrayRCP<size_t> sendSize (1);
          sendSize[0] = numCharsToSend;
          RCP<CommRequest<int> > sendSizeReq =
            isend<int, size_t> (sendSize, 0, sizeTag, comm);
          wait<int> (comm, outArg (sendSizeReq));

          // Send the actual string to Process 0.  We know that the
          // string has length > 0, so it's save to take the address
          // of the first entry.  Make a nonowning ArrayRCP to hold
          // the string.  Process 0 will add a null termination
          // character at the end of the string, after it receives the
          // message.
          ArrayRCP<const char> sendData (&stringToSend[0], 0, numCharsToSend, false);
          RCP<CommRequest<int> > sendDataReq =
            isend<int, char> (sendData, 0, dataTag, comm);
          wait<int> (comm, outArg (sendDataReq));
        }
      } // for each process rank p other than 0
    } // extreme verbosity level (print the whole matrix)
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const Teuchos::Comm<int> >
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getComm() const
  {
    return graph_.getComm();
  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNode() const
  {
    return graph_.getNode();

  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumCols() const
  {
    return graph_.getGlobalNumCols();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeNumCols() const
  {
    return graph_.getNodeNumCols();
  }

  template<class Scalar, class LO, class GO, class Node>
  GO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getIndexBase() const
  {
    return graph_.getIndexBase();
  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumEntries() const
  {
    return graph_.getGlobalNumEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeNumEntries() const
  {
    return graph_.getNodeNumEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNumEntriesInGlobalRow (GO globalRow) const
  {
    return graph_.getNumEntriesInGlobalRow(globalRow);
  }

  template<class Scalar, class LO, class GO, class Node>
  global_size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalNumDiags() const
  {
    return getGlobalNumDiags();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNodeNumDiags() const
  {
    return getNodeNumDiags();
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalMaxNumRowEntries() const
  {
    return graph_.getGlobalMaxNumRowEntries();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  hasColMap() const
  {
    return graph_.hasColMap();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isLowerTriangular() const
  {
    return graph_.isLowerTriangular();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isUpperTriangular() const
  {
    return graph_.isUpperTriangular();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isLocallyIndexed() const
  {
    return graph_.isLocallyIndexed();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isGloballyIndexed() const
  {
    return graph_.isGloballyIndexed();
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  isFillComplete() const
  {
    return true;
  }

  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  supportsRowViews() const
  {
    return true;
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalRowCopy (GO GlobalRow,
                    const Teuchos::ArrayView<GO> &Indices,
                    const Teuchos::ArrayView<Scalar> &Values,
                    size_t &NumEntries) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getGlobalRowCopy: "
      "This class doesn't support global matrix indexing.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGlobalRowView (GO GlobalRow,
                    ArrayView<const GO> &indices,
                    ArrayView<const Scalar> &values) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getGlobalRowView: "
      "This class doesn't support global matrix indexing.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowView (LO LocalRow,
                   ArrayView<const LO> &indices,
                   ArrayView<const Scalar> &values) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getGlobalRowView: "
      "This class doesn't support global matrix indexing.");

  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalDiagCopy (Tpetra::Vector<Scalar,LO,GO,Node> &diag) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getLocalDiagCopy: "
      "not implemented.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  leftScale (const Tpetra::Vector<Scalar, LO, GO, Node>& x)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::leftScale: "
      "not implemented.");

  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  rightScale (const Tpetra::Vector<Scalar, LO, GO, Node>& x)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::rightScale: "
      "not implemented.");

  }

  template<class Scalar, class LO, class GO, class Node>
  Teuchos::RCP<const Tpetra::RowGraph<LO, GO, Node> >
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getGraph() const
  {
    return graphRCP_;
  }

  template<class Scalar, class LO, class GO, class Node>
  typename Tpetra::RowMatrix<Scalar, LO, GO, Node>::mag_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getFrobeniusNorm () const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::getFrobeniusNorm: "
      "not implemented.");
  }

} // namespace Experimental
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_INSTANT(S,LO,GO,NODE) \
  template class Experimental::BlockCrsMatrix< S, LO, GO, NODE >;

#endif // TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DEF_HPP
