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

namespace Tpetra {
namespace Experimental {

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix () :
    dist_object_type (Teuchos::rcp (new map_type ())), // nonnull, so DistObject doesn't throw
    graph_ (Teuchos::rcp (new map_type ()), 0), // FIXME (mfh 16 May 2014) no empty ctor yet
    blockSize_ (static_cast<LO> (0)),
    ptr_ (NULL),
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
    ptr_ (NULL), // to be initialized below
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

    ptr_ = graph.getNodeRowPtrs ().getRawPtr ();
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
    ptr_ (NULL), // to be initialized below
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

    ptr_ = graph.getNodeRowPtrs ().getRawPtr ();
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
         scalar_type alpha,
         scalar_type beta) const
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
        A_cur.fill (alpha);
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
          getConstLocalBlockFromInput (vals, pointOffset);
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
      "BlockCrsMatrix::getDiagonalGraph: You must call computeDiagionalGraph() "
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
    Teuchos::Array<Scalar> localMem (blockSize);
    Teuchos::Array<Scalar> localMat (blockSize*blockSize);
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
        little_block_type D_lcl = getNonConstLocalBlockFromInput (Dmat, 0);

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
          little_block_type D_lcl = getNonConstLocalBlockFromInput(Dmat, 0);

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
          getConstLocalBlockFromInput (vals, pointOffset);
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
          getConstLocalBlockFromInput (vals, pointOffset);
        A_old.update (STS::one (), A_new);
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
      vals = val_ + absBlockOffsetStart * offsetPerBlock ();
      numInds = ptr_[localRowInd + 1] - absBlockOffsetStart;
      return 0; // indicates no error
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowCopy (LO LocalRow,
                   const ArrayView<LO> &Indices,
                   const ArrayView<Scalar> &Values,
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
          getConstLocalBlockFromInput (vals, pointOffset);
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
          getConstLocalBlockFromInput (vals, pointOffset);
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
          getConstLocalBlockFromInput (vals, pointOffset);
        A_old.update (STS::one (), A_new);
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
    Teuchos::Array<Scalar> localMem (blockSize);
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
  getConstLocalBlockFromInput (const Scalar* val, const size_t pointOffset) const
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
  getNonConstLocalBlockFromInput (Scalar* val, const size_t pointOffset) const
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
      return getNonConstLocalBlockFromInput (const_cast<Scalar*> (val_),
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
    else
    {
      return little_block_type (NULL, 0, 0, 0);
    }

  }


  template<class Scalar, class LO, class GO, class Node>
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  checkSizes (const Tpetra::SrcDistObject& source)
  {
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    const this_type* src = dynamic_cast<const this_type* > (&source);

    // Clear out the current local error state.
    * (const_cast<this_type*> (this)->localError_) = false;
    *errs_ = Teuchos::null;

    // We don't allow block sizes to be inconsistent across processes,
    // but it's possible due to user error.  It costs an all-reduce to
    // check in the constructor; we want to save that.
    if (src == NULL ||
        src->getBlockSize () != this->getBlockSize () ||
        ! src->graph_.isFillComplete () ||
        ! this->graph_.isFillComplete () ||
        src->graph_.getColMap ().is_null () ||
        this->graph_.getColMap ().is_null ()) {
      * (const_cast<this_type*> (this)->localError_) = true;
      if ((*errs_).is_null ()) {
        *errs_ = Teuchos::rcp (new std::ostringstream ());
      }
    }

    if (src == NULL) {
      **errs_ << "checkSizes: The source object of the Import or Export "
        "must be a BlockCrsMatrix with the same template parameters as the "
        "target object." << std::endl;
    }
    else {
      // Use a string of ifs, not if-elseifs, because we want to know
      // all the errors.
      if (src->getBlockSize () != this->getBlockSize ()) {
        **errs_ << "checkSizes: The source and target objects of the Import or "
               << "Export must have the same block sizes.  The source's block "
               << "size = " << src->getBlockSize () << " != the target's block "
               << "size = " << this->getBlockSize () << "." << std::endl;
      }
      if (! src->graph_.isFillComplete ()) {
        **errs_ << "checkSizes: The source object of the Import or Export is "
          "not fill complete.  Both source and target objects must be fill "
          "complete." << std::endl;
      }
      if (! this->graph_.isFillComplete ()) {
        **errs_ << "checkSizes: The target object of the Import or Export is "
          "not fill complete.  Both source and target objects must be fill "
          "complete." << std::endl;
      }
      if (src->graph_.getColMap ().is_null ()) {
        **errs_ << "checkSizes: The source object of the Import or Export does "
          "not have a column Map.  Both source and target objects must have "
          "column Maps." << std::endl;
      }
      if (this->graph_.getColMap ().is_null ()) {
        **errs_ << "checkSizes: The target object of the Import or Export does "
          "not have a column Map.  Both source and target objects must have "
          "column Maps." << std::endl;
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
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    const this_type* src = dynamic_cast<const this_type* > (&source);
    if (src == NULL) {
      * (this->localError_) = true;
      if ((*errs_).is_null ()) {
        *errs_ = Teuchos::rcp (new std::ostringstream ());
      }
      **errs_ << "copyAndPermute: The source object of the Import or Export is "
        "either not a BlockCrsMatrix, or does not have the right template "
        "parameters.  checkSizes() should have caught this.  "
        "Please report this bug to the Tpetra developers." << std::endl;
      // There's no communication in this method, so it's safe just to
      // return on error.
      return;
    }

    bool srcInvalidRow = false;
    bool dstInvalidCol = false;

    // Copy the initial sequence of rows that are the same.
    //
    // The two graphs might have different column Maps, so we need to
    // do this using global column indices.  This is purely local, so
    // we only need to check for local sameness of the two column
    // Maps.

    const map_type& srcMap = * (src->graph_.getColMap ());
    const map_type& tgtMap = * (this->graph_.getColMap ());
    const bool canUseLocalColumnIndices = srcMap.locallySameAs (tgtMap);
    const size_t numPermute =
      std::min (permuteToLIDs.size (), permuteFromLIDs.size ());

    if (canUseLocalColumnIndices) {
      // Copy local rows that are the "same" in both source and target.
      for (size_t localRow = 0; localRow < numSameIDs; ++localRow) {
        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        // If this call fails, that means the mesh row local index is
        // invalid.  That means the Import or Export is invalid somehow.
        LO err = src->getLocalRowView (localRow, lclSrcCols, vals, numEntries);
        if (err != 0) {
          srcInvalidRow = true;
        }
        else {
          err = this->replaceLocalValues (localRow, lclSrcCols, vals, numEntries);
          if (err != numEntries) {
            dstInvalidCol = true;
          }
        }
      }

      // Copy the "permute" local rows.
      for (size_t k = 0; k < numPermute; ++k) {
        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        LO err = src->getLocalRowView (permuteFromLIDs[k], lclSrcCols, vals, numEntries);
        if (err != 0) {
          srcInvalidRow = true;
        }
        else {
          err = this->replaceLocalValues (permuteFromLIDs[k], lclSrcCols, vals, numEntries);
          if (err != numEntries) {
            dstInvalidCol = true;
          }
        }
      }
    }
    else {
      // Reserve space to store the destination matrix's local column indices.
      Teuchos::Array<LO> lclDstCols (src->graph_.getNodeMaxNumRowEntries ());

      // Copy local rows that are the "same" in both source and target.
      for (size_t localRow = 0; localRow < numSameIDs; ++localRow) {
        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        // If this call fails, that means the mesh row local index is
        // invalid.  That means the Import or Export is invalid somehow.
        LO err = src->getLocalRowView (localRow, lclSrcCols, vals, numEntries);
        if (err != 0) {
          srcInvalidRow = true;
        }
        else {
          // Convert the source matrix's local column indices to the
          // destination matrix's local column indices.
          Teuchos::ArrayView<LO> lclDstColsView = lclDstCols.view (0, numEntries);
          for (LO k = 0; k < numEntries; ++k) {
            lclDstColsView[k] = tgtMap.getLocalElement (srcMap.getGlobalElement (lclSrcCols[k]));
          }
          err = this->replaceLocalValues (localRow, lclDstColsView.getRawPtr (), vals, numEntries);
          if (err != numEntries) {
            dstInvalidCol = true;
          }
        }
      }

      // Copy the "permute" local rows.
      for (size_t k = 0; k < numPermute; ++k) {
        const LO* lclSrcCols;
        Scalar* vals;
        LO numEntries;
        LO err = src->getLocalRowView (permuteFromLIDs[k], lclSrcCols, vals, numEntries);
        if (err != 0) {
          srcInvalidRow = true;
        }
        else {
          // Convert the source matrix's local column indices to the
          // destination matrix's local column indices.
          Teuchos::ArrayView<LO> lclDstColsView = lclDstCols.view (0, numEntries);
          for (LO j = 0; j < numEntries; ++j) {
            lclDstColsView[j] = tgtMap.getLocalElement (srcMap.getGlobalElement (lclSrcCols[j]));
          }
          err = this->replaceLocalValues (permuteFromLIDs[k], lclDstColsView.getRawPtr (), vals, numEntries);
          if (err != numEntries) {
            dstInvalidCol = true;
          }
        }
      }
    }

    // FIXME (mfh 19 May 2014) Would be great to have a way to report
    // errors, without throwing an exception.  The latter would be bad
    // in the case of differing block sizes, just in case those block
    // sizes are inconsistent across processes.  (The latter is not
    // allowed, but possible due to user error.)

    if (srcInvalidRow || dstInvalidCol) {
      * (this->localError_) = true;
      if ((*errs_).is_null ()) {
        *errs_ = Teuchos::rcp (new std::ostringstream ());
      }
      **errs_ << "copyAndPermute: The graph structure of the source of the "
        "Import or Export must be a subset of the graph structure of the "
        "target." << std::endl;
    }
  }

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
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    const bool debug = true;
    const this_type* src = dynamic_cast<const this_type* > (&source);

    // Should have checked for this case in checkSizes().
    TEUCHOS_TEST_FOR_EXCEPTION(
      src == NULL, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "packAndPrepare: The source object of the Import or Export is either "
      "not a BlockCrsMatrix, or does not have the right template parameters.  "
      "checkSizes() should have caught this.  "
      "Please report this bug to the Tpetra developers.");

    if (debug) {
      const int myRank = graph_.getComm ()->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": packAndPrepare" << std::endl;
      std::cerr << os.str ();
    }

    const crs_graph_type& srcGraph = src->graph_;
    const size_t blockSize = static_cast<size_t> (src->getBlockSize ());
    const size_t bytesPerBlock = blockSize * blockSize * sizeof (Scalar);

    // Graphs and matrices are allowed to have a variable number of
    // entries per row.  We could test whether all rows have the same
    // number of entries, but DistObject can only use this
    // optimization if all rows on _all_ processes have the same
    // number of entries.  Rather than do the all-reduce necessary to
    // test for this unlikely case, we tell DistObject (by setting
    // constantNumPackets to zero) to assume that different rows may
    // have different numbers of entries.
    constantNumPackets = 0;

    // Count the number of packets per row.  The number of packets is
    // the number of entries to send.  However, we don't pack entries
    // contiguously.  Instead, we put the column indices first, and
    // then all the block values.  This makes it easier to align the
    // latter to sizeof (Scalar), and reduces the number of memcpy
    // operations necessary to copy in the row's data.  (We have to
    // use memcpy because of ANSI C(++) aliasing rules, which
    // compilers now enforce.)  We pack global column indices, because
    // the source and target matrices are allowed to have different
    // column Maps.

    // Byte alignment of the start of the block values in each row's
    // data.  We use the max in case sizeof(Scalar) < sizeof(GO),
    // e.g., float vs. long long.
    const size_t alignment = std::max (sizeof (Scalar), sizeof (GO));

    // While counting the number of packets per row, compute the total
    // required send buffer ('exports') size, so we can resize it if
    // needed.  Also count the total number of entries to send, as a
    // sanity check for the total buffer size.
    size_t totalBufSize = 0;
    size_t totalNumEnt = 0;

    for (size_t k = 0; k < static_cast<size_t> (exportLIDs.size ()); ++k) {
      const LO localRow = exportLIDs[k];
      const size_t numEnt = srcGraph.getNumEntriesInLocalRow (localRow);
      totalNumEnt += numEnt;

      // If any given LIDs are invalid, the above might return either
      // zero or the invalid size_t value.  If the former, we have no
      // way to tell, but that's OK; it just means the calling process
      // won't pack anything (it has nothing to pack anyway).  If the
      // latter, we replace it with zero (that row is not owned by the
      // calling process, so it has no entries to pack), and set the
      // local error flag.
      if (numEnt == Teuchos::OrdinalTraits<size_t>::invalid ()) {
        numPacketsPerLID[localRow] = static_cast<size_t> (0);
        * (this->localError_) = true;
      } else {
        numPacketsPerLID[localRow] = numEnt;
      }

      // Space for the global column indices in this row.
      totalBufSize += numEnt * sizeof (GO);

      // Padding, so that the block values are Scalar aligned.
      {
        const size_t rem = (numEnt * sizeof (GO)) % alignment;
        const size_t padding = (rem == 0) ?
          static_cast<size_t> (0) :
          alignment - rem;
        totalBufSize += padding;
      }

      // Space for the block values in this row.  We don't include
      // padding for vectorization when packing the blocks.
      totalBufSize += numEnt * bytesPerBlock;

      // In case sizeof(Scalar) < sizeof(GO) (e.g., float vs. long
      // long), we pad the end of the send buffer so that the next row
      // is aligned to sizeof(GO) bytes.
      {
        const size_t rem = (numEnt * bytesPerBlock) % sizeof (GO);
        const size_t padding = (rem == 0) ?
          static_cast<size_t> (0) :
          sizeof (GO) - rem;
        totalBufSize += padding;
      }
    }

    if (debug) {
      const int myRank = graph_.getComm ()->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": packAndPrepare: totalBufSize = "
         << totalBufSize << std::endl;
      std::cerr << os.str ();
    }

    {
      const size_t minSaneBufSize = totalNumEnt * (bytesPerBlock + sizeof (GO));
      // If the test triggers, there's no point in setting the error
      // flag here, because the above code is completely broken.  This
      // probably means that the code below is broken as well.
      TEUCHOS_TEST_FOR_EXCEPTION(
        minSaneBufSize < totalBufSize, std::logic_error, "Tpetra::Experimental"
        "::BlockCrsMatrix: This method must have computed the total send buffer"
        " size incorrectly.  The minimum size in bytes without padding for "
        "alignment is " << minSaneBufSize << ", but this method computed the "
        "size in bytes _with_ padding for alignment as a lesser number " <<
        totalBufSize << ".  Please report this bug to the Tpetra developers.");
    }

    // packAndPrepare is responsible for resizing the send buffer.
    exports.resize (totalBufSize);

    // In the loop below: Current offset in bytes into the 'exports'
    // array; where to start putting the current row's data.
    size_t exportsOffset = 0;

    // Temporary space for global column indices.
    Teuchos::Array<GO> gblColIndsSpace (srcGraph.getNodeMaxNumRowEntries ());

    // Temporary row-major contiguous buffer for a block's entries.
    Teuchos::Array<Scalar> tempBlockSpace (blockSize * blockSize);
    little_block_type tempBlock (tempBlockSpace.getRawPtr (), blockSize, blockSize, 1);

    // Source matrix's column Map.  We verified in checkSizes() that
    // the column Map exists (is not null).
    const map_type& srcColMap = * (srcGraph.getColMap ());

    // Pack the data for each row to send, into the 'exports' buffer.
    // If any given LIDs are invalid, we pack obviously invalid data
    // (e.g., invalid column indices) into the buffer for that LID,
    // and set the local error flag.
    for (size_t lidInd = 0; lidInd < static_cast<size_t> (exportLIDs.size ()); ++lidInd) {
      // Get a view of row exportLIDs[lidInd].
      const LO lclRowInd = exportLIDs[lidInd];
      const LO* lclColInds;
      Scalar* vals;
      LO numEnt;
      const int err = src->getLocalRowView (lclRowInd, lclColInds, vals, numEnt);
      if (err != 0) {
        * (localError_) = true;
        // TODO (mfh 20 May 2014) Report the local error, without
        // printing a line for each bad LID.  It might help to collect
        // all the bad LIDs, but don't print them all if there are too
        // many.
        continue;
      }

      // Convert column indices from local to global.
      Teuchos::ArrayView<GO> gblColInds = gblColIndsSpace.view (0, numEnt);
      for (size_t j = 0; j < static_cast<size_t> (numEnt); ++j) {
        gblColInds[j] = srcColMap.getGlobalElement (lclColInds[j]);
      }

      // Pack the column indices.  Use memcpy to follow ANSI aliasing rules.
      packet_type* exportsStart = &exports[exportsOffset];
      memcpy (exportsStart, gblColInds.getRawPtr (), numEnt * sizeof (GO));
      exportsOffset += numEnt * sizeof (GO);

      // Compute padding so that block values are Scalar aligned.
      {
        const size_t rem = (numEnt * sizeof (GO)) % alignment;
        const size_t padding = (rem == 0) ?
          static_cast<size_t> (0) :
          alignment - rem;
        exportsOffset += padding;
      }

      // Pack the block values in this row.  Always pack row major,
      // for consistency, in case the source and target objects order
      // their blocks differently.
      //
      // In order to follow ANSI aliasing rules that forbid certain
      // kinds of type punning, we copy each block's entries first
      // into a temporary contiguous buffer, and then memcpy into the
      // send buffer.  We could just memcpy one entry at a time, but
      // it's probably faster to avoid the library call.
      const size_t rowStart = src->ptr_[lclRowInd];
      for (size_t j = 0; j < static_cast<size_t> (numEnt); ++j) {
        const_little_block_type srcBlock =
          src->getConstLocalBlockFromAbsOffset (rowStart + j);
        if (static_cast<size_t> (srcBlock.getBlockSize ()) == blockSize) {
          // A block size of zero is a possible error state.  Of
          // course, the actual block size could be zero.  That would
          // be silly, but why shouldn't it be legal?  That's why we
          // check whether the actual block size is different.
          * (localError_) = true;
          // Pack a block of zeros.  It might make sense to pack NaNs,
          // if Scalar implements a NaN value.
          tempBlock.fill (STS::zero ());
        } else {
          tempBlock.assign (srcBlock);
        }
        memcpy (&exports[exportsOffset], tempBlock.getRawPtr (), bytesPerBlock);
        exportsOffset += bytesPerBlock;
      }

      // In case sizeof(Scalar) < sizeof(GO) (e.g., float vs. long
      // long), we pad the end of the send buffer so that the next row
      // is aligned to sizeof(GO) bytes.
      {
        const size_t rem = (numEnt * bytesPerBlock) % sizeof (GO);
        const size_t padding = (rem == 0) ?
          static_cast<size_t> (0) :
          sizeof (GO) - rem;
        exportsOffset += padding;
      }
    } // for each LID (of a row) to send


    {
      const int myRank = graph_.getComm ()->getRank ();
      std::ostringstream os;
      os << "Proc " << myRank << ": packAndPrepare done" << std::endl;
      std::cerr << os.str ();
    }
  }


  template<class Scalar, class LO, class GO, class Node>
  void
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  unpackAndCombine (const Teuchos::ArrayView<const LO> &importLIDs,
                    const Teuchos::ArrayView<const packet_type> &imports,
                    const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                    size_t /* constantNumPackets */, // not worthwhile to use this
                    Tpetra::Distributor& /* distor */,
                    Tpetra::CombineMode CM)
  {
    if (CM != ADD && CM != INSERT && CM != REPLACE && CM != ABSMAX && CM != ZERO) {
      * (this->localError_) = true;
      if ((*errs_).is_null ()) {
        *errs_ = Teuchos::rcp (new std::ostringstream ());
        **errs_ << "unpackAndCombine: Invalid CombineMode value " << CM << ".  "
                << "Valid values include ADD, INSERT, REPLACE, ABSMAX, and ZERO."
                << std::endl;
        // It won't cause deadlock to return here, since this method
        // does not communicate.
        return;
      }
    }

    if (CM == ZERO) {
      return; // nothing to do; no need to combine entries
    }

    const size_t blockSize = this->getBlockSize ();
    const size_t bytesPerBlock = blockSize * blockSize * sizeof (Scalar);
    const size_t numImportLids = importLIDs.size ();

    // Byte alignment of the start of the block values in each row's
    // data.  We use the max in case sizeof(Scalar) < sizeof(GO),
    // e.g., float vs. long long.
    const size_t alignment = std::max (sizeof (Scalar), sizeof (GO));

    // Temporary space to cache local and global column indices.
    // Column indices come in as global indices, in case the source
    // object's column Map differs from the target object's (this's)
    // column Map.  We need the latter space in order to follow ANSI
    // aliasing rules (we don't want to type pun the buffer).
    Teuchos::Array<LO> lclColIndsSpace (graph_.getNodeMaxNumRowEntries ());
    Teuchos::Array<GO> gblColIndsSpace (graph_.getNodeMaxNumRowEntries ());

    // Temporary contiguous buffer for a row's block entries.
    Teuchos::Array<Scalar> tempVals (graph_.getNodeMaxNumRowEntries ());

    // Target matrix's column Map.  Use to convert the global column
    // indices in the receive buffer to local indices.  We verified in
    // checkSizes() that the column Map exists (is not null).
    const map_type& tgtColMap = * (this->graph_.getColMap ());

    // In the loop below: Current offset (in bytes) into the 'imports'
    // array; from whence to unpack data for the current LID (row).
    size_t importsOffset = 0;

    for (size_t lidInd = 0; lidInd < numImportLids; ++lidInd) {
      const size_t numEnt = numPacketsPerLID[lidInd]; // # entries in row
      const LO lclRowInd = importLIDs[lidInd]; // current local row index

      // Lengths and offsets for packed data in the 'imports' array.
      //
      // We pack all column indices first, then all values.  We pad
      // the former so that the latter is aligned to sizeof(Scalar).
      // 'padding' gives the size in bytes of the padding.
      const size_t numIndexBytes = numEnt * sizeof (GO);
      const size_t numBlockBytes = numEnt * bytesPerBlock;

      // Views for storing global and local column indices.
      Teuchos::ArrayView<GO> gblColInds = gblColIndsSpace.view (0, numEnt);
      Teuchos::ArrayView<LO> lclColInds = lclColIndsSpace.view (0, numEnt);

      // Unpack the global column indices from the receive buffer.
      memcpy (gblColInds.getRawPtr (), &imports[importsOffset], numIndexBytes);
      // Convert global column indices to local.
      for (size_t j = 0; j < numEnt; ++j) {
        lclColInds[j] = tgtColMap.getLocalElement (gblColInds[j]);
      }

      // Update the byte offset into the receive buffer.
      importsOffset += numIndexBytes;
      // Block values are aligned to max(sizeof(Scalar), sizeof(GO)).
      // Update the byte offset to account for padding to this alignment.
      {
        const size_t rem = numIndexBytes % alignment;
        const size_t padding = (rem == 0) ?
          static_cast<size_t> (0) :
          alignment - rem;
        importsOffset += padding;
      }

      // Copy out data for _all_ the blocks.  We memcpy into temp
      // storage, in order to follow ANSI aliasing rules that forbid
      // certain kinds of type punning.
      memcpy (tempVals.getRawPtr (), &imports[importsOffset], numBlockBytes);
      // Update the byte offset into the receive buffer.
      importsOffset += numBlockBytes;
      // In case sizeof(Scalar) < sizeof(GO) (e.g., float vs. long
      // long), the receive buffer was padded so that the next row is
      // aligned to sizeof(GO) bytes.
      const size_t rem = numBlockBytes % sizeof (GO);
      const size_t padding = (rem == 0) ?
        static_cast<size_t> (0) :
        sizeof (GO) - rem;
      importsOffset += padding;

      // Combine the incoming data with the matrix's current data.
      LO successCount = 0;
      if (CM == ADD) {
        successCount =
          this->sumIntoLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                    tempVals.getRawPtr (), numEnt);
      } else if (CM == INSERT || CM == REPLACE) {
        successCount =
          this->replaceLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                    tempVals.getRawPtr (), numEnt);
      } else if (CM == ABSMAX) {
        successCount =
          this->absMaxLocalValues (lclRowInd, lclColInds.getRawPtr (),
                                   tempVals.getRawPtr (), numEnt);
      }
      // We've already checked that CM is valid.

      if (static_cast<size_t> (successCount) != numEnt) {
        * (localError_) = true;
      }
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
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
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
