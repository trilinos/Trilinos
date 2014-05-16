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

#ifdef DOXYGEN_USE_ONLY
#  include "Tpetra_Experimental_BlockMultiVector_decl.hpp"
#endif

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
    X_colMap_initialized_ (false),
    Y_rowMap_initialized_ (false),
    columnPadding_ (0), // no padding by default
    rowMajor_ (true) // row major blocks by default
  {}

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    blockSize_ (blockSize),
    ptr_ (NULL), // to be initialized below
    ind_ (NULL), // to be initialized below
    X_colMap_initialized_ (false),
    Y_rowMap_initialized_ (false),
    columnPadding_ (0), // no padding by default
    rowMajor_ (true) // row major blocks by default
  {
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
    val_.resize (graph.getNodeNumEntries () * allocationSizePerBlock ());
  }

  template<class Scalar, class LO, class GO, class Node>
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  BlockCrsMatrix (const crs_graph_type& graph,
                  const map_type& domainPointMap,
                  const map_type& rangePointMap,
                  const LO blockSize) :
    dist_object_type (graph.getMap ()),
    graph_ (graph),
    domainPointMap_ (domainPointMap),
    rangePointMap_ (rangePointMap),
    blockSize_ (blockSize),
    ptr_ (NULL), // to be initialized below
    ind_ (NULL), // to be initialized below
    X_colMap_initialized_ (false),
    Y_rowMap_initialized_ (false),
    columnPadding_ (0), // no padding by default
    rowMajor_ (true) // row major blocks by default
  {
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
    val_.resize (graph.getNodeNumEntries () * blockSize * blockSize);
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
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  replaceLocalValues (const LO localRowInd,
                      const LO colInds[],
                      const Scalar vals[],
                      const LO numColInds) const
  {
    if (! graph_.getRowMap ()->isNodeLocalElement (localRowInd)) {
      // We replaced no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }

    const size_t perBlockSize = static_cast<LO> (allocationSizePerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t hint = 0; // Guess for the offset into allColIndsInRow
    size_t valsPointOffset = 0; // Current offset into vals
    LO validCount = 0; // number of valid column indices in colInds

    for (size_t k = 0; k < numColInds; ++k, valsPointOffset += perBlockSize) {
      const size_t blockOffset =
        findOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (blockOffset != STINV) {
        little_block_type A_old = getNonConstLocalBlockFromOffset (blockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vals, valsPointOffset);
        A_old.assign (A_new);
        hint = blockOffset + 1;
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
    if (! graph_.getRowMap ()->isNodeLocalElement (localRowInd)) {
      // We replaced no values, because the input local row index is
      // invalid on the calling process.  That may not be an error, if
      // numColInds is zero anyway; it doesn't matter.  This is the
      // advantage of returning the number of valid indices.
      return static_cast<LO> (0);
    }

    const size_t perBlockSize = static_cast<LO> (allocationSizePerBlock ());
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    size_t hint = 0; // Guess for the offset into allColIndsInRow
    size_t valsPointOffset = 0; // Current offset into vals
    LO validCount = 0; // number of valid column indices in colInds

    for (size_t k = 0; k < numColInds; ++k, valsPointOffset += perBlockSize) {
      const size_t blockOffset =
        findOffsetOfColumnIndex (localRowInd, colInds[k], hint);
      if (blockOffset != STINV) {
        little_block_type A_old = getNonConstLocalBlockFromOffset (blockOffset);
        const_little_block_type A_new =
          getConstLocalBlockFromInput (vals, valsPointOffset);
        A_old.update (STS::one (), A_new);
        hint = blockOffset + 1;
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "getLocalRowView: Not yet implemented.");
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowOffsets (const LO localRowInd,
                      ptrdiff_t offsets[],
                      const LO colInds[],
                      const LO numColInds) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "getLocalRowOffsets(4 args): Not yet implemented.");
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getLocalRowOffsets (const LO localRowInd,
                      ptrdiff_t offsets[]) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "getLocalRowOffsets(2 args): Not yet implemented.");
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  replaceLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "replaceLocalValuesByOffsets: Not yet implemented.");
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  sumIntoLocalValuesByOffsets (const LO localRowInd,
                               const ptrdiff_t offsets[],
                               const Scalar vals[],
                               const LO numOffsets) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "sumIntoLocalValuesByOffsets: Not yet implemented.");
  }


  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNumEntriesInLocalRow (const LO localRowInd) const
  {
    const size_t numEntInGraph = graph_.getNumEntriesInLocalRow (localRowInd);
    if (numEntInGraph == Teuchos::OrdinalTraits<size_t>::invalid ()) {
      return static_cast<LO> (0); // the calling process doesn't have that row
    } else {
      return getBlockSize () * static_cast<LO> (numEntInGraph);
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
    typedef Tpetra::Import<LO, GO, Node> import_type;
    typedef Tpetra::Export<LO, GO, Node> export_type;
    const Scalar zero = STS::zero ();
    const Scalar one = STS::one ();
    RCP<const import_type> import = graph_.getImporter ();
    // "export" is a reserved C++ keyword, so we can't use it.
    RCP<const export_type> theExport = graph_.getExporter ();

    if (alpha == zero) {
      if (beta == zero) {
        Y.putScalar (zero); // replace Inf or NaN (BLAS rules)
      } else if (beta != one) {
        Y.scale (beta);
      }
    } else { // alpha != 0
      // BMV gets view semantics from the beginning.
      BMV X_colMap;
      if (import.is_null ()) {
        X_colMap = X; // MUST do a shallow copy
      } else {
        if (! X_colMap_initialized_ ||
            X_colMap_.getNumVectors () != X.getNumVectors () ||
            X_colMap_.getBlockSize () != X.getBlockSize ()) {
          X_colMap_ = BMV (* (graph_.getColMap ()), getBlockSize (),
                           static_cast<LO> (X.getNumVectors ()));
        }
        X_colMap_.doImport (X, *import, Tpetra::REPLACE);
        X_colMap = X_colMap_; // MUST do a shallow copy
      }

      BMV Y_rowMap;
      if (import.is_null ()) {
        Y_rowMap = Y; // MUST do a shallow copy
      } else if (! Y_rowMap_initialized_ ||
                 Y_rowMap_.getNumVectors () != Y.getNumVectors () ||
                 Y_rowMap_.getBlockSize () != Y.getBlockSize ()) {
        Y_rowMap_ = BMV (* (graph_.getRowMap ()), getBlockSize (),
                         static_cast<LO> (X.getNumVectors ()));
      }

      localApplyBlockNoTrans (X_colMap, Y_rowMap, alpha, beta);

      if (! theExport.is_null ()) {
        Y.doExport (Y_rowMap, *theExport, Tpetra::REPLACE);
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
    const size_t* const ptr = graph_.getNodeRowPtrs ().getRawPtr ();
    const LO* const ind = graph_.getNodePackedIndices ().getRawPtr ();
    const LO numLocalMeshRows = static_cast<LO> (graph_.getNodeNumRows ());
    const LO numVecs = static_cast<LO> (X.getNumVectors ());

    // If using (new) Kokkos, replace localMem with thread-local
    // memory.  Note that for larger block sizes, this will affect the
    // two-level parallelization.  Look to Stokhos for best practice
    // on making this fast for GPUs.
    const LO blockSize = getBlockSize ();
    Teuchos::Array<Scalar> localMem (blockSize);
    little_vec_type Y_lcl (localMem.getRawPtr (), blockSize, 1);

    if (numVecs == 1) {
      for (LO localMeshRow = 0; localMeshRow < numLocalMeshRows; ++localMeshRow) {
        little_vec_type Y_cur = Y.getLocalBlock (localMeshRow, 0);
        if (beta == zero) {
          Y_lcl.fill (zero);
        } else if (beta == one) {
          Y_lcl.assign (Y_cur);
        } else {
          Y_lcl.assign (Y_cur);
          Y_lcl.scale (beta);
        }
        const size_t meshStart = ptr[localMeshRow];
        const size_t meshEnd = ptr[localMeshRow+1];
        for (size_t meshEntryIndex = meshStart; meshEntryIndex < meshEnd; ++meshEntryIndex) {
          const LO meshCol = ind[meshEntryIndex];
          const_little_block_type A_cur = getConstLocalBlockFromOffset (meshEntryIndex);
          little_vec_type X_cur = X.getLocalBlock (meshCol, 0);
          Y_cur.matvecUpdate (alpha, A_cur, X_cur); // Y_cur += alpha*A_cur*X_cur;
        }
      }
    }
    else {
      for (LO localMeshRow = 0; localMeshRow < numLocalMeshRows; ++localMeshRow) {
        for (LO j = 0; j < numVecs; ++j) {
          little_vec_type Y_cur = Y.getLocalBlock (localMeshRow, j);
          if (beta == zero) {
            Y_lcl.fill (zero);
          } else if (beta == one) {
            Y_lcl.assign (Y_cur);
          } else {
            Y_lcl.assign (Y_cur);
            Y_lcl.scale (beta);
          }
          const size_t meshStart = ptr[localMeshRow];
          const size_t meshEnd = ptr[localMeshRow+1];
          for (size_t meshEntryIndex = meshStart; meshEntryIndex < meshEnd; ++meshEntryIndex) {
            const LO meshCol = ind[meshEntryIndex];
            const_little_block_type A_cur = getConstLocalBlockFromOffset (meshEntryIndex);
            little_vec_type X_cur = X.getLocalBlock (meshCol, j);
            Y_cur.matvecUpdate (alpha, A_cur, X_cur); // Y_cur += alpha*A_cur*X_cur;
          }
        }
      }
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  size_t
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  findOffsetOfColumnIndex (const LO localRowIndex,
                           const LO colIndexToFind,
                           const size_t hint) const
  {
    const size_t startOffset = ptr_[localRowIndex];
    const size_t endOffset = ptr_[localRowIndex+1];
    const size_t numEntriesInRow = endOffset - startOffset;
    const LO* beg = ind_ + startOffset;
    // If the hint was correct, then the hint is the offset to return.
    if (hint < numEntriesInRow && beg[hint] == colIndexToFind) {
      // Always return the _absolute_ offset, not the offset relative
      // to the current row.
      return startOffset + hint;
    }

    // The hint was wrong, so we must search for the given column
    // index in the column indices for the given row.  How we do the
    // search depends on whether the graph's column indices are
    // sorted.
    const LO* end = ind_ + endOffset;
    const LO* ptr = end;
    bool found = true;

    if (graph_.isSorted ()) { // use binary search
      std::pair<const LO*, const LO*> p =
        std::equal_range (beg, end, colIndexToFind);
      if (p.first == p.second) {
        found = false;
      } else {
        ptr = p.first;
      }
    } else { // use linear search
      ptr = std::find (beg, end, colIndexToFind);
      if (ptr == end) {
        found = false;
      }
    }

    if (found) {
      return startOffset + static_cast<size_t> (ptr - beg);
    } else {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  LO
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  allocationSizePerBlock () const
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
  getConstLocalBlockFromOffset (const size_t blockOffset) const
  {
    if (blockOffset >= ptr_[graph_.getNodeNumRows ()]) {
      // An empty block signifies an error.  We don't expect to see
      // this error in correct code, but it's helpful for avoiding
      // segfaults or memory corruption in case there is a bug.
      return const_little_block_type (NULL, 0, 0, 0);
    } else {
      const size_t pointOffset = blockOffset * allocationSizePerBlock ();
      return getConstLocalBlockFromInput (val_.getRawPtr (), pointOffset);
    }
  }

  template<class Scalar, class LO, class GO, class Node>
  typename BlockCrsMatrix<Scalar, LO, GO, Node>::little_block_type
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  getNonConstLocalBlockFromOffset (const size_t blockOffset) const
  {
    if (blockOffset >= ptr_[graph_.getNodeNumRows ()]) {
      // An empty block signifies an error.  We don't expect to see
      // this error in correct code, but it's helpful for avoiding
      // segfaults or memory corruption in case there is a bug.
      return little_block_type (NULL, 0, 0, 0);
    } else {
      const size_t pointOffset = blockOffset * allocationSizePerBlock ();
      return getNonConstLocalBlockFromInput (const_cast<Scalar*> (val_.getRawPtr ()), pointOffset);
    }
  }

} // namespace Experimental
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra::Experimental namespace!
//
#define TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_INSTANT(S,LO,GO,NODE) \
  template class BlockCrsMatrix< S, LO, GO, NODE >;

#endif // TPETRA_EXPERIMENTAL_BLOCKCRSMATRIX_DEF_HPP
