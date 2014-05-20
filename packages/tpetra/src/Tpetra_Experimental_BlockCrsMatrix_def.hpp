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
    X_colMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    Y_rowMap_ (new Teuchos::RCP<BMV> ()), // ptr to a null ptr
    columnPadding_ (0), // no padding by default
    rowMajor_ (true), // row major blocks by default
    localError_ (false)
  {}

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
    localError_ (false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

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
    localError_ (false)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph_.isSorted (), std::invalid_argument, "Tpetra::Experimental::"
      "BlockCrsMatrix constructor: The input CrsGraph does not have sorted "
      "rows (isSorted() is false).  This class assumes sorted rows.");

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

    for (size_t k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
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

    for (size_t k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
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

    for (size_t k = 0; k < numColInds; ++k, pointOffset += perBlockSize) {
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

    for (size_t k = 0; k < numColInds; ++k) {
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

    for (size_t k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
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

    for (size_t k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
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

    for (size_t k = 0; k < numOffsets; ++k, pointOffset += perBlockSize) {
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
  LO
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
        X_colMap = **X_colMap_; // MUST do a shallow copy
      }

      BMV Y_rowMap;
      if (import.is_null ()) {
        Y_rowMap = Y; // MUST do a shallow copy
      } else if ((*Y_rowMap_).is_null () ||
                 (**Y_rowMap_).getNumVectors () != Y.getNumVectors () ||
                 (**Y_rowMap_).getBlockSize () != Y.getBlockSize ()) {
        *Y_rowMap_ = rcp (new BMV (* (graph_.getRowMap ()), getBlockSize (),
                                   static_cast<LO> (X.getNumVectors ())));
        Y_rowMap = **Y_rowMap_; // MUST do a shallow copy
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
          Y_cur.matvecUpdate (alpha, A_cur, X_cur); // Y_cur += alpha*A_cur*X_cur;
        }
      }
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
            Y_cur.matvecUpdate (alpha, A_cur, X_cur); // Y_cur += alpha*A_cur*X_cur;
          }
        }
      }
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
    if (hint < numEntriesInRow && ind_[absStartOffset] == colIndexToFind) {
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
  bool
  BlockCrsMatrix<Scalar, LO, GO, Node>::
  checkSizes (const Tpetra::SrcDistObject& source)
  {
    typedef BlockCrsMatrix<Scalar, LO, GO, Node> this_type;
    const this_type* src = dynamic_cast<const this_type* > (&source);

    // FIXME (mfh 19 May 2014) Would be great to have a way to report
    // errors, without throwing an exception.  The latter would be bad
    // in the case of differing block sizes, just in case those block
    // sizes are inconsistent across processes.  (The latter is not
    // allowed, but possible due to user error.)

    // TEUCHOS_TEST_FOR_EXCEPTION(
    //   src == NULL, std::invalid_argument, "Tpetra::Experimental::"
    //   "BlockCrsMatrix::checkSizes: The source object of the Import or Export "
    //   "must be a BlockCrsMatrix with the same template parameters as the "
    //   "target object.");
    // TEUCHOS_TEST_FOR_EXCEPTION(
    //   src->getBlockSize () != this->getBlockSize (), std::invalid_argument,
    //   "Tpetra::Experimental::BlockCrsMatrix::checkSizes: "
    //   "The source and target objects of the Import or Export must have the "
    //   "same block sizes.  The source's block size = " << src->getBlockSize ()
    //   << " != the target's block size = " << this->getBlockSize () << ".");

    if (src == NULL || src->getBlockSize () != this->getBlockSize ()) {
      const_cast<this_type*> (this)->localError_ = true;
      return false;
    }
    else {
      const_cast<this_type*> (this)->localError_ = false;

      // FIXME (mfh 19 May 2014) Return false for now, as a way to
      // signal that the DistObject functions haven't been tested yet.
      // Returning false tells the Import or Export not to proceed.
      // Once we finish and test those functions, we need to change this
      // to return true.
      return false;
    }
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      src == NULL, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "copyAndPermute: The source object of the Import or Export is either "
      "not a BlockCrsMatrix, or does not have the right template parameters.  "
      "checkSizes() should have caught this.  "
      "Please report this bug to the Tpetra developers.");

    bool srcInvalidRow = false;
    bool dstInvalidCol = false;
    for (size_t localRow = 0; localRow < numSameIDs; ++localRow) {
      const LO* localCols;
      Scalar* vals;
      LO numEntries;
      // If this call fails, that means the mesh row local index is
      // invalid.  That means the Import or Export is invalid somehow.
      LO err = src->getLocalRowView (localRow, localCols, vals, numEntries);
      if (err != 0) {
        srcInvalidRow = true;
      }
      else {
        err = this->replaceLocalValues (localRow, localCols, vals, numEntries);
        if (err != numEntries) {
          dstInvalidCol = true;
        }
      }
    }

    const size_t numPermute =
      std::min (permuteToLIDs.size (), permuteFromLIDs.size ());
    for (size_t k = 0; k < numPermute; ++k) {
      const LO* localCols;
      Scalar* vals;
      LO numEntries;
      LO err = src->getLocalRowView (permuteFromLIDs[k], localCols, vals, numEntries);
      if (err != 0) {
        srcInvalidRow = true;
      }
      else {
        err = this->replaceLocalValues (permuteFromLIDs[k], localCols, vals, numEntries);
        if (err != numEntries) {
          dstInvalidCol = true;
        }
      }
    }

    if (srcInvalidRow || dstInvalidCol) {
      localError_ = true;
    }

    // FIXME (mfh 19 May 2014) Would be great to have a way to report
    // errors, without throwing an exception.  The latter would be bad
    // in the case of differing block sizes, just in case those block
    // sizes are inconsistent across processes.  (The latter is not
    // allowed, but possible due to user error.)

    // TEUCHOS_TEST_FOR_EXCEPTION(
    //   srcInvalidRow || dstInvalidCol, std::runtime_error,
    //   "Tpetra::Experimental::BlockCrsMatrix::copyAndPermute: The graph "
    //   "structure of the source of the Import or Export must be a subset of the "
    //   "graph structure of the target.");
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
    const this_type* src = dynamic_cast<const this_type* > (&source);
    // Should have checked for this case in checkSizes().
    TEUCHOS_TEST_FOR_EXCEPTION(
      src == NULL, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
      "packAndPrepare: The source object of the Import or Export is either "
      "not a BlockCrsMatrix, or does not have the right template parameters.  "
      "checkSizes() should have caught this.  "
      "Please report this bug to the Tpetra developers.");

    const crs_graph_type& srcGraph = src->graph_;

    // Compute the number of packets per row.  Each packet includes
    // the _global_ column index and the block entry.  (We need the
    // latter, because the two matrices may have different column
    // Maps; this is allowed.)  While doing this, compute the total
    // required send buffer ('exports') size, so we can resize it if
    // needed.

    // Count the number of entries in each row to pack.
    //
    // Graphs and matrices are allowed to have a variable number of
    // entries per row.  We could test whether all rows have the same
    // number of entries, but DistObject can only use this
    // optimization if all rows on _all_ processes have the same
    // number of entries.  Rather than do the all-reduce necessary to
    // test for this unlikely case, we tell DistObject to assume that
    // different rows may have different numbers of entries.
    constantNumPackets = 0; // variable number of packets (entries) per row
    size_t totalNumPackets = 0;
    for (size_t k = 0; k < exportLIDs.size (); ++k) {
      const LO localRow = exportLIDs[k];
      const size_t numEnt = srcGraph.getNumEntriesInLocalRow (localRow);
      numPacketsPerLID[localRow] = numEnt;
      totalNumPackets += numEnt;
    }
    exports.resize (totalNumPackets); // resize the send buffer

    const size_t blockSize = static_cast<size_t> (src->getBlockSize ());

    // Compute the size in bytes of each packet.
    //
    // Each packet contains a global (column) index, and the
    // corresponding block's entries.  We don't pack any padding in
    // the block.  Packed blocks are stored in row-major order,
    // regardless of the matrix's storage format.  This lets us Import
    // or Export between two matrices with different block storage
    // formats or padding specifications, though it's unlikely we
    // would ever want to exercise this case.
    const size_t packetSize =
      sizeof (GO) + blockSize * blockSize * sizeof (Scalar);

    // Pack each block entry to send into the 'exports' buffer.
    // First pack the column index, then the block's entries.
    size_t exportsOffset = 0;
    // If any given LIDs are invalid, we pack obviously invalid data
    // (e.g., invalid column indices) into the buffer for that LID.
    //
    // TODO (mfh 17 May 2014) It would also be wise to set a local
    // error flag, in case some processes don't get the invalid column
    // indices.  That way, we could propagate the error state on the
    // next all-reduce.
    for (size_t k = 0; k < exportLIDs.size (); ++k, exportsOffset += packetSize) {
      //
      // Get a view of row exportLIDs[k].
      //
      const LO localRow = exportLIDs[k];
      const LO* localColInds;
      Scalar* vals;
      LO numEntries;
      const int err = src->getLocalRowView (localRow, localColInds, vals, numEntries);
      //
      // Pack the entries in the row.
      //
      packet_type* exportsStart = &exports[exportsOffset];
      // Don't access ptr_[localRow] if localRow is an invalid LID.
      const size_t rowStart = (err == 0) ? src->ptr_[localRow] : static_cast<size_t> (0);
      for (LO j = 0; j < numEntries; ++j) {
        const GO gblCol = (err == 0) ?
          Teuchos::OrdinalTraits<GO>::invalid () :
          rowMeshMap_.getGlobalElement (localColInds[j]);
        // memcpy is the only safe function to use, given ANSI aliasing rules.
        memcpy (exportsStart, &gblCol, sizeof (GO));
        // FIXME (mfh 17 May 2014) Does this break ANSI aliasing rules?
        // FIXME (mfh 17 May 2014) Do Scalar values need to be aligned?
        Scalar* tgtPtr = reinterpret_cast<Scalar*> (exportsStart + sizeof (GO));

        // Pack row major, regardless of the state of the source or
        // target matrix, for consistency.
        little_block_type tgtBlk (tgtPtr, blockSize, blockSize, 1);
        if (err == 0) {
          const_little_block_type srcBlk =
            src->getConstLocalBlockFromAbsOffset (rowStart + j);
          tgtBlk.assign (srcBlk);
        }
        else {
          // NOTE (mfh 17 May 2014) It might be a good idea to use
          // NaNs here instead, if Scalar supports them, to indicate
          // invalid data.  However, we've already set an invalid
          // column index above.
          tgtBlk.fill (STS::zero ());
        }
      }
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
    using Teuchos::ArrayView;
    using Teuchos::av_reinterpret_cast;

    if (CM == ZERO) {
      return; // nothing to do; no need to combine entries
    }

    const size_t blockSize = this->getBlockSize ();
    const size_t numImportLids = importLIDs.size ();

    // Temporary space to cache local column indices.  Column indices
    // come in as global indices, in case the source object's column
    // Map differs from the target object's (this's) column Map.
    Teuchos::Array<LO> lclColInds (graph_.getNodeMaxNumRowEntries ());

    // Current offset (in bytes) into the 'imports' array.  We unpack
    // data from the 'imports' array into the target object.
    size_t importsOffset = 0;
    for (size_t importLidIndex = 0; importLidIndex < numImportLids; ++importLidIndex) {
      // Number of entries in the row.
      const size_t numEntries = numPacketsPerLID[importLidIndex];
      // The current (local) row index.
      const LO importLid = importLIDs[importLidIndex];

      // Lengths and offsets for packed data in the 'imports' array.
      //
      // We pack all column indices first, then all values.  We pad
      // the former so that the latter is aligned to sizeof(Scalar).
      // 'padding' gives the size in bytes of the padding.
      const size_t numIndexBytes = numEntries * sizeof (GO);
      const size_t numBlockBytes = numEntries * blockSize * sizeof (Scalar);
      const size_t padding = sizeof (Scalar) - (numIndexBytes % sizeof (Scalar));
      const size_t absBlockOffset = importsOffset + numIndexBytes + padding;

      ArrayView<const GO> gblColInds =
        av_reinterpret_cast<const GO> (imports.view (importsOffset, numIndexBytes));
      ArrayView<const Scalar> vals =
        av_reinterpret_cast<const Scalar> (imports.view (absBlockOffset, numBlockBytes));

      ArrayView<LO> lclColIndsView = lclColInds.view (0, numEntries);
      for (size_t k = 0; k < numEntries; ++k) {
        lclColIndsView[k] = rowMeshMap_.getLocalElement (gblColInds[k]);
      }

      LO successCount = 0;
      if (CM == ADD) {
        successCount =
          this->sumIntoLocalValues (importLid, lclColIndsView.getRawPtr (),
                                    vals.getRawPtr (), numEntries);
      } else if (CM == INSERT || CM == REPLACE) {
        successCount =
          this->replaceLocalValues (importLid, lclColIndsView.getRawPtr (),
                                    vals.getRawPtr (), numEntries);
      } else if (CM == ABSMAX) {
        successCount =
          this->absMaxLocalValues (importLid, lclColIndsView.getRawPtr (),
                                   vals.getRawPtr (), numEntries);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error, "Tpetra::Experimental::BlockCrsMatrix::"
          "unpackAndCombine: Invalid CombineMode value " << CM << ".  Valid "
          "values include ADD, INSERT, REPLACE, ABSMAX, and ZERO.");
      }

      if (successCount != numEntries) {
        localError_ = true;
      }
      importsOffset += numIndexBytes + padding + numBlockBytes;
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
