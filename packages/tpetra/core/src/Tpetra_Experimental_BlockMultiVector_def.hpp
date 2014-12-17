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

#ifndef TPETRA_EXPERIMENTAL_BLOCKMULTIVECTOR_DEF_HPP
#define TPETRA_EXPERIMENTAL_BLOCKMULTIVECTOR_DEF_HPP

#include "Tpetra_Experimental_BlockMultiVector_decl.hpp"

namespace Tpetra {
namespace Experimental {

template<class Scalar, class LO, class GO, class Node>
typename BlockMultiVector<Scalar, LO, GO, Node>::mv_type
BlockMultiVector<Scalar, LO, GO, Node>::
getMultiVectorView ()
{
  // Make sure that mv_ has view semantics.
  mv_.setCopyOrView (Teuchos::View);
  // Now the one-argument copy constructor will make a shallow copy,
  // and those view semantics will persist in all of its offspring.
  return mv_;
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const BlockMultiVector<Scalar, LO, GO, Node> >
BlockMultiVector<Scalar, LO, GO, Node>::
getBlockMultiVectorFromSrcDistObject (const Tpetra::SrcDistObject& src)
{
  typedef BlockMultiVector<Scalar, LO, GO, Node> BMV;
  const BMV* src_bmv = dynamic_cast<const BMV*> (&src);
  TEUCHOS_TEST_FOR_EXCEPTION(
    src_bmv == NULL, std::invalid_argument, "Tpetra::Experimental::"
    "BlockMultiVector: The source object of an Import or Export to a "
    "BlockMultiVector, must also be a BlockMultiVector.");
  return Teuchos::rcp (src_bmv, false); // nonowning RCP
}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const map_type& meshMap,
                  const LO blockSize,
                  const LO numVecs) :
  dist_object_type (Teuchos::rcp (new map_type (meshMap))), // shallow copy
  meshMap_ (meshMap),
  pointMap_ (makePointMap (meshMap, blockSize)),
  mv_ (Teuchos::rcpFromRef (pointMap_), numVecs), // nonowning RCP is OK, since pointMap_ won't go away
  mvData_ (mv_.get1dViewNonConst ().getRawPtr ()),
  blockSize_ (blockSize)
{
  // Make sure that mv_ has view semantics.
  mv_.setCopyOrView (Teuchos::View);
}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const map_type& meshMap,
                  const map_type& pointMap,
                  const LO blockSize,
                  const LO numVecs) :
  dist_object_type (Teuchos::rcp (new map_type (meshMap))), // shallow copy
  meshMap_ (meshMap),
  pointMap_ (pointMap),
  mv_ (Teuchos::rcpFromRef (pointMap_), numVecs),
  mvData_ (mv_.get1dViewNonConst ().getRawPtr ()),
  blockSize_ (blockSize)
{
  // Make sure that mv_ has view semantics.
  mv_.setCopyOrView (Teuchos::View);
}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector (const mv_type& X_mv,
                  const map_type& meshMap,
                  const LO blockSize) :
  dist_object_type (Teuchos::rcp (new map_type (meshMap))), // shallow copy
  meshMap_ (meshMap),
  mvData_ (NULL), // just for now
  blockSize_ (blockSize)
{
  using Teuchos::RCP;

  if (X_mv.getCopyOrView () == Teuchos::View) {
    // The input MultiVector has view semantics, so assignment just
    // does a shallow copy.
    mv_ = X_mv;
  }
  else if (X_mv.getCopyOrView () == Teuchos::Copy) {
    // The input MultiVector has copy semantics.  We can't change
    // that, but we can make a view of the input MultiVector and
    // change the view to have view semantics.  (That sounds silly;
    // shouldn't views always have view semantics? but remember that
    // "view semantics" just governs the default behavior of the copy
    // constructor and assignment operator.)
    RCP<const mv_type> X_view_const;
    const size_t numCols = X_mv.getNumVectors ();
    if (numCols == 0) {
      Teuchos::Array<size_t> cols (0); // view will have zero columns
      X_view_const = X_mv.subView (cols ());
    } else { // Range1D is an inclusive range
      X_view_const = X_mv.subView (Teuchos::Range1D (0, numCols-1));
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      X_view_const.is_null (), std::logic_error, "Tpetra::Experimental::"
      "BlockMultiVector constructor: X_mv.subView(...) returned null.  This "
      "should never happen.  Please report this bug to the Tpetra developers.");

    // It's perfectly OK to cast away const here.  Those view methods
    // should be marked const anyway, because views can't change the
    // allocation (just the entries themselves).
    RCP<mv_type> X_view = Teuchos::rcp_const_cast<mv_type> (X_view_const);
    X_view->setCopyOrView (Teuchos::View);
    TEUCHOS_TEST_FOR_EXCEPTION(
      X_view->getCopyOrView () != Teuchos::View, std::logic_error, "Tpetra::"
      "Experimental::BlockMultiVector constructor: We just set a MultiVector "
      "to have view semantics, but it claims that it doesn't have view "
      "semantics.  This should never happen.  "
      "Please report this bug to the Tpetra developers.");
    mv_ = *X_view; // MultiVector::operator= does a shallow copy here
  }

  // At this point, mv_ has been assigned, so we can ignore X_mv.
  Teuchos::RCP<const map_type> pointMap = mv_.getMap ();
  if (! pointMap.is_null ()) {
    pointMap_ = *pointMap; // Map::operator= also does a shallow copy
  }
  mvData_ = mv_.get1dViewNonConst ().getRawPtr ();
}

template<class Scalar, class LO, class GO, class Node>
BlockMultiVector<Scalar, LO, GO, Node>::
BlockMultiVector () :
  dist_object_type (Teuchos::null),
  mvData_ (NULL),
  blockSize_ (0)
{
  // Make sure that mv_ has view semantics.
  mv_.setCopyOrView (Teuchos::View);
}

template<class Scalar, class LO, class GO, class Node>
typename BlockMultiVector<Scalar, LO, GO, Node>::map_type
BlockMultiVector<Scalar, LO, GO, Node>::
makePointMap (const map_type& meshMap, const LO blockSize)
{
  typedef Tpetra::global_size_t GST;
  typedef typename Teuchos::ArrayView<const GO>::size_type size_type;

  const GST gblNumMeshMapInds =
    static_cast<GST> (meshMap.getGlobalNumElements ());
  const size_t lclNumMeshMapIndices =
    static_cast<size_t> (meshMap.getNodeNumElements ());
  const GST gblNumPointMapInds =
    gblNumMeshMapInds * static_cast<GST> (blockSize);
  const size_t lclNumPointMapInds =
    lclNumMeshMapIndices * static_cast<size_t> (blockSize);
  const GO indexBase = meshMap.getIndexBase ();

  if (meshMap.isContiguous ()) {
    return map_type (gblNumPointMapInds, lclNumPointMapInds, indexBase,
                     meshMap.getComm (), meshMap.getNode ());
  }
  else {
    // "Hilbert's Hotel" trick: multiply each process' GIDs by
    // blockSize, and fill in.  That ensures correctness even if the
    // mesh Map is overlapping.
    Teuchos::ArrayView<const GO> lclMeshGblInds = meshMap.getNodeElementList ();
    const size_type lclNumMeshGblInds = lclMeshGblInds.size ();
    Teuchos::Array<GO> lclPointGblInds (lclNumPointMapInds);
    for (size_type g = 0; g < lclNumMeshGblInds; ++g) {
      const GO meshGid = lclMeshGblInds[g];
      const GO pointGidStart = indexBase + (meshGid - indexBase) * static_cast<GO> (blockSize);
      const size_type offset = g * static_cast<size_type> (blockSize);
      for (LO k = 0; k < blockSize; ++k) {
        const GO pointGid = pointGidStart + static_cast<GO> (k);
        lclPointGblInds[offset + static_cast<size_type> (k)] = pointGid;
      }
    }

    return map_type (gblNumPointMapInds, lclPointGblInds (), indexBase,
                     meshMap.getComm (), meshMap.getNode ());
  }
}


template<class Scalar, class LO, class GO, class Node>
void
BlockMultiVector<Scalar, LO, GO, Node>::
replaceLocalValuesImpl (const LO localRowIndex,
                        const LO colIndex,
                        const Scalar vals[]) const
{
  little_vec_type X_dst = getLocalBlock (localRowIndex, colIndex);
  const LO strideX = 1;
  const_little_vec_type X_src (vals, getBlockSize (), strideX);
  X_dst.assign (X_src);
}


template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
replaceLocalValues (const LO localRowIndex,
                    const LO colIndex,
                    const Scalar vals[]) const
{
  if (! meshMap_.isNodeLocalElement (localRowIndex)) {
    return false;
  } else {
    replaceLocalValuesImpl (localRowIndex, colIndex, vals);
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
replaceGlobalValues (const GO globalRowIndex,
                     const LO colIndex,
                     const Scalar vals[]) const
{
  const LO localRowIndex = meshMap_.getLocalElement (globalRowIndex);
  if (localRowIndex == Teuchos::OrdinalTraits<LO>::invalid ()) {
    return false;
  } else {
    replaceLocalValuesImpl (localRowIndex, colIndex, vals);
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
void
BlockMultiVector<Scalar, LO, GO, Node>::
sumIntoLocalValuesImpl (const LO localRowIndex,
                        const LO colIndex,
                        const Scalar vals[]) const
{
  little_vec_type X_dst = getLocalBlock (localRowIndex, colIndex);
  const LO strideX = 1;
  const_little_vec_type X_src (vals, getBlockSize (), strideX);

  X_dst.update (STS::one (), X_src);
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
sumIntoLocalValues (const LO localRowIndex,
                    const LO colIndex,
                    const Scalar vals[]) const
{
  if (! meshMap_.isNodeLocalElement (localRowIndex)) {
    return false;
  } else {
    sumIntoLocalValuesImpl (localRowIndex, colIndex, vals);
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
sumIntoGlobalValues (const GO globalRowIndex, const LO colIndex, const Scalar vals[]) const {
  const LO localRowIndex = meshMap_.getLocalElement (globalRowIndex);
  if (localRowIndex == Teuchos::OrdinalTraits<LO>::invalid ()) {
    return false;
  } else {
    sumIntoLocalValuesImpl (localRowIndex, colIndex, vals);
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
getLocalRowView (const LO localRowIndex, const LO colIndex, Scalar*& vals) const
{
  if (! meshMap_.isNodeLocalElement (localRowIndex)) {
    return false;
  } else {
    little_vec_type X_ij = getLocalBlock (localRowIndex, colIndex);
    vals = X_ij.getRawPtr ();
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
bool
BlockMultiVector<Scalar, LO, GO, Node>::
getGlobalRowView (const GO globalRowIndex, const LO colIndex, Scalar*& vals) const
{
  const LO localRowIndex = meshMap_.getLocalElement (globalRowIndex);
  if (localRowIndex == Teuchos::OrdinalTraits<LO>::invalid ()) {
    return false;
  } else {
    little_vec_type X_ij = getLocalBlock (localRowIndex, colIndex);
    vals = X_ij.getRawPtr ();
    return true;
  }
}

template<class Scalar, class LO, class GO, class Node>
typename BlockMultiVector<Scalar, LO, GO, Node>::little_vec_type
BlockMultiVector<Scalar, LO, GO, Node>::
getLocalBlock (const LO localRowIndex,
               const LO colIndex) const
{
  if (! isValidLocalMeshIndex (localRowIndex)) {
    return little_vec_type (NULL, 0, 0);
  } else {
    const LO strideX = this->getStrideX ();
    const size_t blockSize = getBlockSize ();
    const size_t offset = colIndex * this->getStrideY () +
      localRowIndex * blockSize * strideX;
    return little_vec_type (this->getRawPtr () + offset, blockSize, strideX);
  }
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const typename BlockMultiVector<Scalar, LO, GO, Node>::mv_type>
BlockMultiVector<Scalar, LO, GO, Node>::
getMultiVectorFromSrcDistObject (const Tpetra::SrcDistObject& src)
{
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  // The source object of an Import or Export must be a
  // BlockMultiVector or MultiVector (a Vector is a MultiVector).  Try
  // them in that order; one must succeed.  Note that the target of
  // the Import or Export calls checkSizes in DistObject's doTransfer.
  typedef BlockMultiVector<Scalar, LO, GO, Node> this_type;
  const this_type* srcBlkVec = dynamic_cast<const this_type*> (&src);
  if (srcBlkVec == NULL) {
    const mv_type* srcMultiVec = dynamic_cast<const mv_type*> (&src);
    if (srcMultiVec == NULL) {
      // FIXME (mfh 05 May 2014) Tpetra::MultiVector's operator=
      // currently does a shallow copy by default.  This is why we
      // return by RCP, rather than by value.
      return rcp (new mv_type ());
    } else { // src is a MultiVector
      return rcp (srcMultiVec, false); // nonowning RCP
    }
  } else { // src is a BlockMultiVector
    return rcpFromRef (srcBlkVec->mv_); // nonowning RCP
  }
}

template<class Scalar, class LO, class GO, class Node>
bool BlockMultiVector<Scalar, LO, GO, Node>::
checkSizes (const Tpetra::SrcDistObject& src)
{
  return ! getMultiVectorFromSrcDistObject (src).is_null ();
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
copyAndPermute (const Tpetra::SrcDistObject& src,
                size_t numSameIDs,
                const Teuchos::ArrayView<const LO>& permuteToLIDs,
                const Teuchos::ArrayView<const LO>& permuteFromLIDs)
{
  const char tfecfFuncName[] = "copyAndPermute";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
    ": permuteToLIDs and permuteFromLIDs must have the same size."
    << std::endl << "permuteToLIDs.size() = " << permuteToLIDs.size ()
    << " != permuteFromLIDs.size() = " << permuteFromLIDs.size () << ".");

  typedef BlockMultiVector<Scalar, LO, GO, Node> BMV;
  Teuchos::RCP<const BMV> srcAsBmvPtr = getBlockMultiVectorFromSrcDistObject (src);
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    srcAsBmvPtr.is_null (), std::invalid_argument,
    ": The source of an Import or Export to a BlockMultiVector "
    "must also be a BlockMultiVector.");
  const BMV& srcAsBmv = *srcAsBmvPtr;

  // FIXME (mfh 23 Apr 2014) This implementation is sequential and
  // assumes UVM.

  const LO numVecs = getNumVectors ();
  const LO numSame = static_cast<LO> (numSameIDs);
  for (LO j = 0; j < numVecs; ++j) {
    for (LO lclRow = 0; lclRow < numSame; ++lclRow) {
      getLocalBlock (lclRow, j).assign (srcAsBmv.getLocalBlock (lclRow, j));
    }
  }

  // FIXME (mfh 20 June 2012) For an Export with duplicate GIDs on the
  // same process, this merges their values by replacement of the last
  // encountered GID, not by the specified merge rule (such as ADD).
  const LO numPermuteLIDs = static_cast<LO> (permuteToLIDs.size ());
  for (LO j = 0; j < numVecs; ++j) {
    for (LO k = numSame; k < numPermuteLIDs; ++k) {
      getLocalBlock (permuteToLIDs[k], j).assign (srcAsBmv.getLocalBlock (permuteFromLIDs[k], j));
    }
  }
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
packAndPrepare (const Tpetra::SrcDistObject& src,
                const Teuchos::ArrayView<const LO>& exportLIDs,
                Teuchos::Array<packet_type>& exports,
                const Teuchos::ArrayView<size_t>& /* numPacketsPerLID */,
                size_t& constantNumPackets,
                Tpetra::Distributor& /* distor */)
{
  typedef BlockMultiVector<Scalar, LO, GO, Node> BMV;
  typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
  const char tfecfFuncName[] = "packAndPrepare";

  Teuchos::RCP<const BMV> srcAsBmvPtr = getBlockMultiVectorFromSrcDistObject (src);
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    srcAsBmvPtr.is_null (), std::invalid_argument,
    ": The source of an Import or Export to a BlockMultiVector "
    "must also be a BlockMultiVector.");
  const BMV& srcAsBmv = *srcAsBmvPtr;

  const LO numVecs = getNumVectors ();
  const LO blockSize = getBlockSize ();

  // Number of things to pack per LID is the block size, times the
  // number of columns.  Input LIDs correspond to the mesh points, not
  // the degrees of freedom (DOFs).
  constantNumPackets =
    static_cast<size_t> (blockSize) * static_cast<size_t> (numVecs);
  const size_type numMeshLIDs = exportLIDs.size ();

  const size_type requiredExportsSize = numMeshLIDs *
    static_cast<size_type> (blockSize) * static_cast<size_type> (numVecs);
  exports.resize (requiredExportsSize);

  try {
    size_type curExportPos = 0;
    for (size_type meshLidIndex = 0; meshLidIndex < numMeshLIDs; ++meshLidIndex) {
      for (LO j = 0; j < numVecs; ++j, curExportPos += blockSize) {
        const LO meshLid = exportLIDs[meshLidIndex];
        Scalar* const curExportPtr = &exports[curExportPos];
        little_vec_type X_dst (curExportPtr, blockSize, 1);
        little_vec_type X_src = srcAsBmv.getLocalBlock (meshLid, j);

        X_dst.assign (X_src);
      }
    }
  } catch (std::exception& e) {
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      true, std::logic_error, ": Oh no!  packAndPrepare on Process "
      << meshMap_.getComm ()->getRank () << " raised the following exception: "
      << e.what ());
  }
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
unpackAndCombine (const Teuchos::ArrayView<const LO>& importLIDs,
                  const Teuchos::ArrayView<const packet_type>& imports,
                  const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                  size_t constantNumPackets,
                  Tpetra::Distributor& distor,
                  Tpetra::CombineMode CM)
{
  typedef typename Teuchos::ArrayView<const LO>::size_type size_type;
  const char tfecfFuncName[] = "unpackAndCombine";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    CM != ADD && CM != REPLACE && CM != INSERT && CM != ABSMAX && CM != ZERO,
    std::invalid_argument, ": Invalid CombineMode: " << CM << ".  Valid "
    "CombineMode values are ADD, REPLACE, INSERT, ABSMAX, and ZERO.");

  if (CM == ZERO) {
    return; // Combining does nothing, so we don't have to combine anything.
  }

  // Number of things to pack per LID is the block size.
  // Input LIDs correspond to the mesh points, not the DOFs.
  const size_type numMeshLIDs = importLIDs.size ();
  const LO blockSize = getBlockSize ();
  const LO numVecs = getNumVectors ();

  const size_type requiredImportsSize = numMeshLIDs *
    static_cast<size_type> (blockSize) * static_cast<size_type> (numVecs);
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    imports.size () < requiredImportsSize, std::logic_error,
    ": imports.size () = " << imports.size ()
    << " < requiredImportsSize = " << requiredImportsSize << ".");

  size_type curImportPos = 0;
  for (size_type meshLidIndex = 0; meshLidIndex < numMeshLIDs; ++meshLidIndex) {
    for (LO j = 0; j < numVecs; ++j, curImportPos += blockSize) {
      const LO meshLid = importLIDs[meshLidIndex];
      const Scalar* const curImportPtr = &imports[curImportPos];

      const_little_vec_type X_src (curImportPtr, blockSize, 1);
      little_vec_type X_dst = getLocalBlock (meshLid, j);

      if (CM == INSERT || CM == REPLACE) {
        X_dst.assign (X_src);
      } else if (CM == ADD) {
        X_dst.update (STS::one (), X_src);
      } else if (CM == ABSMAX) {
        X_dst.absmax (X_src);
      }
    }
  }
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
putScalar (const Scalar& val)
{
  getMultiVectorView ().putScalar (val);
}

template<class Scalar, class LO, class GO, class Node>
void BlockMultiVector<Scalar, LO, GO, Node>::
scale (const Scalar& val)
{
  getMultiVectorView ().scale (val);
}

} // namespace Experimental
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_EXPERIMENTAL_BLOCKMULTIVECTOR_INSTANT(S,LO,GO,NODE) \
  template class Experimental::BlockMultiVector< S, LO, GO, NODE >;

#endif // TPETRA_EXPERIMENTAL_BLOCKMULTIVECTOR_DEF_HPP
