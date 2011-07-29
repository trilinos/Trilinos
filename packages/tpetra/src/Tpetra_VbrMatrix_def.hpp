//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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

#ifndef TPETRA_VBRMATRIX_DEF_HPP
#define TPETRA_VBRMATRIX_DEF_HPP

#include <algorithm>
#include <sstream>

#include <Kokkos_NodeHelpers.hpp>
#include <Tpetra_Vector.hpp>

#ifdef DOXYGEN_USE_ONLY
#include "Tpetra_VbrMatrix_decl.hpp"
#endif

namespace Tpetra {

template<class Scalar,class Node>
void fill_device_ArrayRCP(Teuchos::RCP<Node>& node, Teuchos::ArrayRCP<Scalar>& ptr, Scalar value)
{
  Kokkos::ReadyBufferHelper<Node> rbh(node);
  Kokkos::InitOp<Scalar> wdp;
  wdp.alpha = value;
  rbh.begin();
  wdp.x = rbh.addNonConstBuffer(ptr);
  rbh.end();
  node->template parallel_for<Kokkos::InitOp<Scalar> >(0, ptr.size(), wdp);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::VbrMatrix(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype)
 : DistObject<char, LocalOrdinal, GlobalOrdinal, Node>(convertBlockMapToPointMap(*blkRowMap)),
   blkGraph_(Teuchos::rcp(new BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>(blkRowMap, maxNumEntriesPerRow, pftype))),
   constBlkGraph_(blkGraph_),
   lclMatrix_(blkRowMap->getNodeNumBlocks(), blkRowMap->getPointMap()->getNode()),
   pbuf_values1D_(),
   pbuf_indx_(),
   lclMatOps_(blkRowMap->getPointMap()->getNode()),
   importer_(),
   exporter_(),
   importedVec_(),
   exportedVec_(),
   data_2D_(Teuchos::rcp(new Teuchos::Array<RowGlobalCols>)),
   nonlocal_data_(),
   is_fill_completed_(false),
   is_storage_optimized_(false)
{
  //The graph of this VBR matrix is a BlockCrsGraph, which is a CrsGraph where
  //each entry in the graph corresponds to a block-entry in the matrix.
  //That is, you can think of a VBR matrix as a Crs matrix of dense
  //submatrices...
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::VbrMatrix(const Teuchos::RCP<const BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node> > &blkGraph)
 : DistObject<char, LocalOrdinal, GlobalOrdinal, Node>(convertBlockMapToPointMap(*blkGraph->getBlockRowMap())),
   blkGraph_(),
   constBlkGraph_(blkGraph),
   lclMatrix_(blkGraph->getBlockRowMap()->getNodeNumBlocks(),
              blkGraph->getBlockRowMap()->getPointMap()->getNode()),
   pbuf_values1D_(),
   pbuf_indx_(),
   lclMatOps_(blkGraph->getBlockRowMap()->getPointMap()->getNode()),
   importer_(),
   exporter_(),
   importedVec_(),
   exportedVec_(),
   data_2D_(Teuchos::rcp(new Teuchos::Array<RowGlobalCols>)),
   nonlocal_data_(),
   is_fill_completed_(false),
   is_storage_optimized_(false)
{
  //The graph of this VBR matrix is a BlockCrsGraph, which is a CrsGraph where
  //each entry in the graph corresponds to a block-entry in the matrix.
  //That is, you can think of a VBR matrix as a Crs matrix of dense
  //submatrices...

  TEST_FOR_EXCEPTION(blkGraph->isFillComplete() == false, std::runtime_error,
   "Tpetra::VbrMatrix::VbrMatrix(BlockCrsGraph) ERROR, this constructor requires graph.isFillComplete()==true.");

  createImporterExporter();

  is_fill_completed_ = true;

  optimizeStorage();

  fillLocalMatrix();
  fillLocalMatVec();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::~VbrMatrix()
{
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getDomainMap() const
{
  return getBlockDomainMap()->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRangeMap() const
{
  return getBlockRangeMap()->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
bool
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isFillComplete() const
{
  return is_fill_completed_;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
template <class DomainScalar, class RangeScalar>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const
{
  TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
    "Tpetra::VbrMatrix::multiply ERROR, multiply may only be called after fillComplete has been called.");

  const Kokkos::MultiVector<Scalar,Node> *lclX = &X.getLocalMV();
  Kokkos::MultiVector<Scalar,Node>        *lclY = &Y.getLocalMVNonConst();

  lclMatOps_.template multiply<Scalar,Scalar>(trans,alpha,*lclX,beta,*lclY);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
template<class DomainScalar, class RangeScalar>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>& Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node>& X, Teuchos::ETransp trans) const
{
  TEST_FOR_EXCEPTION(X.isConstantStride() == false || Y.isConstantStride() == false, std::runtime_error,
        "Tpetra::VbrMatrix::solve(X,Y): X and Y must be constant stride.");

  TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
    "Tpetra::VbrMatrix::solve ERROR, solve may only be called after fillComplete has been called.");

  TEST_FOR_EXCEPTION(constBlkGraph_->isUpperTriangular()==false && constBlkGraph_->isLowerTriangular()==false, std::runtime_error,
    "Tpetra::VbrMatrix::solve ERROR, matrix must be either upper or lower triangular.");

  const Kokkos::MultiVector<RangeScalar,Node> *lclY = &Y.getLocalMV();
  Kokkos::MultiVector<DomainScalar,Node>      *lclX = &X.getLocalMVNonConst();

  Teuchos::EUplo triang = constBlkGraph_->isUpperTriangular() ? Teuchos::UPPER_TRI : Teuchos::LOWER_TRI;
  Teuchos::EDiag diag = (constBlkGraph_->getNodeNumBlockDiags() < constBlkGraph_->getNodeNumBlockRows()) ? Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG;

  lclMatOps_.template solve<DomainScalar,RangeScalar>(trans,triang,diag, *lclY,*lclX);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::updateImport(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X) const
{
  if (importer_ != Teuchos::null) {
    if (importedVec_ == Teuchos::null || importedVec_->getNumVectors() != X.getNumVectors()) {
      importedVec_ = Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(getBlockColMap()->getPointMap(), X.getNumVectors()));
    }

    importedVec_->doImport(X, *importer_, REPLACE);
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::updateExport(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  if (exporter_ != Teuchos::null) {
    if (exportedVec_ == Teuchos::null || exportedVec_->getNumVectors() != Y.getNumVectors()) {
      exportedVec_ = Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(getBlockColMap()->getPointMap(), Y.getNumVectors()));
    }

    exportedVec_->doExport(Y, *exporter_, REPLACE);
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::apply(
         const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
               Teuchos::ETransp trans,
               Scalar alpha,
               Scalar beta) const
{
  if (trans == Teuchos::NO_TRANS) {
    updateImport(X);
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xref = importedVec_ == Teuchos::null ? X : *importedVec_;
    this->template multiply<Scalar,Scalar>(Xref, Y, trans, alpha, beta);
  }
  else if (trans == Teuchos::TRANS || trans == Teuchos::CONJ_TRANS) {
    updateImport(Y);
    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Yref = importedVec_ == Teuchos::null ? Y : *importedVec_;
    this->template multiply<Scalar,Scalar>(X, Yref, trans, alpha, beta);
    if (importedVec_ != Teuchos::null) {
      Y.doExport(*importedVec_, *importer_, ADD);
    }
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::applyInverse(
         const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
               Teuchos::ETransp trans) const
{
  if (trans == Teuchos::NO_TRANS) {
    updateImport(Y);
    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xref = importedVec_ == Teuchos::null ? X : *importedVec_;
    this->template solve<Scalar,Scalar>(Y, Xref, trans);
  }
  else if (trans == Teuchos::TRANS || trans == Teuchos::CONJ_TRANS) {
    updateImport(Y);
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Yref = importedVec_ == Teuchos::null ? Y : *importedVec_;
    this->template solve<Scalar,Scalar>(Yref, X, trans);
    if (importedVec_ != Teuchos::null) {
      X.doExport(*importedVec_, *importer_, ADD);
    }
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
bool
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasTransposeApply() const
{
  return true;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getBlockRowMap() const
{
  return constBlkGraph_->getBlockRowMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getPointRowMap() const
{
  return constBlkGraph_->getBlockRowMap()->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getBlockColMap() const
{
  return constBlkGraph_->getBlockColMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getBlockDomainMap() const
{
  return constBlkGraph_->getBlockDomainMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getBlockRangeMap() const
{
  return constBlkGraph_->getBlockRangeMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getPointColMap() const
{
  return constBlkGraph_->getBlockColMap()->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalBlockRowView(
         GlobalOrdinal globalBlockRow,
         LocalOrdinal& numPtRows,
         Teuchos::ArrayView<const GlobalOrdinal>& blockCols,
         Teuchos::Array<LocalOrdinal>& ptColsPerBlockCol,
         Teuchos::Array<Teuchos::ArrayRCP<const Scalar> >& blockEntries) const
{
  TEST_FOR_EXCEPTION(isFillComplete(), std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockRowView internal ERROR, isFillComplete() is required to be false.");

  typedef typename Teuchos::ArrayView<const GlobalOrdinal>::size_type Tsize_t;

  LocalOrdinal localRow = getBlockRowMap()->getLocalBlockID(globalBlockRow);
  numPtRows = getBlockRowMap()->getLocalBlockSize(localRow);
  constBlkGraph_->getGlobalBlockRowView(globalBlockRow, blockCols);
  ptColsPerBlockCol.resize(blockCols.size());
  blockEntries.resize(blockCols.size());
  for(Tsize_t i=0; i<blockCols.size(); ++i) {
    getGlobalBlockEntryView(globalBlockRow, blockCols[i],
                           numPtRows, ptColsPerBlockCol[i], blockEntries[i]);
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalBlockRowView(
         LocalOrdinal localBlockRow,
         LocalOrdinal& numPtRows,
         Teuchos::ArrayView<const LocalOrdinal>& blockCols,
         Teuchos::Array<LocalOrdinal>& ptColsPerBlockCol,
         Teuchos::ArrayRCP<const Scalar>& blockEntries) const
{
  TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockRowView internal ERROR, isFillComplete() is required to be true.");

  typedef typename Teuchos::ArrayView<const GlobalOrdinal>::size_type Tsize_t;
  typedef Teuchos::ArrayRCP<const size_t>       Host_View;
  typedef Teuchos::ArrayRCP<const LocalOrdinal> Host_View_LO;

  numPtRows = getBlockRowMap()->getLocalBlockSize(localBlockRow);
  constBlkGraph_->getLocalBlockRowView(localBlockRow, blockCols);
  ptColsPerBlockCol.resize(blockCols.size());
  LocalOrdinal num_scalars_in_block_row = 0;
  for(Tsize_t i=0; i<blockCols.size(); ++i) {
    ptColsPerBlockCol[i] = getBlockColMap()->getLocalBlockSize(blockCols[i]);
    num_scalars_in_block_row += numPtRows*ptColsPerBlockCol[i];
  }

  Teuchos::RCP<Node> node = getNode();
  Host_View bpr = constBlkGraph_->getNodeRowOffsets();
  const size_t bindx_offset = bpr[localBlockRow];
  Host_View_LO indx = node->template viewBuffer<LocalOrdinal>(1, pbuf_indx_+bindx_offset);
  const LocalOrdinal offset = indx[0];
  blockEntries = node->template viewBuffer<Scalar>(num_scalars_in_block_row, pbuf_values1D_+offset);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalBlockEntryViewNonConst(
    GlobalOrdinal globalBlockRow,
    GlobalOrdinal globalBlockCol,
    LocalOrdinal& numPtRows,
    LocalOrdinal& numPtCols,
    Teuchos::ArrayRCP<Scalar>& blockEntry)
{
  //Return a non-const, read-write view of a block-entry (as an ArrayRCP),
  //creating/allocating the block-entry if it doesn't already exist, (but only
  //if fillComplete hasn't been called yet, and if the arguments numPtRows
  //and numPtCols have been set on input).

  LocalOrdinal localBlockRow = getBlockRowMap()->getLocalBlockID(globalBlockRow);

  if (localBlockRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid()) {
    VbrUtils::getGlobalBlockEntryViewNonConst(nonlocal_data_,
                   globalBlockRow, globalBlockCol,
                   numPtRows, numPtCols, blockEntry);
    return;
  }

  if (is_storage_optimized_) {
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockEntryViewNonConst internal ERROR, storage is optimized but isFillComplete() is false.");

    LocalOrdinal localBlockCol = getBlockColMap()->getLocalBlockID(globalBlockCol);

    getLocalBlockEntryViewNonConst(localBlockRow, localBlockCol,
                                   numPtRows, numPtCols, blockEntry);
    return;
  }

  //If we get to here, fillComplete hasn't been called yet, and the matrix data
  //is stored in un-packed '2D' form.

  if (data_2D_->size() == 0) {
    LocalOrdinal numBlockRows = constBlkGraph_->getNodeNumBlockRows();
    data_2D_->resize(numBlockRows);
  }

  RowGlobalCols& blkrow = (*data_2D_)[localBlockRow];

  typename RowGlobalCols::iterator col_iter = blkrow.find(globalBlockCol);

  if (col_iter != blkrow.end()) {
    blockEntry = col_iter->second;
  }
  else {
    //blockEntry doesn't already exist, so we will create it.

    //make sure block-size is specified:
    TEST_FOR_EXCEPTION(numPtRows==0 || numPtCols==0, std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockEntryViewNonConst ERROR: creating block-entry, but numPtRows and/or numPtCols is 0.");

    Teuchos::RCP<Node> node = getNode();
    size_t blockSize = numPtRows*numPtCols;
    blockEntry = Teuchos::arcp(new Scalar[blockSize], 0, blockSize);
    std::fill(blockEntry.begin(), blockEntry.end(), 0);
    blkrow.insert(std::make_pair(globalBlockCol, blockEntry));
    blkGraph_->insertGlobalIndices(globalBlockRow, Teuchos::ArrayView<GlobalOrdinal>(&globalBlockCol, 1));
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalBlockEntryView(
      GlobalOrdinal globalBlockRow,
      GlobalOrdinal globalBlockCol,
      LocalOrdinal& numPtRows,
      LocalOrdinal& numPtCols,
      Teuchos::ArrayRCP<const Scalar>& blockEntry) const
{
  //This method returns a const-view of a block-entry (as an ArrayRCP).
  //Throws an exception if the block-entry doesn't already exist.

  if (is_storage_optimized_) {
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockEntryView internal ERROR, storage is optimized but isFillComplete() is false.");

    LocalOrdinal localBlockRow = getBlockRowMap()->getLocalBlockID(globalBlockRow);
    LocalOrdinal localBlockCol = getBlockColMap()->getLocalBlockID(globalBlockCol);
    getLocalBlockEntryView(localBlockRow, localBlockCol, numPtRows, numPtCols, blockEntry);
    return;
  }

  TEST_FOR_EXCEPTION(data_2D_->size() == 0, std::runtime_error,
    "Tpetra::VbrMatrix::getGlobalBlockEntryView ERROR, matrix storage not yet allocated, can't return a const view.");

  //this acts as a range-check for globalBlockRow:
  LocalOrdinal localBlockRow = getBlockRowMap()->getLocalBlockID(globalBlockRow);
  TEST_FOR_EXCEPTION( localBlockRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
     std::runtime_error,
     "Tpetra::VbrMatrix::getGlobalBlockEntryView, globalBlockRow not on the local processor.");

  RowGlobalCols& blkrow = (*data_2D_)[localBlockRow];
  typename RowGlobalCols::iterator col_iter = blkrow.find(globalBlockCol);

  if (col_iter == blkrow.end()) {
    throw std::runtime_error("Tpetra::VbrMatrix::getGlobalBlockEntryView ERROR, specified block-entry doesn't exist.");
  }

  numPtRows = getBlockRowMap()->getLocalBlockSize(localBlockRow);

  TEST_FOR_EXCEPTION(numPtRows == 0, std::runtime_error,
    "Tpetra::VbrMatrix::getGlobalBlockEntryView ERROR, numPtRows == 0.");

  blockEntry = col_iter->second;
  numPtCols = blockEntry.size() / numPtRows;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalBlockEntryViewNonConst(
    LocalOrdinal localBlockRow,
    LocalOrdinal localBlockCol,
    LocalOrdinal& numPtRows,                     
    LocalOrdinal& numPtCols,
    Teuchos::ArrayRCP<Scalar>& blockEntry)
{
  typedef Teuchos::ArrayRCP<const LocalOrdinal> Host_View_LO;
  typedef Teuchos::ArrayRCP<const size_t> Host_View;
  typedef typename Host_View_LO::iterator ITER;
  //This method returns a non-constant view of a block-entry (as an ArrayRCP).

  TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
   "Tpetra::VbrMatrix::getLocalBlockEntryViewNonConst ERROR, this method can only be called after fillComplete() has been called.");

  TEST_FOR_EXCEPTION(is_storage_optimized_ == false, std::runtime_error,
   "Tpetra::VbrMatrix::getLocalBlockEntryViewNonConst ERROR, this method can only be called if storage is optimized.");

  Teuchos::RCP<Node> node = getNode();

  Host_View bptr = constBlkGraph_->getNodeRowOffsets();
  LocalOrdinal bindx_offset = bptr[localBlockRow];
  LocalOrdinal length = bptr[localBlockRow+1] - bindx_offset;

  TEST_FOR_EXCEPTION( length < 1, std::runtime_error,
    "Tpetra::VbrMatrix::getLocalBlockEntryViewNonConst ERROR, specified localBlockCol not found in localBlockRow.");

  Host_View_LO bindx = constBlkGraph_->getNodePackedIndices();
  ITER bindx_beg = bindx.begin() + bindx_offset,  
       bindx_end = bindx_beg + length;
  ITER it = std::lower_bound(bindx_beg, bindx_end, localBlockCol);

  TEST_FOR_EXCEPTION(it == bindx_end || *it != localBlockCol, std::runtime_error,
    "Tpetra::VbrMatrix::getLocalBlockEntryViewNonConst ERROR, specified localBlockCol not found.");

  numPtRows = getBlockRowMap()->getLocalBlockSize(localBlockRow);
  numPtCols = getBlockColMap()->getLocalBlockSize(localBlockCol);

  const LocalOrdinal blkSize = numPtRows*numPtCols;
  Host_View_LO indx = node->template viewBuffer<LocalOrdinal>(1,pbuf_indx_+bptr[localBlockRow]+(it-bindx_beg));
  const LocalOrdinal offset = indx[0];
  blockEntry = node->template viewBufferNonConst<Scalar>(Kokkos::ReadWrite, blkSize, pbuf_values1D_ + offset);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalDiagCopy(
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& diag) const
{
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& rowmap = getBlockRowMap();
  TEST_FOR_EXCEPTION(diag.getMap()->isSameAs(*(rowmap->getPointMap())) != true,
    std::runtime_error, "Tpetra::VbrMatrix::getLocalDiagCopy ERROR, vector must be distributed the same as this matrix' row-map.");

  Teuchos::ArrayRCP<Scalar> diag_view = diag.get1dViewNonConst();
  Teuchos::ArrayView<const GlobalOrdinal> blockIDs = rowmap->getNodeBlockIDs();
  size_t offset = 0;
  typedef typename Teuchos::ArrayView<const GlobalOrdinal>::size_type Tsize_t;
  for(Tsize_t i=0; i<blockIDs.size(); ++i) {
    LocalOrdinal localBlockID = rowmap->getLocalBlockID(blockIDs[i]);
    LocalOrdinal blockSize = rowmap->getLocalBlockSize(localBlockID);
    Teuchos::ArrayRCP<const Scalar> blockEntry;
    getGlobalBlockEntryView(blockIDs[i], blockIDs[i], blockSize, blockSize, blockEntry);

    for(LocalOrdinal j=0; j<blockSize; ++j) {
      diag_view[offset++] = blockEntry[j*blockSize+j];
    }
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
bool
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::checkSizes(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source)
{
  bool ok = this->getMap()->getMinAllGlobalIndex() <= source.getMap()->getMinAllGlobalIndex();
  ok = ok && this->getMap()->getMaxAllGlobalIndex() >= source.getMap()->getMaxAllGlobalIndex();
  return ok;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::copyAndPermute(
   const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source,
   size_t numSameIDs,
   const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
   const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs)
{
  typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type Tsize_t;

  const VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>* src_mat_ptr = dynamic_cast<const VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>*>(&source);
  if (src_mat_ptr == NULL) {
    throw std::runtime_error("VbrMatrix::copyAndPermute ERROR, dynamic_cast failed.");
  }
  const VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>& src_mat = *src_mat_ptr;

  Teuchos::Array<LocalOrdinal> src_ptColsPerBlkCol;
  LocalOrdinal numPtRows;

  LocalOrdinal numSame = numSameIDs;
  for(LocalOrdinal i=0; i<numSame; ++i) {
    if (isFillComplete()) {
      Teuchos::ArrayRCP<const Scalar> src_blkEntries;
      Teuchos::ArrayView<const LocalOrdinal> src_blkCols;

      src_mat.getLocalBlockRowView(i, numPtRows, src_blkCols, src_ptColsPerBlkCol, src_blkEntries);
      unsigned offset = 0;
      for(Tsize_t j=0; j<src_blkCols.size(); ++j) {
        LocalOrdinal numPtCols = src_ptColsPerBlkCol[j];
        setLocalBlockEntry(i, src_blkCols[j], numPtRows, numPtCols, numPtCols,
                           src_blkEntries.view(offset, numPtRows*numPtCols));
        offset += numPtRows*numPtCols;
      }
    }
    else {
      GlobalOrdinal gRow = getBlockRowMap()->getGlobalBlockID(i);
      Teuchos::Array<Teuchos::ArrayRCP<const Scalar> > src_blkEntries;
      Teuchos::ArrayView<const GlobalOrdinal> src_blkCols;

      src_mat.getGlobalBlockRowView(gRow, numPtRows, src_blkCols, src_ptColsPerBlkCol, src_blkEntries);
      for(Tsize_t j=0; j<src_blkCols.size(); ++j) {
        LocalOrdinal numPtCols = src_ptColsPerBlkCol[j];
        setGlobalBlockEntry(gRow, src_blkCols[j], numPtRows, numPtCols, numPtRows, src_blkEntries[j].view(0, numPtRows*numPtCols));
      }
    }
  }

  for(Tsize_t i=0; i<permuteToLIDs.size(); ++i) {

    if (isFillComplete()) {
      Teuchos::ArrayRCP<const Scalar> src_blkEntries;
      Teuchos::ArrayView<const LocalOrdinal> src_blkCols;

      src_mat.getLocalBlockRowView(permuteFromLIDs[i], numPtRows, src_blkCols, src_ptColsPerBlkCol, src_blkEntries);
      unsigned offset = 0;
      for(Tsize_t j=0; j<src_blkCols.size(); ++j) {
        LocalOrdinal numPtCols = src_ptColsPerBlkCol[j];
        setLocalBlockEntry(permuteToLIDs[i], src_blkCols[j], numPtRows, numPtCols, numPtCols,
                           src_blkEntries.view(offset, numPtRows*numPtCols));
        offset += numPtRows*numPtCols;
      }
    }
    else {
      GlobalOrdinal gRow = getBlockRowMap()->getGlobalBlockID(permuteToLIDs[i]);
      Teuchos::Array<Teuchos::ArrayRCP<const Scalar> > src_blkEntries;
      Teuchos::ArrayView<const GlobalOrdinal> src_blkCols;

      src_mat.getGlobalBlockRowView(gRow, numPtRows, src_blkCols, src_ptColsPerBlkCol, src_blkEntries);
      for(Tsize_t j=0; j<src_blkCols.size(); ++j) {
        LocalOrdinal numPtCols = src_ptColsPerBlkCol[j];
        setGlobalBlockEntry(gRow, src_blkCols[j], numPtRows, numPtCols, numPtRows, src_blkEntries[j].view(0, numPtRows*numPtCols));
      }
    }
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::packAndPrepare(
   const DistObject<char, LocalOrdinal, GlobalOrdinal, Node>& source,
   const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
   Teuchos::Array<char>& exports,
   const Teuchos::ArrayView<size_t>& numPacketsPerLID,
   size_t& constantNumPackets,
   Distributor& distor)
{
  const VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>* src_mat_ptr = dynamic_cast<const VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>*>(&source);
  if (src_mat_ptr == NULL) {
    typedef VbrUtils::VbrDataDist<LocalOrdinal,GlobalOrdinal,Scalar,Node> VDD;
    VDD* vdd = const_cast<VDD*>(dynamic_cast<const VDD*>(&source));
    if (vdd != NULL) {
      vdd->packAndPrepare(source, exportLIDs, exports, numPacketsPerLID, constantNumPackets, distor);
    }
    else {
      throw std::runtime_error("VbrMatrix::packAndPrepare ERROR, dynamic_cast failed.");
    }
    return;
  }

  const VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>& src_mat = *src_mat_ptr;

  typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type Tsize_t;

  //We will pack each row's data into the exports buffer as follows:
  //[num-block-cols,numPtRows,{list-of-blk-cols},{list-of-ptColsPerBlkCol},{all vals}]
  //so the length of the char exports buffer for a row is:
  //sizeof(LocalOrdinal)*(2+2*(num-block-cols)) + sizeof(Scalar)*numPtRows*sum(numPtCols_i)

  //For each row corresponding to exportLIDs, accumulate the size that it will
  //occupy in the exports buffer:
  size_t total_exports_size = 0;
  for(Tsize_t i=0; i<exportLIDs.size(); ++i) {
    Teuchos::Array<LocalOrdinal> src_ptColsPerBlkCol;
    LocalOrdinal numPtRows;

    LocalOrdinal numBlkCols = 0;
    LocalOrdinal numScalars = 0;

    if (src_mat.isFillComplete()) {
      Teuchos::ArrayRCP<const Scalar> src_blkEntries;
      Teuchos::ArrayView<const LocalOrdinal> src_blkCols;
      src_mat.getLocalBlockRowView(exportLIDs[i], numPtRows,
                            src_blkCols, src_ptColsPerBlkCol, src_blkEntries);
      numBlkCols = src_blkCols.size();
      numScalars = src_blkEntries.size();
    }
    else {//src_mat.isFillComplete()==false
      Teuchos::Array<Teuchos::ArrayRCP<const Scalar> > src_blkEntries;
      GlobalOrdinal gRow = src_mat.getBlockRowMap()->getGlobalBlockID(exportLIDs[i]);
      Teuchos::ArrayView<const GlobalOrdinal> src_blkCols;
      src_mat.getGlobalBlockRowView(gRow, numPtRows, src_blkCols, src_ptColsPerBlkCol, src_blkEntries);
      numBlkCols = src_blkCols.size();
      for(Tsize_t j=0; j<src_blkEntries.size(); ++j) {
        numScalars += src_blkEntries[j].size();
      }
    }

    size_t size_for_this_row = sizeof(GlobalOrdinal)*(2+2*numBlkCols)
                         + sizeof(Scalar)*numScalars;
    numPacketsPerLID[i] = size_for_this_row;
    total_exports_size += size_for_this_row;
  }

  exports.resize(total_exports_size);

  ArrayView<char> avIndsC, avValsC;
  ArrayView<Scalar>        avVals;

  size_t offset = 0;

  for(Tsize_t i=0; i<exportLIDs.size(); ++i) {
    Teuchos::Array<LocalOrdinal> src_ptColsPerBlkCol;
    LocalOrdinal numPtRows;

    LocalOrdinal numBlkCols = 0;
    LocalOrdinal numScalars = 0;

    if (src_mat.isFillComplete()) {
      Teuchos::ArrayRCP<const Scalar> src_blkEntries;
      Teuchos::ArrayView<const LocalOrdinal> src_blkCols;
      src_mat.getLocalBlockRowView(exportLIDs[i], numPtRows, src_blkCols, src_ptColsPerBlkCol, src_blkEntries);
      numBlkCols = src_blkCols.size();
      numScalars = src_blkEntries.size();

      LocalOrdinal num_chars_for_ordinals = (2*numBlkCols+2)*sizeof(LocalOrdinal);
      //get export views
      avIndsC = exports(offset, num_chars_for_ordinals);
      avValsC = exports(offset+ num_chars_for_ordinals, numScalars*sizeof(Scalar));
      ArrayView<GlobalOrdinal> avInds = av_reinterpret_cast<GlobalOrdinal>(avIndsC);
      typename ArrayView<GlobalOrdinal>::iterator ind_it = avInds.begin();

      const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& col_map = src_mat.getBlockColMap();

      //put row info into the buffer views:
      *ind_it++ = numBlkCols;
      *ind_it++ = numPtRows;
      for(Tsize_t j=0; j<src_blkCols.size(); ++j) {
        *ind_it++ = col_map->getGlobalBlockID(src_blkCols[j]);
      }
      std::copy(src_ptColsPerBlkCol.begin(), src_ptColsPerBlkCol.end(), ind_it);
      avVals = av_reinterpret_cast<Scalar>(avValsC);
      std::copy(src_blkEntries.begin(), src_blkEntries.end(), avVals.begin());
    }
    else {//src_mat.isFillComplete()==false
      Teuchos::Array<Teuchos::ArrayRCP<const Scalar> > src_blkEntries;
      GlobalOrdinal gRow = src_mat.getBlockRowMap()->getGlobalBlockID(exportLIDs[i]);
      Teuchos::ArrayView<const GlobalOrdinal> src_blkCols;
      src_mat.getGlobalBlockRowView(gRow, numPtRows, src_blkCols, src_ptColsPerBlkCol, src_blkEntries);
      numBlkCols = src_blkCols.size();
      numScalars = 0;
      for(Tsize_t j=0; j<src_blkEntries.size(); ++j) {
        numScalars += numPtRows*src_ptColsPerBlkCol[j];
      }

      LocalOrdinal num_chars_for_ordinals = (2*numBlkCols+2)*sizeof(GlobalOrdinal);
      //get export views
      avIndsC = exports(offset, num_chars_for_ordinals);
      avValsC = exports(offset+ num_chars_for_ordinals, numScalars*sizeof(Scalar));
      ArrayView<GlobalOrdinal> avInds = av_reinterpret_cast<GlobalOrdinal>(avIndsC);
      typename ArrayView<GlobalOrdinal>::iterator ind_it = avInds.begin();

      //put row info into the buffer views:
      *ind_it++ = numBlkCols;
      *ind_it++ = numPtRows;
      std::copy(src_blkCols.begin(), src_blkCols.end(), ind_it);
      ind_it += src_blkCols.size();
      std::copy(src_ptColsPerBlkCol.begin(), src_ptColsPerBlkCol.end(), ind_it);
      avVals = av_reinterpret_cast<Scalar>(avValsC);
      typename ArrayView<Scalar>::iterator val_it = avVals.begin();
      for(Tsize_t j=0; j<src_blkEntries.size(); ++j) {
        std::copy(src_blkEntries[j].begin(), src_blkEntries[j].end(), val_it);
        val_it += src_blkEntries[j].size();
      }
    }

    size_t size_for_this_row = sizeof(GlobalOrdinal)*(2+2*numBlkCols)
                         + sizeof(Scalar)*numScalars;
    offset += size_for_this_row;
  }

  constantNumPackets = 0;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::unpackAndCombine(
    const Teuchos::ArrayView<const LocalOrdinal>& importLIDs,
    const Teuchos::ArrayView<const char>& imports,
    const Teuchos::ArrayView<size_t>& numPacketsPerLID,
    size_t constantNumPackets,
    Distributor& distor,
    CombineMode CM)
{
  typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type Tsize_t;
 
  if (CM == Tpetra::ABSMAX) {
    std::cout << "Warning, VbrMatrix Import/Export doesn't support combine-mode==ABSMAX; use ADD, INSERT or REPLACE. (REPLACE is being used now.)"<<std::endl;
  }

  size_t offset = 0;
  for(Tsize_t i=0; i<importLIDs.size(); ++i) {
    ArrayView<const char> avC = imports.view(offset, numPacketsPerLID[i]);
    ArrayView<const GlobalOrdinal> avOrds = av_reinterpret_cast<const GlobalOrdinal>(avC);
    GlobalOrdinal gRow = this->getBlockRowMap()->getGlobalBlockID(importLIDs[i]);
    size_t avOrdsOffset = 0;
    LocalOrdinal numBlkCols = avOrds[avOrdsOffset++];
    LocalOrdinal numPtRows = avOrds[avOrdsOffset++];
    LocalOrdinal num_chars_for_ordinals = (2*numBlkCols+2)*sizeof(GlobalOrdinal);
    ArrayView<const char> avValsC = imports.view(offset+num_chars_for_ordinals, numPacketsPerLID[i]-num_chars_for_ordinals);
    ArrayView<const Scalar> avVals = av_reinterpret_cast<const Scalar>(avValsC);

    size_t avValsOffset = 0;
    for(Tsize_t j=0; j<numBlkCols; ++j) {
      GlobalOrdinal blkCol = avOrds[avOrdsOffset+j];
      LocalOrdinal numPtCols = avOrds[avOrdsOffset+numBlkCols+j];
      LocalOrdinal LDA = numPtRows;
      if (CM == Tpetra::ADD) {
        sumIntoGlobalBlockEntry(gRow, blkCol, numPtRows, numPtCols, LDA,
                         avVals.view(avValsOffset, numPtRows*numPtCols));
      }
      else if (CM == Tpetra::INSERT || CM == Tpetra::REPLACE) {
        setGlobalBlockEntry(gRow, blkCol, numPtRows, numPtCols, LDA,
                           avVals.view(avValsOffset, numPtRows*numPtCols));
      }
      else {
//    std::cout << "Warning, VbrMatrix Import/Export doesn't support combine-mode==ABSMAX; use ADD, INSERT or REPLACE. (REPLACE is being used now.)"<<std::endl;
        setGlobalBlockEntry(gRow, blkCol, numPtRows, numPtCols, LDA,
                           avVals.view(avValsOffset, numPtRows*numPtCols));
      }
      avValsOffset += numPtRows*numPtCols;
    }
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalBlockEntryView(
      LocalOrdinal localBlockRow,
      LocalOrdinal localBlockCol,
      LocalOrdinal& numPtRows,
      LocalOrdinal& numPtCols,
      Teuchos::ArrayRCP<const Scalar>& blockEntry) const
{
  typedef Teuchos::ArrayRCP<const size_t>       Host_View;
  typedef Teuchos::ArrayRCP<const LocalOrdinal> Host_View_LO;
  typedef typename Host_View_LO::iterator ITER;
  //This method returns a constant view of a block-entry (as an ArrayRCP).

  TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
   "Tpetra::VbrMatrix::getLocalBlockEntryView ERROR, this method can only be called after fillComplete() has been called.");

  TEST_FOR_EXCEPTION(is_storage_optimized_ == false, std::runtime_error,
   "Tpetra::VbrMatrix::getLocalBlockEntryView ERROR, this method can only be called if storage is optimized.");

  Teuchos::RCP<Node> node = getNode();

  Host_View bpr = constBlkGraph_->getNodeRowOffsets();
  const size_t bindx_offset = bpr[localBlockRow];
  const size_t length       = bpr[localBlockRow+1] - bindx_offset;

  Host_View_LO bindx = constBlkGraph_->getNodePackedIndices();
  ITER bindx_beg = bindx.begin()+bindx_offset, 
       bindx_end = bindx_beg + length;
  ITER it = std::lower_bound(bindx_beg, bindx_end, localBlockCol);

  TEST_FOR_EXCEPTION(it == bindx_end || *it != localBlockCol, std::runtime_error,
    "Tpetra::VbrMatrix::getLocalBlockEntryView ERROR, specified localBlockCol not found.");

  numPtRows = getBlockRowMap()->getLocalBlockSize(localBlockRow);
  numPtCols = getBlockColMap()->getLocalBlockSize(localBlockCol);

  const LocalOrdinal blkSize = numPtRows*numPtCols;
  Host_View_LO indx = node->template viewBuffer<LocalOrdinal>(1, pbuf_indx_+bpr[localBlockRow]+(it-bindx_beg));
  const LocalOrdinal offset = indx[0];
  blockEntry = node->template viewBuffer<Scalar>(blkSize, pbuf_values1D_ + offset);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
Teuchos::RCP<Node>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNode() const
{
  return getBlockRowMap()->getPointMap()->getNode();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::putScalar(Scalar s)
{
  Teuchos::RCP<Node> node = getNode();

  if (isFillComplete()) {
    fill_device_ArrayRCP(node, pbuf_values1D_, s);
  }
  else {
    typedef typename Teuchos::Array<RowGlobalCols>::size_type Tsize_t;
    Teuchos::Array<RowGlobalCols>& rows = *data_2D_;
    Tsize_t numBlockRows = rows.size();
    for(Tsize_t r=0; r<numBlockRows; ++r) {
      typename RowGlobalCols::iterator
        col_iter = rows[r].begin(),
        col_end  = rows[r].end();
      for(; col_iter != col_end; ++col_iter) {
        Teuchos::ArrayRCP<Scalar>& vals = col_iter->second;
        std::fill(vals.begin(), vals.end(), s);
      }
    }
  }
}

//-------------------------------------------------------------------
template<class ArrayType1,class ArrayType2,class Ordinal>
void set_array_values(ArrayType1& dest, ArrayType2& src, Ordinal rows, Ordinal cols, Ordinal stride, Tpetra::CombineMode mode)
{
  size_t offset = 0;
  size_t src_offset = 0;

  if (mode == ADD) {
    for(Ordinal col=0; col<cols; ++col) {
      for(Ordinal row=0; row<rows; ++row) {
        dest[offset++] += src[src_offset + row];
      }
      src_offset += stride;
    }
  }
  else {
    for(Ordinal col=0; col<cols; ++col) {
      for(Ordinal row=0; row<rows; ++row) {
        dest[offset++] = src[src_offset + row];
      }
      src_offset += stride;
    }
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  LocalOrdinal blkRowSize = blockEntry.numRows();
  LocalOrdinal blkColSize = blockEntry.numCols();
  getGlobalBlockEntryViewNonConst(globalBlockRow,globalBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  //now copy the incoming block-entry into internal storage:
  const Scalar* inputvalues = blockEntry.values();
  set_array_values(internalBlockEntry, inputvalues, blkRowSize, blkColSize, blockEntry.stride(), REPLACE);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  LocalOrdinal blkRowSize = blockEntry.numRows();
  LocalOrdinal blkColSize = blockEntry.numCols();
  getLocalBlockEntryViewNonConst(localBlockRow,localBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  //now copy the incoming block-entry into internal storage:
  const Scalar* inputvalues = blockEntry.values();
  set_array_values(internalBlockEntry, inputvalues, blkRowSize, blkColSize, blockEntry.stride(), REPLACE);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  LocalOrdinal blkRowSize = blockEntry.numRows();
  LocalOrdinal blkColSize = blockEntry.numCols();
  getGlobalBlockEntryViewNonConst(globalBlockRow,globalBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  //now sum (add) the incoming block-entry into internal storage:
  const Scalar* inputvalues = blockEntry.values();
  set_array_values(internalBlockEntry, inputvalues, blkRowSize, blkColSize, blockEntry.stride(), ADD);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, const Teuchos::SerialDenseMatrix<LocalOrdinal,Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  LocalOrdinal blkRowSize = blockEntry.numRows();
  LocalOrdinal blkColSize = blockEntry.numCols();
  getLocalBlockEntryViewNonConst(localBlockRow,localBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  //now sum (add) the incoming block-entry into internal storage:
  const Scalar* inputvalues = blockEntry.values();
  set_array_values(internalBlockEntry, inputvalues, blkRowSize, blkColSize, blockEntry.stride(), ADD);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  getGlobalBlockEntryViewNonConst(globalBlockRow,globalBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  LocalOrdinal blk_size = blockEntry.size();
  TEST_FOR_EXCEPTION(blkColSize*LDA > blk_size, std::runtime_error,
    "Tpetra::VbrMatrix::setGlobalBlockEntry ERROR, inconsistent block-entry sizes.");

  //copy the incoming block-entry into internal storage:
  set_array_values(internalBlockEntry, blockEntry, blkRowSize, blkColSize, LDA, REPLACE);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  getLocalBlockEntryViewNonConst(localBlockRow,localBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  LocalOrdinal blk_size = blockEntry.size();
  TEST_FOR_EXCEPTION(blkColSize*LDA > blk_size, std::runtime_error,
    "Tpetra::VbrMatrix::setLocalBlockEntry ERROR, inconsistent block-entry sizes.");

  //copy the incoming block-entry into internal storage:
  set_array_values(internalBlockEntry, blockEntry, blkRowSize, blkColSize, LDA, REPLACE);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  getGlobalBlockEntryViewNonConst(globalBlockRow,globalBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  LocalOrdinal blk_size = blockEntry.size();
  TEST_FOR_EXCEPTION(blkColSize*LDA > blk_size, std::runtime_error,
    "Tpetra::VbrMatrix::setGlobalBlockEntry ERROR, inconsistent block-entry sizes.");

  //copy the incoming block-entry into internal storage:
  set_array_values(internalBlockEntry, blockEntry, blkRowSize, blkColSize, LDA, ADD);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sumIntoLocalBlockEntry(LocalOrdinal localBlockRow, LocalOrdinal localBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  getLocalBlockEntryViewNonConst(localBlockRow,localBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  LocalOrdinal blk_size = blockEntry.size();
  TEST_FOR_EXCEPTION(blkColSize*LDA > blk_size, std::runtime_error,
    "Tpetra::VbrMatrix::setLocalBlockEntry ERROR, inconsistent block-entry sizes.");

  //copy the incoming block-entry into internal storage:
  set_array_values(internalBlockEntry, blockEntry, blkRowSize, blkColSize, LDA, ADD);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::optimizeStorage()
{
  typedef Teuchos::ArrayRCP<const size_t> Host_View;
  typedef Teuchos::ArrayRCP<const LocalOrdinal> Host_View_LO;
  typedef typename Teuchos::ArrayRCP<size_t>::size_type Tsize_t;

  if (is_storage_optimized_ == true) return;

  TEST_FOR_EXCEPTION(constBlkGraph_->isFillComplete() != true, std::runtime_error,
    "Tpetra::VbrMatrix::optimizeStorage ERROR, isFillComplete() is false, required to be true before optimizeStorage is called.");

  size_t num_block_nonzeros = constBlkGraph_->getNodeNumBlockEntries();

  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& rowmap = constBlkGraph_->getBlockRowMap();
  const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& colmap = constBlkGraph_->getBlockColMap();

  const Teuchos::RCP<Node>& node = getBlockRowMap()->getPointMap()->getNode();

  //need to count the number of point-entries:
  size_t num_point_entries = 0;
  Host_View bptr = constBlkGraph_->getNodeRowOffsets();

  if (bptr.size() == 0) {
    is_storage_optimized_ = true;
    return;
  }

  Host_View_LO bindx = constBlkGraph_->getNodePackedIndices();

  for(Tsize_t r=0; r<bptr.size()-1; ++r) {
    size_t rbeg = bptr[r];
    size_t rend = bptr[r+1];

    LocalOrdinal rsize = rowmap->getLocalBlockSize(r);

    for(size_t j=rbeg; j<rend; ++j) {
      LocalOrdinal csize = colmap->getLocalBlockSize(bindx[j]);
      num_point_entries += rsize*csize;
    }
  }

  pbuf_indx_ = node->template allocBuffer<LocalOrdinal>(num_block_nonzeros+1);
  pbuf_values1D_ = node->template allocBuffer<Scalar>(num_point_entries);

  Teuchos::ArrayRCP<LocalOrdinal> view_indx = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, num_block_nonzeros+1, pbuf_indx_);
  Teuchos::ArrayRCP<Scalar> view_values1D = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly, num_point_entries, pbuf_values1D_);

  bool have_2D_data = data_2D_->size() > 0;
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  size_t ioffset = 0;
  size_t offset = 0;
  for(Tsize_t r=0; r<bptr.size()-1; ++r) {
    LocalOrdinal rsize = rowmap->getLocalBlockSize(r);

    RowGlobalCols* blk_row = have_2D_data ? &((*data_2D_)[r]) : NULL;

    for(size_t c=bptr[r]; c<bptr[r+1]; ++c) {
      view_indx[ioffset++] = offset;

      LocalOrdinal csize = colmap->getLocalBlockSize(bindx[c]);
      Tsize_t blkSize = rsize*csize;

      if (blk_row != NULL) {
        GlobalOrdinal global_col = colmap->getGlobalBlockID(bindx[c]);
        typename RowGlobalCols::iterator iter = blk_row->find(global_col);
        TEST_FOR_EXCEPTION(iter == blk_row->end(), std::runtime_error,
          "Tpetra::VbrMatrix::optimizeStorage ERROR, global_col not found in row.");
  
        Teuchos::ArrayRCP<Scalar> vals = iter->second;
        for(Tsize_t n=0; n<blkSize; ++n) {
          view_values1D[offset+n] = vals[n];
        }
      }
      else {
        for(Tsize_t n=0; n<blkSize; ++n) view_values1D[offset+n] = zero;
      }
      offset += blkSize;
    }
  }
  view_indx[ioffset] = offset;

  //Final step: release memory for the "2D" storage:
  data_2D_ = Teuchos::null;

  is_storage_optimized_ = true;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillLocalMatrix()
{
  //We insist that optimzeStorage has already been called.
  //We don't care whether this function (fillLocalMatrix()) is being
  //called for the first time or not.
  TEST_FOR_EXCEPTION(is_storage_optimized_ != true, std::runtime_error,
    "Tpetra::VbrMatrix::fillLocalMatrix ERROR, optimizeStorage is required to have already been called.");

  Teuchos::ArrayRCP<const size_t      > nodeRowOffsets = constBlkGraph_->getNodeRowOffsets();
  Teuchos::ArrayRCP<const LocalOrdinal> nodePackedInds = constBlkGraph_->getNodePackedIndices();
  if (Node::isHostNode == false) {
    RCP<Node> node = getRangeMap()->getNode();
    //
    Teuchos::ArrayRCP<size_t> dev_nodeRowOffsets = node->template allocBuffer<size_t>(nodeRowOffsets.size());
    node->copyToBuffer(nodeRowOffsets.size(), nodeRowOffsets(), dev_nodeRowOffsets);
    nodeRowOffsets = dev_nodeRowOffsets;
    //
    Teuchos::ArrayRCP<LocalOrdinal> dev_nodePackedInds = node->template allocBuffer<LocalOrdinal>(nodePackedInds.size());
    node->copyToBuffer(nodePackedInds.size(), nodePackedInds(), dev_nodePackedInds);
    nodePackedInds = dev_nodePackedInds;
  }
  lclMatrix_.setPackedValues(pbuf_values1D_,
                             getBlockRowMap()->getNodeFirstPointInBlocks_Device(),
                             getBlockColMap()->getNodeFirstPointInBlocks_Device(),
                             nodeRowOffsets,
                             nodePackedInds,
                             pbuf_indx_);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillLocalMatVec()
{
  //We insist that optimzeStorage has already been called.
  //We don't care whether this function (fillLocalMatVec()) is being
  //called for the first time or not.
  TEST_FOR_EXCEPTION(is_storage_optimized_ != true, std::runtime_error,
    "Tpetra::VbrMatrix::fillLocalMatrix ERROR, optimizeStorage is required to have already been called.");

  lclMatOps_.initializeValues(lclMatrix_);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::createImporterExporter()
{
  typedef typename Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > PtMap;

  const PtMap& ptDomMap = (constBlkGraph_->getBlockDomainMap())->getPointMap();
  const PtMap& ptRngMap = (constBlkGraph_->getBlockRangeMap())->getPointMap();
  const PtMap& ptRowMap = (constBlkGraph_->getBlockRowMap())->getPointMap();
  const PtMap& ptColMap = (constBlkGraph_->getBlockColMap())->getPointMap();

  if (!ptDomMap->isSameAs(*ptColMap)) {
    importer_ = Teuchos::rcp(new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(ptDomMap, ptColMap));
  }
  if (!ptRngMap->isSameAs(*ptRowMap)) {
    exporter_ = Teuchos::rcp(new Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node>(ptRowMap, ptRngMap));
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap)
{
  if (isFillComplete()) return;

  globalAssemble();

  blkGraph_->fillComplete(blockDomainMap, blockRangeMap);

  createImporterExporter();

  is_fill_completed_ = true;

  optimizeStorage();

  fillLocalMatrix();
  fillLocalMatVec();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete()
{
  //In this case, our block-row-map will also be our domain-map and
  //range-map.

  fillComplete(getBlockRowMap(), getBlockRowMap());
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::globalAssemble()
{
  if (getPointRowMap()->getComm()->getSize() == 1) {
    return;
  }

  //nonlocal_data_ contains data that belongs on other processors.
  //We'll refer to this as overlapping data, and create an 'overlapMap'
  //that describes the layout of this overlapping data.
  Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >
     overlapMap = VbrUtils::createOverlapMap(nonlocal_data_, *getBlockRowMap());

  if (overlapMap->getGlobalNumBlocks() == 0) {
    return;
  }

  //VbrDataDist is a wrapper that makes our overlapping data behave
  //like a DistObject so that the overlapping data can be exported
  //to the owning processors in a standard way.
  VbrUtils::VbrDataDist<LocalOrdinal,GlobalOrdinal,Scalar,Node>
      vdd(nonlocal_data_, *overlapMap, *this->getBlockRowMap());

  //Create an exporter where the source map is our overlapMap and the
  //target map is the rowmap of our VbrMatrix.
  Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>
      importer( vdd.getMap(), convertBlockMapToPointMap(*getBlockRowMap()));

  //Communicate the overlapping data to the owning processors and add it
  //into this VbrMatrix.
  this->doImport(vdd, importer, Tpetra::ADD);

  //Zero out the overlapping data so it can be re-populated and re-assembled
  //in future calls to globalAssemble.
  VbrUtils::zeroEntries(nonlocal_data_);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::string
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isFillComplete()) {
    oss << "{status = fill complete, global num block rows = " << getBlockRowMap()->getGlobalNumBlocks() << "}";
  }
  else {
    oss << "{status = fill not complete, global num block rows = " << getBlockRowMap()->getGlobalNumBlocks() << "}";
  }
  return oss.str();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == Teuchos::VERB_DEFAULT) vl = Teuchos::VERB_LOW;

  if (vl == Teuchos::VERB_NONE) return;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = getBlockRowMap()->getPointMap()->getComm();
  const int myProc = comm->getRank();

  if (myProc == 0) out << "VbrMatrix::describe is under construction..." << std::endl;
}

}//namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_VBRMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class VbrMatrix< SCALAR , LO , GO , NODE >;

#endif //TPETRA_VBRMATRIX_DEF_HPP

