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

#ifndef TPETRA_CRSMATRIX_DEF_HPP
#define TPETRA_CRSMATRIX_DEF_HPP

// FINISH: need to check that fill is active before performing a number of the methods here; adding to the tests currently

// TODO: row-wise insertion of entries in globalAssemble() may be more efficient
// TODO: consider maintaining sorted entries at all times and leaning heavily on STL set_intersect/set_union methods for all insert/replace/suminto

#include <Kokkos_NodeHelpers.hpp>
#include <Kokkos_NodeTrace.hpp>

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

  template <class Ordinal, class Scalar>
  CrsIJV<Ordinal,Scalar>::CrsIJV() {}

  template <class Ordinal, class Scalar>
  CrsIJV<Ordinal,Scalar>::CrsIJV(Ordinal row, Ordinal col, const Scalar &val) {
    i = row;
    j = col;
    v = val;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::CrsMatrix(
                                          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          size_t maxNumEntriesPerRow, 
                                          ProfileType pftype)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , lclMatOps_(rowMap->getNode())
  {
    try {
      myGraph_ = rcp( new Graph(rowMap,maxNumEntriesPerRow,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    staticGraph_ = myGraph_;
    lclMatrix_.setOwnedGraph(myGraph_->getLocalGraphNonConst());
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );
    resumeFill();
    //
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::CrsMatrix(
                                          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, 
                                          ProfileType pftype)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , lclMatOps_(rowMap->getNode())
  {
    try {
      myGraph_ = rcp( new Graph(rowMap,NumEntriesPerRowToAlloc,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    staticGraph_ = myGraph_;
    lclMatrix_.setOwnedGraph(myGraph_->getLocalGraphNonConst());
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );
    resumeFill();
    //
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::CrsMatrix(
                                          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, 
                                          size_t maxNumEntriesPerRow, 
                                          ProfileType pftype)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , lclMatOps_(rowMap->getNode())
  {
    try {
      myGraph_ = rcp( new Graph(rowMap,colMap,maxNumEntriesPerRow,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    staticGraph_ = myGraph_;
    lclMatrix_.setOwnedGraph(myGraph_->getLocalGraphNonConst());
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );
    resumeFill();
    //
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::CrsMatrix(
                                          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, 
                                          const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, 
                                          ProfileType pftype)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , lclMatOps_(rowMap->getNode())
  {
    try {
      myGraph_ = rcp( new Graph(rowMap,colMap,NumEntriesPerRowToAlloc,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    staticGraph_ = myGraph_;
    lclMatrix_.setOwnedGraph(myGraph_->getLocalGraphNonConst());
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );
    resumeFill();
    //
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::CrsMatrix(const RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &graph)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(graph->getRowMap())
  , staticGraph_(graph)
  , lclMatOps_(graph->getNode())
  {
    const std::string tfecfFuncName("CrsMatrix(graph)");
    TEST_FOR_EXCEPTION_CLASS_FUNC(staticGraph_ == null, std::runtime_error, ": specified pointer is null.");
    // we prohibit the case where the graph is not yet filled
    TEST_FOR_EXCEPTION_CLASS_FUNC( staticGraph_->isFillComplete() == false, std::runtime_error, 
        ": specified graph is not fill-complete. You must fillComplete() the graph before using it to construct a CrsMatrix.");
    lclMatrix_.setStaticGraph(staticGraph_->getLocalGraph());
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );
    // the graph has entries, and the matrix should have entries as well, set to zero. no need or point in lazy allocating in this case.
    // first argument doesn't actually matter, because the graph is allocated.
    allocateValues( LocalIndices, GraphAlreadyAllocated );
    resumeFill();
    //
    checkInternalState();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::~CrsMatrix() {
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Comm<int> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getComm() const {
    return getCrsGraph()->getComm();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Node>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNode() const {
    return lclMatOps_.getNode();
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


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                    Internal utility methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::allocateValues(ELocalGlobal lg, GraphAllocationStatus gas) {
    // allocate values and, optionally, ask graph to allocate indices
#ifdef HAVE_TPETRA_DEBUG
    // if the graph is already allocated, then gas should be GraphAlreadyAllocated
    // otherwise, gas should be GraphNotYetAllocated
    // this method is for internal use only. debug checks occur outside. only do them here if debugging is enabled.
    std::string err("::allocateValues(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION((gas == GraphAlreadyAllocated) != staticGraph_->indicesAreAllocated(), std::logic_error, typeName(*this) << err);
    // if the graph is unallocated, then it better be a matrix-owned graph
    TEST_FOR_EXCEPTION(staticGraph_->indicesAreAllocated() == false && myGraph_ == null, std::logic_error, typeName(*this) << err);
#endif
    if (gas == GraphNotYetAllocated) {
      myGraph_->allocateIndices(lg);
    }
    // ask graph to allocate our values, with the same structure
    // this will allocate values2D_ one way or the other
    staticGraph_->template allocateValues<Scalar>(values1D_, values2D_);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::pushToLocalMatrix() {
    std::string err("::pushToLocalMatrix(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(lclMatrix_.isFinalized() == true, std::logic_error, typeName(*this) << err);
    // fill local graph
    if (getProfileType() == StaticProfile) {
      lclMatrix_.set1DValues(values1D_);
      values1D_ = null;
    }
    else {
      lclMatrix_.set2DValues(values2D_);
      values2D_ = null;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::pullFromLocalMatrix() {
    std::string err("::pullFromLocalMatrix(): Internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(lclMatrix_.isFinalized() == false, std::logic_error, typeName(*this) << err);
    // get new data from local matrix
    if (lclMatrix_.is1DStructure()) {
      lclMatrix_.get1DValues(values1D_);
    }
    else {
      lclMatrix_.get2DValues(values2D_);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillLocalMatrix(OptimizeOption os) {
    // if not static graph: fill local graph
    if (isStaticGraph() == false) {
      myGraph_->pushToLocalGraph();
    }
    pushToLocalMatrix();
    // finalize local matrix(with os) (this will finalize local graph if not static)
    const bool optStorage = (os == DoOptimizeStorage);
    lclMatrix_.finalize( optStorage );
    // get the data back from the local objects
    if (isStaticGraph() == false) {
      myGraph_->pullFromLocalGraph();
    }
    pullFromLocalMatrix();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillLocalSparseOps() 
  {
    lclMatOps_.initializeStructure(staticGraph_->getLocalGraph());
    lclMatOps_.initializeValues(lclMatrix_);
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
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::insertLocalValues(
                                                         LocalOrdinal localRow, 
                         const ArrayView<const LocalOrdinal> &indices,
                         const ArrayView<const Scalar>       &values) 
  {
    const std::string tfecfFuncName("insertLocalValues()");
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false,                            std::runtime_error, " requires that fill is active.");
    TEST_FOR_EXCEPTION_CLASS_FUNC( isStaticGraph() == true,                            std::runtime_error, " cannot insert indices with static graph; use replaceLocalValues() instead.");
    TEST_FOR_EXCEPTION_CLASS_FUNC( myGraph_->isGloballyIndexed() == true,              std::runtime_error, ": graph indices are global; use insertGlobalValues().");
    TEST_FOR_EXCEPTION_CLASS_FUNC( hasColMap() == false,                               std::runtime_error, " cannot insert local indices without a column map.");
    TEST_FOR_EXCEPTION_CLASS_FUNC( values.size() != indices.size(),                    std::runtime_error, ": values.size() must equal indices.size().");
    TEST_FOR_EXCEPTION_CLASS_FUNC( getRowMap()->isNodeLocalElement(localRow) == false, std::runtime_error, ": row does not belong to this node.");
    if (myGraph_->indicesAreAllocated() == false) {
      allocateValues(LocalIndices, GraphNotYetAllocated);
    }
    // use column map to filter the entries:
    Array<LocalOrdinal> f_inds(indices);
    Array<Scalar>       f_vals(values);
    typename Graph::SLocalGlobalNCViews inds_ncview;
    inds_ncview.linds = f_inds();
    const size_t numFilteredEntries = myGraph_->template filterIndicesAndValues<LocalIndices,Scalar>(inds_ncview, f_vals());
    if (numFilteredEntries > 0) {
      RowInfo rowInfo = myGraph_->getRowInfo(localRow);
      const size_t curNumEntries = rowInfo.numEntries;
      const size_t newNumEntries = curNumEntries + numFilteredEntries;
      if (newNumEntries > rowInfo.allocSize) {
        TEST_FOR_EXCEPTION_CLASS_FUNC(getProfileType() == StaticProfile, std::runtime_error, ": new indices exceed statically allocated graph structure.");
        TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
            "::insertLocalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
        // update allocation only as much as necessary
        rowInfo = myGraph_->template updateAllocAndValues<LocalIndices,Scalar>(rowInfo, newNumEntries, values2D_[localRow]);
      }
      typename Graph::SLocalGlobalViews inds_view;
      inds_view.linds = f_inds(0,numFilteredEntries);
      myGraph_->template insertIndicesAndValues<LocalIndices,LocalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), f_vals.begin());
#ifdef HAVE_TPETRA_DEBUG
      {
        const size_t chkNewNumEntries = myGraph_->getNumEntriesInLocalRow(localRow);
        TEST_FOR_EXCEPTION_CLASS_FUNC(chkNewNumEntries != newNumEntries, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
      }
#endif
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() == false, std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::insertGlobalValues(
                                                        GlobalOrdinal globalRow, 
                       const ArrayView<const GlobalOrdinal> &indices,
                       const ArrayView<const Scalar>        &values) 
  {
    const std::string tfecfFuncName("insertGlobalValues()");
    TEST_FOR_EXCEPTION_CLASS_FUNC(isStaticGraph() == true,         std::runtime_error, ": matrix was constructed with static graph. Cannot insert new entries.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(), std::runtime_error, ": values.size() must equal indices.size().");
    if (myGraph_->indicesAreAllocated() == false) {
      allocateValues(GlobalIndices, GraphNotYetAllocated);
    }
    const LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    typename Graph::SLocalGlobalViews         inds_view;
    ArrayView<const Scalar> vals_view;
    if (lrow != LOT::invalid()) {
      // if we have a column map, use it to filter the entries.
      Array<GlobalOrdinal> filtered_indices;
      Array<Scalar>        filtered_values;
      if (hasColMap()) {
        typename Graph::SLocalGlobalNCViews inds_ncview;
        ArrayView<Scalar> vals_ncview;
        // filter indices and values through the column map
        filtered_indices.assign(indices.begin(), indices.end());
        filtered_values.assign(values.begin(), values.end());
        inds_ncview.ginds = filtered_indices();
        const size_t numFilteredEntries = myGraph_->template filterIndicesAndValues<GlobalIndices,Scalar>(inds_ncview,filtered_values());
        inds_view.ginds = filtered_indices(0,numFilteredEntries);
        vals_view       = filtered_values(0,numFilteredEntries);
      }
      else {
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
          TEST_FOR_EXCEPTION_CLASS_FUNC(getProfileType() == StaticProfile, std::runtime_error, ": new indices exceed statically allocated graph structure.");
          TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
              "::insertGlobalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
          // update allocation only as much as necessary
          rowInfo = myGraph_->template updateAllocAndValues<GlobalIndices,Scalar>(rowInfo, newNumEntries, values2D_[lrow]);
        }
        if (isGloballyIndexed()) {
          myGraph_->template insertIndicesAndValues<GlobalIndices,GlobalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), vals_view.begin());
        }
        else {
          myGraph_->template insertIndicesAndValues<GlobalIndices,LocalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), vals_view.begin());
        }
#ifdef HAVE_TPETRA_DEBUG
        {
          const size_t chkNewNumEntries = myGraph_->getNumEntriesInLocalRow(lrow);
          TEST_FOR_EXCEPTION_CLASS_FUNC(chkNewNumEntries != newNumEntries, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
        }
#endif
      }
    }
    else {
      typename ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
      typename ArrayView<const Scalar       >::iterator val =  values.begin();
      nonlocals_[globalRow].reserve( nonlocals_[globalRow].size() + indices.size() );
      for (; val != values.end(); ++val, ++ind) {
        nonlocals_[globalRow].push_back(std::make_pair(*ind, *val));
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::replaceLocalValues(      
                                        LocalOrdinal localRow,
                                        const ArrayView<const LocalOrdinal> &indices,
                                        const ArrayView<const Scalar>        &values) {
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const std::string tfecfFuncName("replaceLocalValues()");
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false,        std::runtime_error, " requires that fill is active.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(), std::runtime_error, ": values.size() must equal indices.size().");
    bool isLocalRow = getRowMap()->isNodeLocalElement(localRow);
    TEST_FOR_EXCEPTION_CLASS_FUNC(hasColMap() == false,            std::runtime_error, ": cannot replace local indices without a column map.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(isLocalRow == false,             std::runtime_error, ": specified local row does not belong to this processor.");
    // 
    RowInfo rowInfo = staticGraph_->getRowInfo(localRow);
    if (indices.size() > 0) {
      if (isLocallyIndexed() == true) {
        typename Graph::SLocalGlobalViews inds_view;
        inds_view.linds = indices;
        staticGraph_->template transformValues<LocalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), secondArg<Scalar,Scalar>());
      }
      else if (isGloballyIndexed() == true) {
        // must convert to global indices
        const Map<LocalOrdinal,GlobalOrdinal,Node> &colMap = *getColMap();
        Array<GlobalOrdinal> gindices(indices.size());
        typename ArrayView<const LocalOrdinal>::iterator lindit = indices.begin();
        typename Array<GlobalOrdinal>::iterator          gindit = gindices.begin();
        while (lindit != indices.end()) {
          // no need to filter: if it doesn't exist, it will be mapped to invalid(), which will not be found in the graph. 
          *gindit++ = colMap.getGlobalElement(*lindit++);
        }
        typename Graph::SLocalGlobalViews inds_view;
        inds_view.ginds = gindices();
        staticGraph_->template transformValues<GlobalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), secondArg<Scalar,Scalar>());
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::replaceGlobalValues(      
                                        GlobalOrdinal globalRow, 
                                        const ArrayView<const GlobalOrdinal> &indices,
                                        const ArrayView<const Scalar>        &values) {
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const std::string tfecfFuncName("replaceGlobalValues()");
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false,         std::runtime_error, " requires that fill is active.");
    TEST_FOR_EXCEPTION_CLASS_FUNC( values.size() != indices.size(), std::runtime_error, " values.size() must equal indices.size().");
    const LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION_CLASS_FUNC( lrow == LOT::invalid(), std::runtime_error,          ": specified global row does not belong to this processor.");
    // 
    RowInfo rowInfo = staticGraph_->getRowInfo(lrow);
    if (indices.size() > 0) {
      if (isLocallyIndexed() == true) {
        // must convert global indices to local indices
        const Map<LocalOrdinal,GlobalOrdinal,Node> &colMap = *getColMap();
        Array<LocalOrdinal> lindices(indices.size());
        typename ArrayView<const GlobalOrdinal>::iterator gindit = indices.begin();
        typename Array<LocalOrdinal>::iterator            lindit = lindices.begin();
        while (gindit != indices.end()) {
          // no need to filter: if it doesn't exist, it will be mapped to invalid(), which will not be found in the graph. 
          *lindit++ = colMap.getLocalElement(*gindit++);
        }
        typename Graph::SLocalGlobalViews inds_view;
        inds_view.linds = lindices();
        staticGraph_->template transformValues<LocalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), secondArg<Scalar,Scalar>());
      }
      else if (isGloballyIndexed() == true) {
        typename Graph::SLocalGlobalViews inds_view;
        inds_view.ginds = indices;
        staticGraph_->template transformValues<GlobalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), secondArg<Scalar,Scalar>());
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sumIntoGlobalValues(GlobalOrdinal globalRow, 
                         const ArrayView<const GlobalOrdinal> &indices,
                         const ArrayView<const Scalar>        &values) 
  {
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const std::string tfecfFuncName("sumIntoGlobalValues()");
    TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(), std::runtime_error, ": values.size() must equal indices.size().");
    const LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION_CLASS_FUNC(lrow == LOT::invalid(),          std::runtime_error, ": specified global row does not belong to this processor.");
    // 
    RowInfo rowInfo = staticGraph_->getRowInfo(lrow);
    if (indices.size() > 0) {
      if (isLocallyIndexed() == true) {
        // must convert global indices to local indices
        const Map<LocalOrdinal,GlobalOrdinal,Node> &colMap = *getColMap();
        Array<LocalOrdinal> lindices(indices.size());
        typename ArrayView<const GlobalOrdinal>::iterator gindit = indices.begin();
        typename Array<LocalOrdinal>::iterator            lindit = lindices.begin();
        while (gindit != indices.end()) {
          // no need to filter: if it doesn't exist, it will be mapped to invalid(), which will not be found in the graph. 
          *lindit++ = colMap.getLocalElement(*gindit++);
        }
        typename Graph::SLocalGlobalViews inds_view;
        inds_view.linds = lindices();
        staticGraph_->template transformValues<LocalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), std::plus<Scalar>());
      }
      else if (isGloballyIndexed() == true) {
        typename Graph::SLocalGlobalViews inds_view;
        inds_view.ginds = indices;
        staticGraph_->template transformValues<GlobalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), std::plus<Scalar>());
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sumIntoLocalValues(LocalOrdinal localRow, 
                         const ArrayView<const LocalOrdinal>  &indices,
                         const ArrayView<const Scalar>        &values) 
  {
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const std::string tfecfFuncName("sumIntoLocalValues()");
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false,                           std::runtime_error, " requires that fill is active.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(values.size() != indices.size(),                    std::runtime_error, ": values.size() must equal indices.size().");
    TEST_FOR_EXCEPTION_CLASS_FUNC(getRowMap()->isNodeLocalElement(localRow) == false, std::runtime_error, ": specified local row does not belong to this processor.");
    // 
    RowInfo rowInfo = staticGraph_->getRowInfo(localRow);
    if (indices.size() > 0) {
      if (isGloballyIndexed() == true) {
        // must convert local indices to global indices
        const Map<LocalOrdinal,GlobalOrdinal,Node> &colMap = *getColMap();
        Array<GlobalOrdinal> gindices(indices.size());
        typename ArrayView<const LocalOrdinal>::iterator lindit = indices.begin();
        typename Array<GlobalOrdinal>::iterator          gindit = gindices.begin();
        while (lindit != indices.end()) {
          // no need to filter: if it doesn't exist, it will be mapped to invalid(), which will not be found in the graph. 
          *gindit++ = colMap.getGlobalElement(*lindit++);
        }
        typename Graph::SLocalGlobalViews inds_view;
        inds_view.ginds = gindices();
        staticGraph_->template transformValues<GlobalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), std::plus<Scalar>());
      }
      else if (isLocallyIndexed() == true) {
        typename Graph::SLocalGlobalViews inds_view;
        inds_view.linds = indices;
        staticGraph_->template transformValues<LocalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), std::plus<Scalar>());
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ArrayView<const Scalar> CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getView(RowInfo rowinfo) const
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
    const std::string tfecfFuncName("getLocalRowCopy()");
    TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed()==true && hasColMap()==false, std::runtime_error,
        ": local indices cannot be produced.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(getRowMap()->isNodeLocalElement(localRow) == false, std::runtime_error,
        ": specified row (==" << localRow << ") is not valid on this node.");
    const RowInfo rowinfo = staticGraph_->getRowInfo(localRow);
    numEntries = rowinfo.numEntries;
    TEST_FOR_EXCEPTION_CLASS_FUNC(static_cast<size_t>(indices.size()) < numEntries || static_cast<size_t>(values.size()) < numEntries, 
        std::runtime_error, ": size of indices,values must be sufficient to store the specified row.");
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
      TEST_FOR_EXCEPTION_CLASS_FUNC( staticGraph_->indicesAreAllocated() == true, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
#endif
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
    TEST_FOR_EXCEPTION_CLASS_FUNC(lrow == LOT::invalid(), std::runtime_error, ": globalRow does not belong to this node.");
    const RowInfo rowinfo = staticGraph_->getRowInfo(lrow);
    numEntries = rowinfo.numEntries;
    TEST_FOR_EXCEPTION_CLASS_FUNC(static_cast<size_t>(indices.size()) < numEntries || static_cast<size_t>(values.size()) < numEntries, 
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
      TEST_FOR_EXCEPTION_CLASS_FUNC( staticGraph_->indicesAreAllocated() == true, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
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
    TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed() == true, std::runtime_error, ": local indices cannot be provided.");
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
    TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)indices.size() != getNumEntriesInLocalRow(localRow) || indices.size() != values.size(), 
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
    TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() == true, std::runtime_error, ": global indices cannot be provided.");
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
    TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)indices.size() != getNumEntriesInGlobalRow(globalRow) || indices.size() != values.size(), 
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
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false, std::runtime_error, " requires that fill is active.");
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
    const std::string tfecfFuncName("scale()");
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false, std::runtime_error, " requires that fill is active.");
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
    TEST_FOR_EXCEPTION_CLASS_FUNC(isFillComplete() == false, std::runtime_error, " until fillComplete() has been called.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(dvec.getMap()->isSameAs(*getRowMap()) == false, std::runtime_error, ": dvec must have the same map as the CrsMatrix.");
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
    TEST_FOR_EXCEPTION_CLASS_FUNC(numDiagFound != getNodeNumDiags(), std::logic_error, ": logic error. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::globalAssemble() 
  {
    using Teuchos::SerialDenseMatrix;
    using std::pair;
    using std::make_pair;
    typedef typename std::map<GlobalOrdinal,Array<pair<GlobalOrdinal,Scalar> > >::const_iterator NLITER;
    typedef typename Array<pair<GlobalOrdinal,Scalar> >::const_iterator NLRITER;
    const int numImages = getComm()->getSize();
    const int myImageID = getComm()->getRank();
    const std::string tfecfFuncName("globalAssemble()");
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *getRowMap()->getComm() );
#endif
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false, std::runtime_error, " requires that fill is active.");
    // Determine if any nodes have global entries to share
    size_t MyNonlocals = nonlocals_.size(), 
           MaxGlobalNonlocals;
    Teuchos::reduceAll<int,size_t>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,
      outArg(MaxGlobalNonlocals));
    if (MaxGlobalNonlocals == 0) return;  // no entries to share

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
        TEST_FOR_EXCEPTION_CLASS_FUNC(gblerror, std::runtime_error, ": non-local entries correspond to invalid rows.");
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
        TEST_FOR_EXCEPTION_CLASS_FUNC(sendIDs[numSends] != id, std::logic_error, ": internal logic error. Contact Tpetra team.");
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
    TEST_FOR_EXCEPTION_CLASS_FUNC(Teuchos::as<typename Array<int>::size_type>(numSends) != sendIDs.size(), 
        std::logic_error, ": internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    nonlocals_.clear();

    ////////////////////////////////////////////////////////////////////////////////////// 
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    ////////////////////////////////////////////////////////////////////////////////////// 
    // perform non-blocking sends: send sizes to our recipients
    Array<RCP<Teuchos::CommRequest> > sendRequests;
    for (size_t s=0; s < numSends ; ++s) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back( Teuchos::isend<int,size_t>(*getComm(),rcpFromRef(sendSizes[s]),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<RCP<Teuchos::CommRequest> > recvRequests;
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
    for (typename Array<CrsIJV<GlobalOrdinal,Scalar> >::const_iterator ijv = IJVRecvBuffer.begin(); ijv != IJVRecvBuffer.end(); ++ijv)
    {
      try {
        insertGlobalValues(ijv->i, tuple(ijv->j), tuple(ijv->v));
      }
      catch (std::runtime_error &e) {
        std::ostringstream omsg;
        omsg << e.what() << std::endl
          << "caught in globalAssemble() in " << __FILE__ << ":" << __LINE__ << std::endl ;
        throw std::runtime_error(omsg.str());
      }
    }

    // WHEW! THAT WAS TIRING!
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::resumeFill() {
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *getRowMap()->getComm() );
#endif
    if (isStaticGraph() == false) {
      myGraph_->resumeFill();
    }
    clearGlobalConstants();
    lclMatrix_.clear();
    lclMatOps_.clear();
    fillComplete_ = false;
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION( isFillActive() == false || isFillComplete() == true, std::logic_error,
        typeName(*this) << "::resumeFill(): Violated stated post-conditions. Please contact Tpetra team.");
#endif
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
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(OptimizeOption os) {
    fillComplete(getRowMap(),getRowMap(),os);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(
                                            const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, 
                                            const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, 
                                            OptimizeOption os) 
  {
    const std::string tfecfFuncName("fillComplete()");
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false || isFillComplete() == true, std::runtime_error, ": Matrix fill state must be active.");
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *getRowMap()->getComm() );
#endif
    // allocate if unallocated
    if (getCrsGraph()->indicesAreAllocated() == false) {
      // allocate global, in case we do not have a column map
      allocateValues( GlobalIndices, GraphNotYetAllocated );
    }
    // global assemble
    if (getComm()->getSize() > 1) {
      globalAssemble();
    }
    else {
      TEST_FOR_EXCEPTION_CLASS_FUNC(nonlocals_.size() > 0, std::runtime_error, ": cannot have non-local entries on a serial run. Invalid entry was submitted to the CrsMatrix.");
    }
    //
    // if we're not allowed to change a static graph, then we can't call optimizeStorage() on it.
    // then we can't call fillComplete() on the matrix with DoOptimizeStorage
    // throw a warning. this is an unfortunate late evaluation; however, we couldn't know when we received the graph 
    // that the user would try to optimize the storage later on.
    if (os == DoOptimizeStorage && isStaticGraph() && staticGraph_->isStorageOptimized() == false) {
      TPETRA_ABUSE_WARNING(true,std::runtime_error,
          "::fillComplete(): requested optimized storage, but static graph does not have optimized storage. Ignoring request to optimize storage.");
      os = DoNotOptimizeStorage;
    }
    //
    if (isStaticGraph() == true) {
      TEST_FOR_EXCEPTION_CLASS_FUNC((staticGraph_->getDomainMap() != getDomainMap()) || (staticGraph_->getRangeMap() != getRangeMap()), std::runtime_error,
          ": domain map and range map do not match maps in existing graph, and the graph cannot be changed because it was specified during matrix construction.");
    }
    else {
      // set domain/range map: may clear the import/export objects
      myGraph_->setDomainRangeMaps(domainMap, rangeMap);
      // make column map
      if (myGraph_->hasColMap() == false) {
        myGraph_->makeColMap();
      }
      // make indices local
      if (myGraph_->isGloballyIndexed() == true) {
        myGraph_->makeIndicesLocal();
      }
      // sort entries
      sortEntries();
      // merge entries
      mergeRedundantEntries();
      // make import/export objects
      myGraph_->makeImportExport();
      // compute global constants
      myGraph_->computeGlobalConstants();
      myGraph_->fillComplete_ = true;
      myGraph_->checkInternalState();
    }
    computeGlobalConstants();
    // fill local objects; will fill and finalize local graph if appropriate
    fillLocalMatrix(os);
    fillLocalSparseOps();
    //
    fillComplete_ = true;
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == true || isFillComplete() == false, std::logic_error,
        ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    //
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sortEntries() 
  {
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
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
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
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
    TEST_FOR_EXCEPTION( isFillComplete() == false, std::runtime_error, 
        typeName(*this) << "::apply(): cannot call apply() until fillComplete() has been called.");
    sameScalarMultiplyOp_->apply(X,Y,mode,alpha,beta);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class DomainScalar, class RangeScalar>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::multiply(
                                        const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                                              MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                                              Teuchos::ETransp mode, RangeScalar alpha, RangeScalar beta) const 
  {
    using Teuchos::NO_TRANS;
    const std::string tfecfFuncName("multiply()");
    typedef ScalarTraits<RangeScalar> RST;
    const Kokkos::MultiVector<DomainScalar,Node> *lclX = &X.getLocalMV();
    Kokkos::MultiVector<RangeScalar,Node>        *lclY = &Y.getLocalMVNonConst();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION_CLASS_FUNC(mode == NO_TRANS && X.getMap() != getColMap() && *X.getMap() != *getColMap(), std::runtime_error, " X is not distributed according to the appropriate map.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(mode != NO_TRANS && X.getMap() != getRowMap() && *X.getMap() != *getRowMap(), std::runtime_error, " X is not distributed according to the appropriate map.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(mode == NO_TRANS && Y.getMap() != getRowMap() && *Y.getMap() != *getRowMap(), std::runtime_error, " Y is not distributed according to the appropriate map.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(mode != NO_TRANS && Y.getMap() != getColMap() && *Y.getMap() != *getColMap(), std::runtime_error, " Y is not distributed according to the appropriate map.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(!isFillComplete(),                                              std::runtime_error, " until fillComplete() has been called.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(X.getNumVectors() != Y.getNumVectors(),                         std::runtime_error, ": X and Y must have the same number of vectors.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(X.isConstantStride() == false || Y.isConstantStride() == false, std::runtime_error, ": X and Y must be constant stride.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(lclX==lclY,                                                     std::runtime_error, ": X and Y cannot share data.");
#endif
    //
    // Call the matvec
    if (beta == RST::zero()) {
      // Y = alpha*op(M)*X with overwrite semantics
      lclMatOps_.template multiply<DomainScalar,RangeScalar>(mode, alpha, *lclX, *lclY);
    }
    else {
      // Y = alpha*op(M) + beta*Y
      lclMatOps_.template multiply<DomainScalar,RangeScalar>(mode, alpha, *lclX, beta, *lclY);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class DomainScalar, class RangeScalar>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::solve(
                                    const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>  &Y, 
                                          MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                                          Teuchos::ETransp mode) const 
  {
    using Teuchos::NO_TRANS;
    const std::string tfecfFuncName("solve()");
    const Kokkos::MultiVector<RangeScalar,Node> *lclY = &Y.getLocalMV();
    Kokkos::MultiVector<DomainScalar,Node>      *lclX = &X.getLocalMVNonConst();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION_CLASS_FUNC(!isFillComplete(),                                              std::runtime_error, " until fillComplete() has been called.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(X.getNumVectors() != Y.getNumVectors(),                         std::runtime_error, ": X and Y must have the same number of vectors.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(X.isConstantStride() == false || Y.isConstantStride() == false, std::runtime_error, ": X and Y must be constant stride.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(isUpperTriangular() == false && isLowerTriangular() == false,   std::runtime_error, ": can only solve() triangular matrices.");
    TEST_FOR_EXCEPTION_CLASS_FUNC(ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS,      std::logic_error, " does not currently support transposed solve for complex scalar types.");
#endif
    //
    // Call the solve
    Teuchos::EDiag diag = ( getNodeNumDiags() < getNodeNumRows() ? Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG );
    if (mode == Teuchos::NO_TRANS) {
      if (isUpperTriangular()) {
        lclMatOps_.template solve<DomainScalar,RangeScalar>(Teuchos::NO_TRANS, Teuchos::UPPER_TRI, diag, *lclY, *lclX);
      }
      else {
        lclMatOps_.template solve<DomainScalar,RangeScalar>(Teuchos::NO_TRANS, Teuchos::LOWER_TRI, diag, *lclY, *lclX);
      }
    }
    else {
      if (isUpperTriangular()) {
        lclMatOps_.template solve<DomainScalar,RangeScalar>(Teuchos::CONJ_TRANS, Teuchos::UPPER_TRI, diag, *lclY, *lclX);
      }
      else {
        lclMatOps_.template solve<DomainScalar,RangeScalar>(Teuchos::CONJ_TRANS, Teuchos::LOWER_TRI, diag, *lclY, *lclX);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  template <class T>
  RCP<CrsMatrix<T,LocalOrdinal,GlobalOrdinal,Node> > 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::convert() const
  {
    const std::string tfecfFuncName("convert()");
    // FINISH: we have a problem here: converted matrices will be statically allocated, and therefore will not benefit from first touch 
    // allocation. must address this in the future.
    TEST_FOR_EXCEPTION_CLASS_FUNC(isFillComplete() == false, std::runtime_error, ": fill must be complete.");
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
    // we don't choose here; we have to abide by the existing graph
    const OptimizeOption oo = (this->isStorageOptimized() == true ? DoOptimizeStorage : DoNotOptimizeStorage);    
    newmat->fillComplete(this->getDomainMap(), this->getRangeMap(), oo);
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
    TEST_FOR_EXCEPTION_CLASS_FUNC( staticGraph_ == null,                                             std::logic_error, err);
    TEST_FOR_EXCEPTION_CLASS_FUNC( myGraph_ != null && myGraph_ != staticGraph_,                     std::logic_error, err);
    // if matrix is fill complete, then graph must be fill complete
    TEST_FOR_EXCEPTION_CLASS_FUNC( fillComplete_ == true && staticGraph_->isFillComplete() == false, std::logic_error, err);
    // if matrix is storage optimized, it should have a 1D allocation 
    TEST_FOR_EXCEPTION_CLASS_FUNC( isStorageOptimized() == true && values2D_ != null,                std::logic_error, err);
    // if matrix/graph are static profile, then 2D allocation should not be present
    TEST_FOR_EXCEPTION_CLASS_FUNC( getProfileType() == StaticProfile  && values2D_ != null,          std::logic_error, err);
    // if matrix/graph are dynamic profile, then 1D allocation should not be present
    TEST_FOR_EXCEPTION_CLASS_FUNC( getProfileType() == DynamicProfile && values1D_ != null,          std::logic_error, err);
    // if values are allocated and they are non-zero in number, then one of the allocations should be present
    TEST_FOR_EXCEPTION_CLASS_FUNC( staticGraph_->indicesAreAllocated() 
                        && staticGraph_->getNodeAllocationSize() > 0 && staticGraph_->getNodeNumRows() > 0
                        && values2D_ == null && values1D_ == null,                                   std::logic_error, err);
    // we can nae have both a 1D and 2D allocation
    TEST_FOR_EXCEPTION_CLASS_FUNC( values1D_ != null && values2D_ != null,                           std::logic_error, err);
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


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,Teuchos::as<size_t>(11)) + 2;
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
        out << "Global max number of entries = " << getGlobalMaxNumRowEntries() << std::endl;
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


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::checkSizes(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node> & source)
  {
    // It's not clear what kind of compatibility checks on sizes can be performed here.
    // Epetra_CrsGraph doesn't check any sizes for compatibility.

    // right now, we'll only support import/exporting between CrsMatrix<Scalar>
    // if the source dist object isn't CrsMatrix or some offspring, flag this operation as incompatible.
    try  {
      const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> & A = dynamic_cast<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &>(source);
      (void)A;
    }
    catch (...) {
      return false;
    }
    return true;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::copyAndPermute(
                          const DistObject<char, LocalOrdinal,GlobalOrdinal,Node> & source,
                          size_t numSameIDs,
                          const ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    // this should succeed, because we already tested compatibility in checkSizes()
    const std::string tfecfFuncName("copyAndPermute()");
    const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> & src_mat = dynamic_cast<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &>(source);
    TEST_FOR_EXCEPTION_CLASS_FUNC(permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
        ": permuteToLIDs and permuteFromLIDs must have the same size.");
    const bool src_is_locally_indexed = src_mat.isLocallyIndexed();

    // do numSame: copy the first numSame row from the source to *this
    // specifically, copy rows corresponding to Local Elements 0,numSame-1
    Array<GlobalOrdinal> row_indices;
    Array<Scalar>        row_values;
    LocalOrdinal mylid = 0;
    for (size_t i=0; i<numSameIDs; ++i, ++mylid) {
      // get Global ID for this row
      GlobalOrdinal gid = src_mat.getMap()->getGlobalElement(mylid);
      if (src_is_locally_indexed) {
        const size_t row_length = src_mat.getNumEntriesInGlobalRow(gid);
        row_indices.resize( row_length );
        row_values.resize( row_length );
        size_t check_row_length = 0;
        src_mat.getGlobalRowCopy(gid, row_indices(), row_values(), check_row_length);
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION_CLASS_FUNC(row_length != check_row_length, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
#endif
        insertGlobalValues( gid, row_indices(), row_values() );
      }
      else {
        ArrayView<const GlobalOrdinal> row_inds; 
        ArrayView<const Scalar>        row_vals; 
        src_mat.getGlobalRowView(gid, row_inds, row_vals);
        insertGlobalValues( gid, row_inds(), row_vals() );
      }
    }

    // handle the permuted rows.
    for (size_t p=0; p<(size_t)permuteToLIDs.size(); ++p) {
      const GlobalOrdinal  mygid =   this->getMap()->getGlobalElement(permuteToLIDs[p]);
      const GlobalOrdinal srcgid = src_mat.getMap()->getGlobalElement(permuteFromLIDs[p]);
      if (src_is_locally_indexed) {
        const size_t row_length = src_mat.getNumEntriesInGlobalRow(srcgid);
        row_indices.resize( row_length );
        row_values.resize( row_length );
        size_t check_row_length = 0;
        src_mat.getGlobalRowCopy(srcgid, row_indices(), row_values(), check_row_length);
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION_CLASS_FUNC(row_length != check_row_length, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
#endif
        insertGlobalValues( mygid, row_indices(), row_values() );
      }
      else {
        ArrayView<const GlobalOrdinal> row_inds;
        ArrayView<const Scalar>        row_vals;
        src_mat.getGlobalRowView( srcgid, row_inds, row_vals);
        insertGlobalValues( mygid, row_inds(), row_vals());
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::packAndPrepare(
                          const DistObject<char, LocalOrdinal,GlobalOrdinal,Node> & source,
                          const ArrayView<const LocalOrdinal> &exportLIDs,
                          Array<char> &exports,
                          const ArrayView<size_t> & numPacketsPerLID,
                          size_t& constantNumPackets,
                          Distributor &distor)
  {

    const std::string tfecfFuncName("packAndPrepare()");
    TEST_FOR_EXCEPTION_CLASS_FUNC(exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
        ": exportLIDs and numPacketsPerLID must have the same size.");
    // this should succeed, because we already tested compatibility in checkSizes() and performed this cast in packAndPrepare()
    const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> & src_mat = dynamic_cast<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &>(source);
    const bool src_is_locally_indexed = src_mat.isLocallyIndexed();
    constantNumPackets = 0;

    // first, set the contents of numPacketsPerLID, and accumulate a total-num-packets:
    // grab the max row size, while we're at it. may need it below.
    // Subtle: numPacketsPerLID is for byte-packets, so it needs to be multiplied
    const size_t SizeOfOrdValPair = sizeof(GlobalOrdinal)+sizeof(Scalar);
    size_t totalNumEntries = 0;
    size_t maxExpRowLength = 0;
    for (size_t i=0; i<(size_t)exportLIDs.size(); ++i) {
      GlobalOrdinal expGID = src_mat.getMap()->getGlobalElement(exportLIDs[i]);
      const size_t row_length = src_mat.getNumEntriesInGlobalRow(expGID);
      numPacketsPerLID[i] = row_length * SizeOfOrdValPair;
      totalNumEntries += row_length;
      maxExpRowLength = (row_length > maxExpRowLength ? row_length : maxExpRowLength);
    }

    // Need to do the following:
    // [inds_row0 vals_row0 inds_row1 vals_row1 ... inds_rowN vals_rowN]
    if (totalNumEntries > 0) {
      // exports is an array of char (bytes). it needs room for all of the indices and values
      const size_t totalNumBytes = totalNumEntries * SizeOfOrdValPair;
      exports.resize(totalNumBytes);

      ArrayView<char> avIndsC, avValsC;
      ArrayView<GlobalOrdinal> avInds;
      ArrayView<Scalar>        avVals;

      // now loop again and pack rows of indices into exports:
      // if global indices exist in the source, then we can use view semantics
      // otherwise, we are forced to use copy semantics (for the indices; for simplicity, we'll use them for values as well)
      size_t curOffsetInBytes = 0;
      if (src_is_locally_indexed) {
        Array<GlobalOrdinal> row_inds(maxExpRowLength);
        Array<Scalar>        row_vals(maxExpRowLength);
        for (size_t i=0; i<(size_t)exportLIDs.size(); ++i) {
          // get copy
          const GlobalOrdinal GID = src_mat.getMap()->getGlobalElement(exportLIDs[i]);
          size_t rowSize;
          src_mat.getGlobalRowCopy(GID, row_inds(), row_vals(), rowSize);
          // get export views
          avIndsC = exports(curOffsetInBytes,rowSize*sizeof(GlobalOrdinal));
          avValsC = exports(curOffsetInBytes+rowSize*sizeof(GlobalOrdinal),rowSize*sizeof(Scalar));
          avInds = av_reinterpret_cast<GlobalOrdinal>(avIndsC);
          avVals = av_reinterpret_cast<Scalar       >(avValsC);
          // copy
          std::copy( row_inds.begin(), row_inds.begin()+rowSize, avInds.begin());
          std::copy( row_vals.begin(), row_vals.begin()+rowSize, avVals.begin());
          curOffsetInBytes += SizeOfOrdValPair * rowSize;
        }
      }
      else {
        ArrayView<const GlobalOrdinal> row_inds;
        ArrayView<const Scalar>        row_vals;
        for (size_t i=0; i<(size_t)exportLIDs.size(); ++i) {
          // get view
          const GlobalOrdinal GID = src_mat.getMap()->getGlobalElement(exportLIDs[i]);
          src_mat.getGlobalRowView(GID, row_inds, row_vals);
          const size_t rowSize = (size_t)row_inds.size();
          // get export views
          avIndsC = exports(curOffsetInBytes,rowSize*sizeof(GlobalOrdinal));
          avValsC = exports(curOffsetInBytes+rowSize*sizeof(GlobalOrdinal),rowSize*sizeof(Scalar));
          avInds = av_reinterpret_cast<GlobalOrdinal>(avIndsC);
          avVals = av_reinterpret_cast<Scalar       >(avValsC);
          // copy
          std::copy( row_inds.begin(), row_inds.end(), avInds.begin());
          std::copy( row_vals.begin(), row_vals.end(), avVals.begin());
          curOffsetInBytes += SizeOfOrdValPair * rowSize;
        }
      }
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION_CLASS_FUNC(curOffsetInBytes != totalNumBytes, std::logic_error,
          ": Internal logic error. Please contact Tpetra team.");
#endif
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::unpackAndCombine(
                            const ArrayView<const LocalOrdinal> &importLIDs,
                            const ArrayView<const char> &imports,
                            const ArrayView<size_t> &numPacketsPerLID,
                            size_t constantNumPackets,
                            Distributor & /* distor */,
                            CombineMode /* CM */)
  {
    // We are not checking the value of the CombineMode input-argument.
    // Any incoming column-indices are inserted into the target graph. In this context, CombineMode values
    // of ADD vs INSERT are equivalent. What is the meaning of REPLACE for CrsGraph? If a duplicate column-index
    // is inserted, it will be compressed out when fillComplete is called.
    // NOTE: I have added a note to the Tpetra todo list to revisit this discussion. CGB, 6/18/2010

    const std::string tfecfFuncName("unpackAndCombine()");
    TEST_FOR_EXCEPTION_CLASS_FUNC(importLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
        ": importLIDs and numPacketsPerLID must have the same size.");

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
      for (size_t i=0; i<(size_t)importLIDs.size(); ++i) {
        // get row info
        const LocalOrdinal LID = importLIDs[i];
        const GlobalOrdinal myGID = this->getMap()->getGlobalElement(LID);
        const size_t rowSize = numPacketsPerLID[i] / SizeOfOrdValPair;
        //Needs to be in here in case of zero lenght rows.
        //If not, the lines following the if statement error out
        //if the row length is zero. KLN 13/06/2011
        if(rowSize == 0){
          continue;
        }
        // get import views
        avIndsC = imports(curOffsetInBytes,rowSize*sizeof(GlobalOrdinal));
        avValsC = imports(curOffsetInBytes+rowSize*sizeof(GlobalOrdinal),rowSize*sizeof(Scalar));
        avInds = av_reinterpret_cast<const GlobalOrdinal>(avIndsC);
        avVals = av_reinterpret_cast<const Scalar       >(avValsC);
        // do insert
        insertGlobalValues(myGID, avInds(), avVals());
        curOffsetInBytes += rowSize * SizeOfOrdValPair;
      }
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION_CLASS_FUNC(curOffsetInBytes != totalNumBytes, std::logic_error,
          ": Internal logic error. Please contact Tpetra team.");
#endif
    }
  }

    

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                         Deprecated methods                              //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  // DEPRECATED
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalRowView(
                                GlobalOrdinal globalRow, 
                                ArrayRCP<const GlobalOrdinal> &indices,
                                ArrayRCP<const Scalar>        &values) const 
  {
    const std::string tfecfFuncName("getGlobalRowView()");
    TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() == true, std::runtime_error, ": global indices do not exist; call getLocalRowView().");
    const LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION_CLASS_FUNC(lrow == LOT::invalid(), std::runtime_error, ": globalRow (== " << globalRow << ") does not belong to this node.");
    const RowInfo rowinfo = staticGraph_->getRowInfo(lrow);
    if (values1D_ != null && rowinfo.numEntries > 0) {
      values  =                values1D_.persistingView(rowinfo.offset1D,rowinfo.numEntries);
      indices = staticGraph_->gblInds1D_.persistingView(rowinfo.offset1D,rowinfo.numEntries);
    }
    else if (values2D_ != null && rowinfo.numEntries > 0) {
      values  =                values2D_[lrow].persistingView(0,rowinfo.numEntries);
      indices = staticGraph_->gblInds2D_[lrow].persistingView(0,rowinfo.numEntries);
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  // DEPRECATED
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalRowView(
                                LocalOrdinal localRow, 
                                ArrayRCP<const LocalOrdinal> &indices,
                                ArrayRCP<const Scalar>        &values) const 
  {
    const std::string tfecfFuncName("getLocalRowView()");
    TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed() == true, std::runtime_error, ": local indices do not exist; call getGlobalRowView().");
    TEST_FOR_EXCEPTION_CLASS_FUNC(getRowMap()->isNodeLocalElement(localRow) == false, std::runtime_error, ": localRow (== " << localRow << ") is not valid on this node.");
    const RowInfo rowinfo = staticGraph_->getRowInfo(localRow);
    if (values1D_ != null && rowinfo.numEntries > 0) {
      values  =                values1D_.persistingView(rowinfo.offset1D,rowinfo.numEntries);
      indices = staticGraph_->lclInds1D_.persistingView(rowinfo.offset1D,rowinfo.numEntries);
    }
    else if (values2D_ != null && rowinfo.numEntries > 0) {
      values  =                values2D_[localRow].persistingView(0,rowinfo.numEntries);
      indices = staticGraph_->lclInds2D_[localRow].persistingView(0,rowinfo.numEntries);
    }
    return;
  }

  /////////////////////////////////////////////////////////////////////////////
  // DEPRECATED
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::optimizeStorage() {
    // provided only for backwards compatibility
    // previous semantics required that fillComplete() had been called.
    const std::string tfecfFuncName("optimizeStorage()");
    TEST_FOR_EXCEPTION_CLASS_FUNC(isFillComplete() == false, std::runtime_error, 
        " requires that fillComplete() has already been called.");
    if (isStorageOptimized() == false) {
      resumeFill();
      fillComplete(DoOptimizeStorage);
    }
  }


} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsMatrix< SCALAR , LO , GO , NODE >;

#endif
