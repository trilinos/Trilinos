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

// TODO: row-wise insertion of entries in globalAssemble() may be more efficient

#include <Kokkos_NodeHelpers.hpp>
#include <Kokkos_NodeTrace.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          size_t maxNumEntriesPerRow, 
                                          ProfileType pftype)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , lclMatrix_(rowMap->getNodeNumElements(), rowMap->getNode())
  , lclMatVec_(rowMap->getNode())
  , lclMatSolve_(rowMap->getNode())
  , constructedWithOptimizedGraph_(false)
  , fillComplete_(false)
  {
    try {
      myGraph_ = Teuchos::rcp( new CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,maxNumEntriesPerRow,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );

    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, 
                                          ProfileType pftype)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , lclMatrix_(rowMap->getNodeNumElements(), rowMap->getNode())
  , lclMatVec_(rowMap->getNode())
  , lclMatSolve_(rowMap->getNode())
  , constructedWithOptimizedGraph_(false)
  , fillComplete_(false)
  {
    try {
      myGraph_ = Teuchos::rcp( new CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,NumEntriesPerRowToAlloc,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );

    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, 
                                          size_t maxNumEntriesPerRow, 
                                          ProfileType pftype)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , lclMatrix_(rowMap->getNodeNumElements(), rowMap->getNode())
  , lclMatVec_(rowMap->getNode())
  , lclMatSolve_(rowMap->getNode())
  , constructedWithOptimizedGraph_(false)
  , fillComplete_(false)
  {
    try {
      myGraph_ = Teuchos::rcp( new CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,colMap,maxNumEntriesPerRow,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );

    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, 
                                          const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, 
                                          ProfileType pftype)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(rowMap)
  , lclMatrix_(rowMap->getNodeNumElements(), rowMap->getNode())
  , lclMatVec_(rowMap->getNode())
  , lclMatSolve_(rowMap->getNode())
  , constructedWithOptimizedGraph_(false)
  , fillComplete_(false)
  {
    try {
      myGraph_ = Teuchos::rcp( new CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,colMap,NumEntriesPerRowToAlloc,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );

    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(const Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > &graph)
  : DistObject<char, LocalOrdinal,GlobalOrdinal,Node>(graph->getRowMap())
  , staticGraph_(graph)
  , lclMatrix_(graph->getRowMap()->getNodeNumElements(), graph->getRowMap()->getNode())
  , lclMatVec_(graph->getNode())
  , lclMatSolve_(graph->getNode())
  , fillComplete_(false) 
  {
    TEST_FOR_EXCEPTION(staticGraph_ == Teuchos::null, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(graph): specified pointer is null.");
    // we prohibit the case where the graph is not yet filled
    // we do not require optimized storage, but check to ensure that optimizeStorage() isn't called between now and later.
    TEST_FOR_EXCEPTION( staticGraph_->isFillComplete() == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::CrsMatrix(graph): specified graph is not fill-complete. You must fillComplete() the graph before using it to construct a CrsMatrix.");
    constructedWithOptimizedGraph_ = staticGraph_->isStorageOptimized();

    // it is okay to create this now; this will prevent us from having to check for it on every call to apply()
    // we will use a non-owning rcp to wrap *this; this is safe as long as we do not shared sameScalarMultiplyOp_ with anyone, 
    // which would allow it to persist past the destruction of *this
    sameScalarMultiplyOp_ = createCrsMatrixMultiplyOp<Scalar>( rcp(this,false).getConst() );

    // the graph has entries, and the matrix should have entries as well, set to zero. no need or point in lazy allocating in this case.
    // first argument doesn't actually matter, due to the third
    allocateValues( CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateLocal, Teuchos::ScalarTraits<Scalar>::zero(), GraphAlreadyAllocated );

    checkInternalState();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~CrsMatrix() {
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Teuchos::Comm<int> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComm() const {
    return getCrsGraph()->getComm();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::RCP<Node>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNode() const {
    return lclMatVec_.getNode();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  ProfileType CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getProfileType() const {
    return getCrsGraph()->getProfileType();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isFillComplete() const {
    return fillComplete_; 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isFillResumed() const {
    return !fillComplete_; 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isStorageOptimized() const {
    return getCrsGraph()->isStorageOptimized();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isLocallyIndexed() const {
    return getCrsGraph()->isLocallyIndexed();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isGloballyIndexed() const {
    return getCrsGraph()->isGloballyIndexed();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasColMap() const {
    return getCrsGraph()->hasColMap();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumEntries() const {
    return getCrsGraph()->getGlobalNumEntries();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeNumEntries() const {
    return getCrsGraph()->getNodeNumEntries();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumRows() const {
    return getCrsGraph()->getGlobalNumRows(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumCols() const { 
    return getCrsGraph()->getGlobalNumCols(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeNumRows() const { 
    return getCrsGraph()->getNodeNumRows(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeNumCols() const { 
    return getCrsGraph()->getNodeNumCols(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumDiags() const { 
    return getCrsGraph()->getGlobalNumDiags(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeNumDiags() const { 
    return getCrsGraph()->getNodeNumDiags(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const { 
    return getCrsGraph()->getNumEntriesInGlobalRow(globalRow); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumEntriesInLocalRow(LocalOrdinal localRow) const { 
    return getCrsGraph()->getNumEntriesInLocalRow(localRow);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalMaxNumRowEntries() const { 
    return getCrsGraph()->getGlobalMaxNumRowEntries(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeMaxNumRowEntries() const { 
    return getCrsGraph()->getNodeMaxNumRowEntries(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getIndexBase() const { 
    return getRowMap()->getIndexBase(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRowMap() const { 
    return getCrsGraph()->getRowMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getColMap() const {
    return getCrsGraph()->getColMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getDomainMap() const { 
    return getCrsGraph()->getDomainMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRangeMap() const { 
    return getCrsGraph()->getRangeMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGraph() const { 
    if (staticGraph_ != Teuchos::null) return staticGraph_;
    return myGraph_;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getCrsGraph() const { 
    if (staticGraph_ != Teuchos::null) return staticGraph_;
    return myGraph_;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isLowerTriangular() const { 
    return getCrsGraph()->isLowerTriangular(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isUpperTriangular() const { 
    return getCrsGraph()->isUpperTriangular(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isStaticGraph() const { 
    return (staticGraph_ != Teuchos::null);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApply() const {
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
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::allocateValues(typename CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateLocalGlobal lg, Scalar alpha, GraphAllocationStatus gas) {
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    // allocate values and, optionally, ask graph to allocate indices
#ifdef HAVE_TPETRA_DEBUG
    // if the graph is already allocated, then gas should be GraphAlreadyAllocated
    // otherwise, gas should be GraphNotYetAllocated
    TEST_FOR_EXCEPTION((gas == GraphAlreadyAllocated) != graph->indicesAreAllocated(), std::logic_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): Internal logic error. Please contact Tpetra team.");
    // if the graph is unallocated, then it better be a matrix-owned graph
    TEST_FOR_EXCEPTION(graph->indicesAreAllocated() == false && myGraph_ == Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): Internal logic error. Please contact Tpetra team.");
#endif
    if (gas == GraphNotYetAllocated) {
      myGraph_->allocateIndices(lg);
    }
    const size_t nlrs = getRowMap()->getNodeNumElements(),
                  nta = graph->getNodeAllocationSize(),
           numEntries = graph->getNodeNumEntries();
    // do we even have anything to allocate?
    if (nta > 0) {
      Teuchos::RCP<Node> node = lclMatrix_.getNode();
      ////////////////////////////////////////
      if (getProfileType() == StaticProfile) {
        //
        //  STATIC ALLOCATION PROFILE
        //
        // determine how many entries to allocate and setup offsets into 1D arrays
        values1D_ = node->template allocBuffer<Scalar>(nta);
        // init values if the graph already has valid entries
        if (numEntries > 0) {
          Kokkos::ReadyBufferHelper<Node> rbh(node);
          Kokkos::InitOp<Scalar> wdp;
          wdp.alpha = alpha;
          rbh.begin();
          wdp.x   = rbh.addNonConstBuffer(values1D_);
          rbh.end();
          node->template parallel_for<Kokkos::InitOp<Scalar> >(0,nta,wdp);
          KOKKOS_NODE_TRACE("CrsMatrix::allocateValues()")
        }
        else {
          KOKKOS_NODE_TRACE("CrsMatrix::allocateValues()")
          // a simple WriteOnly view will suffice
        }
      }
      else {
        //
        //  DYNAMIC ALLOCATION PROFILE
        //
        Kokkos::InitOp<Scalar> wdp;
        wdp.alpha = alpha;
        // allocate array of buffers
        values2D_ = Teuchos::arcp< Teuchos::ArrayRCP<Scalar> >(nlrs);
        bool someRowWasInitialized = false;
        Kokkos::ReadyBufferHelper<Node> rbh(node);
        for (size_t r=0; r<nlrs; ++r) {
          // this call to getNumAllocatedEntries() is cheap for the DynamicProfile case
          const size_t ntarow = graph->getNumAllocatedEntriesInLocalRow(r),
                        nErow = graph->getNumEntriesInLocalRow(r);
          if (ntarow > 0) {
            // allocate values for this row
            values2D_[r] = node->template allocBuffer<Scalar>(ntarow);
            // initi values in parallel, if the graph already has valid entries
            if (nErow > 0) {
              rbh.begin();
              wdp.x   = rbh.addNonConstBuffer(values2D_[r]);
              rbh.end();
              node->template parallel_for<Kokkos::InitOp<Scalar> >(0,ntarow,wdp);
              someRowWasInitialized = true;
            }
          }
        }
      }
    } // num to allocate > 0
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::checkInternalState() const {
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::RCP<Node> node = getNode();
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    using Teuchos::null;
    std::string err = Teuchos::typeName(*this) + "::checkInternalState(): Likely internal logic error. Please contact Tpetra team.";
    // check the internal state of this data structure
    // this is called by numerous state-changing methods, in a debug build, to ensure that the object 
    // always remains in a valid state

    // we must have a single graph, static or matrix-owned
    TEST_FOR_EXCEPTION( staticGraph_ == null && myGraph_ == null,                                      std::logic_error, err );
    TEST_FOR_EXCEPTION( staticGraph_ != null && myGraph_ != null,                                      std::logic_error, err );
    // constructedWithOptimizedGraph_ should only be true if matrix was constructed with a graph, in which case staticGraph_ should be non-null
    TEST_FOR_EXCEPTION( constructedWithOptimizedGraph_ == true && staticGraph_ == null,                std::logic_error, err ); 
    // if matrix is fill complete, then graph must be fill complete
    TEST_FOR_EXCEPTION( fillComplete_ == true && graph->isFillComplete() == false,             std::logic_error, err );
    // if matrix is storage optimized, it should have a 1D allocation 
    TEST_FOR_EXCEPTION( isStorageOptimized() == true && values2D_ != Teuchos::null,               std::logic_error, err );
    // if matrix/graph are static profile, then 2D allocation should not be present
    TEST_FOR_EXCEPTION( graph->getProfileType() == StaticProfile  && values2D_ != Teuchos::null, std::logic_error, err );
    // if matrix/graph are dynamic profile, then 1D allocation should not be present
    TEST_FOR_EXCEPTION( graph->getProfileType() == DynamicProfile && values1D_ != Teuchos::null, std::logic_error, err );
    // if values are allocated and they are non-zero in number, then one of the allocations should be present
    TEST_FOR_EXCEPTION( graph->indicesAreAllocated() && graph->getNodeAllocationSize() > 0 && values2D_ == Teuchos::null && values1D_ == Teuchos::null,
                        std::logic_error, err );
    // we can nae have both a 1D and 2D allocation
    TEST_FOR_EXCEPTION( values1D_ != Teuchos::null && values2D_ != Teuchos::null, std::logic_error, err );
    // compare matrix allocations against graph allocations
    if (graph->indicesAreAllocated() && graph->getNodeAllocationSize() > 0) {
      if (graph->getProfileType() == StaticProfile) {
        if (graph->isLocallyIndexed()) {
          TEST_FOR_EXCEPTION( values1D_.size() != graph->pbuf_lclInds1D_.size(),  std::logic_error, err );
        }
        else {
          TEST_FOR_EXCEPTION( values1D_.size() != graph->pbuf_gblInds1D_.size(),  std::logic_error, err );
        }
      } 
      else { // graph->getProfileType() == DynamicProfile
        if (graph->isLocallyIndexed()) {
          for (size_t r=0; r < getNodeNumRows(); ++r) {
            TEST_FOR_EXCEPTION( (values2D_[r] == Teuchos::null) != (graph->pbuf_lclInds2D_[r] == Teuchos::null), std::logic_error, err );
            if (values2D_[r] != Teuchos::null && graph->pbuf_lclInds2D_[r] != Teuchos::null) {
              TEST_FOR_EXCEPTION( values2D_[r].size() != graph->pbuf_lclInds2D_[r].size(), std::logic_error, err );
            }
          }
        }
        else {
          for (size_t r=0; r < getNodeNumRows(); ++r) {
            TEST_FOR_EXCEPTION( (values2D_[r] == Teuchos::null) != (graph->pbuf_gblInds2D_[r] == Teuchos::null), std::logic_error, err );
            if (values2D_[r] != Teuchos::null && graph->pbuf_gblInds2D_[r] != Teuchos::null) {
              TEST_FOR_EXCEPTION( values2D_[r].size() != graph->pbuf_gblInds2D_[r].size(), std::logic_error, err );
            }
          }
        }
      }
    }
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::ArrayRCP<const Scalar> 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getFullView(size_t myRow, RowInfo sizeInfo) const {
#ifdef HAVE_TPETRA_DEBUG
    std::string err = Teuchos::typeName(*this) + "::getFullView(): Internal logic error. Please contact Tpetra team.";
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(myRow) == false, std::logic_error, err);
#endif
    Teuchos::ArrayRCP<const Scalar> values = Teuchos::null;
    // sizeInfo indicates the allocation size for this row, whether it has actually been allocated or not
    if (sizeInfo.allocSize > 0 && getCrsGraph()->indicesAreAllocated()) {
      Teuchos::RCP<Node> node = getNode();
      if (getCrsGraph()->getProfileType() == StaticProfile) {
        KOKKOS_NODE_TRACE("CrsMatrix::getFullView()")
        values = node->template viewBuffer<Scalar>(sizeInfo.allocSize, values1D_ + sizeInfo.offset1D);
      }
      else {  // dynamic profile
        KOKKOS_NODE_TRACE("CrsMatrix::getFullView()")
        values = node->template viewBuffer<Scalar>(sizeInfo.allocSize, values2D_[myRow]);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(getCrsGraph()->indicesAreAllocated() && getCrsGraph()->getNodeAllocationSize() > 0 && values != Teuchos::null && static_cast<size_t>(values.size()) != sizeInfo.allocSize, std::logic_error, err);
#endif
    return values;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::ArrayRCP<Scalar> 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getFullViewNonConst(size_t myRow, RowInfo sizeInfo) {
#ifdef HAVE_TPETRA_DEBUG
    std::string err = Teuchos::typeName(*this) + "::getFullViewNonConst(): Internal logic error. Please contact Tpetra team.";
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(myRow) == false, std::logic_error, err);
#endif
    Teuchos::ArrayRCP<Scalar> values = Teuchos::null;
    // sizeInfo indicates the allocation size for this row, whether it has actually been allocated or not
    if (sizeInfo.allocSize > 0 && getCrsGraph()->indicesAreAllocated()) {
      Teuchos::RCP<Node> node = getNode();
      // if there are no valid entries, then this view can be constructed WriteOnly
      Kokkos::ReadWriteOption rw = (sizeInfo.numEntries == 0 ? Kokkos::WriteOnly : Kokkos::ReadWrite);
      if (getCrsGraph()->getProfileType() == StaticProfile) {
        KOKKOS_NODE_TRACE("CrsMatrix::getFullViewNonConst()")
        values = node->template viewBufferNonConst<Scalar>(rw, sizeInfo.allocSize, values1D_ + sizeInfo.offset1D);
      }
      else {  // dynamic profile
        KOKKOS_NODE_TRACE("CrsMatrix::getFullViewNonConst()")
        values = node->template viewBufferNonConst<Scalar>(rw, sizeInfo.allocSize, values2D_[myRow]);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(getCrsGraph()->indicesAreAllocated() && getCrsGraph()->getNodeAllocationSize() > 0 && values != Teuchos::null && static_cast<size_t>(values.size()) != sizeInfo.allocSize, std::logic_error, err);
#endif
    return values;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::updateAllocation(size_t lrow, size_t allocSize) {
    using Teuchos::ArrayRCP;
    Teuchos::RCP<Node> node = getNode();
    RowInfo sizeInfo = getCrsGraph()->getRowInfo(lrow);
    // allocate a larger space for row "lrow"
    // copy any existing data from previous allocation to new allocation
    // update sizes
    // 
    // if we already have views of the data, we will create a new view and do the copy on the host
    // otherwise, don't create a view, and do the copy on the device
    // 
    if (values2D_ == Teuchos::null) {
      values2D_ = Teuchos::arcp< ArrayRCP<Scalar> >(getNodeNumRows());
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPT( getRowMap()->isNodeLocalElement(lrow) == false );
    TEST_FOR_EXCEPT( values2D_[lrow] != Teuchos::null && allocSize < static_cast<size_t>(values2D_[lrow].size()) );
    TEST_FOR_EXCEPT( allocSize == 0 );
#endif
    ArrayRCP<Scalar> old_alloc, new_row;
    old_alloc = values2D_[lrow];
    values2D_[lrow] = node->template allocBuffer<Scalar>(allocSize);
    if (sizeInfo.numEntries) {
      KOKKOS_NODE_TRACE("CrsMatrix::updateAllocation()")
      node->template copyBuffers<Scalar>(sizeInfo.numEntries,old_alloc,values2D_[lrow]);
    }
    old_alloc = Teuchos::null;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillLocalMatrix() {
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    lclMatrix_.clear();
    if (isStorageOptimized()) {
      // fill packed matrix; it is okay for values1D_ to be null; the matrix will flag itself as empty
      lclMatrix_.setPackedValues(values1D_);
    }
    else if (graph->getProfileType() == StaticProfile) {
      if (values1D_ != Teuchos::null) {
        const size_t nlrs = getNodeNumRows();
        for (size_t r=0; r < nlrs; ++r) {
          RowInfo sizeInfo = graph->getRowInfo(r);
          Teuchos::ArrayRCP<const Scalar> rowvals;
          if (sizeInfo.numEntries > 0) {
            rowvals = values1D_.persistingView(sizeInfo.offset1D, sizeInfo.numEntries);
            lclMatrix_.set2DValues(r,rowvals);
          }
        }
      }
    }
    else if (graph->getProfileType() == DynamicProfile) {
      if (values2D_ != Teuchos::null) {
        const size_t nlrs = getNodeNumRows();
        for (size_t r=0; r < nlrs; ++r) {
          RowInfo sizeInfo = graph->getRowInfo(r);
          Teuchos::ArrayRCP<const Scalar> rowvals = values2D_[r];
          if (sizeInfo.numEntries > 0) {
            rowvals = rowvals.persistingView(0,sizeInfo.numEntries);
            lclMatrix_.set2DValues(r,rowvals);
          }
        }
      }
    }

    // submit local matrix and local graph to lclMatVec_ and lclMatSolve_
    // lclMatVec_ and lclMatSolve_ are permitted to view, but we don't care whether they do or not
    lclMatVec_.clear();
    lclMatSolve_.clear();
    Teuchos::DataAccess ret;
    ret = lclMatVec_.initializeStructure( graph->lclGraph_, Teuchos::View );
    ret = lclMatVec_.initializeValues( lclMatrix_, Teuchos::View );
    if (isLowerTriangular() || isUpperTriangular()) {
      ret = lclMatSolve_.initializeStructure( graph->lclGraph_, Teuchos::View );
      ret = lclMatSolve_.initializeValues( lclMatrix_, Teuchos::View );
    }
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
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::insertLocalValues(
                                                         LocalOrdinal localRow, 
                         const Teuchos::ArrayView<const LocalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>       &values) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): matrix was constructed with static graph; cannot insert new entries.");
    TEST_FOR_EXCEPTION(isStorageOptimized() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): cannot insert new values after optimizeStorage() has been called.");
    TEST_FOR_EXCEPTION(myGraph_->isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): graph indices are global; use insertGlobalValues().");
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): cannot insert local indices without a column map; ");
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): values.size() must equal indices.size().");
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(localRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): row does not belong to this node.");
    Teuchos::Array<LocalOrdinal> finds;
    Teuchos::Array<Scalar>       fvals;
    finds.reserve(indices.size());
    fvals.reserve(values.size());
    // use column map to filter the entries:
    const Map<LocalOrdinal,GlobalOrdinal,Node> &cmap = *getColMap();
    for (size_t i=0; i < static_cast<size_t>(indices.size()); ++i) {
      if (cmap.isNodeLocalElement(indices[i])) {
        finds.push_back(indices[i]);
        fvals.push_back(values[i]);
      }
    }
    if (finds.size() > 0) {
      if (myGraph_->indicesAreAllocated() == false) {
        allocateValues(CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateLocal, ST::zero(), GraphNotYetAllocated);
      }
      //
      ArrayRCP<LocalOrdinal> rowindsview;
      ArrayRCP<Scalar>       rowvalsview;
      RowInfo sizeInfo = myGraph_->getFullLocalViewNonConst(localRow, rowindsview);
      rowvalsview = getFullViewNonConst(localRow, sizeInfo);
      const size_t newSize = sizeInfo.numEntries + finds.size();
      if (newSize > sizeInfo.allocSize) {
        TEST_FOR_EXCEPTION(myGraph_->getProfileType() == StaticProfile, std::runtime_error,
            Teuchos::typeName(*this) << "::insertLocalValues(): new indices exceed statically allocated graph structure.");
        TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
            "::insertLocalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
        // update allocation only as much as necessary
        myGraph_->updateLocalAllocation(localRow,newSize);
        updateAllocation(localRow,newSize);
        // get new views; inefficient, but acceptible in this already inefficient case
        sizeInfo = myGraph_->getFullLocalViewNonConst(localRow, rowindsview);
        rowvalsview = getFullViewNonConst(localRow, sizeInfo);
      }
      // add new indices, values to graph and matrix rows
      myGraph_->insertLocalIndicesViaView( localRow, finds, rowindsview + sizeInfo.numEntries );
      std::copy( fvals.begin(), fvals.end(), rowvalsview + sizeInfo.numEntries);
#ifdef HAVE_TPETRA_DEBUG
      {
        RowInfo sizeInfoNew = myGraph_->getRowInfo(localRow);
        TEST_FOR_EXCEPTION(sizeInfoNew.numEntries != newSize, std::logic_error,
            Teuchos::typeName(*this) << "::insertLocalValues(): Internal logic error. Please contact Tpetra team.");
      }
#endif
      // 
      rowindsview = Teuchos::null;
      rowvalsview = Teuchos::null;
      // checkInternalState();
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::insertGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>        &values) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): matrix was constructed with static graph. Cannot insert new entries.");
    TEST_FOR_EXCEPTION(myGraph_->isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): graph indices are local; use insertLocalValues().");
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): values.size() must equal indices.size().");
    const LocalOrdinal myRow = getRowMap()->getLocalElement(globalRow);
    if (myRow != LOT::invalid()) {
      // if we have a column map, use it to filter the entries.
      // only filter if this is our row.
      Teuchos::Array<GlobalOrdinal> finds_is_temporary;
      Teuchos::Array<Scalar>        fvals_is_temporary;
      Teuchos::ArrayView<const GlobalOrdinal> findices = indices;
      Teuchos::ArrayView<const Scalar       > fvalues  = values;
      if (hasColMap()) {
        // filter indices and values through the column map
        finds_is_temporary.reserve(indices.size());
        fvals_is_temporary.reserve(values.size());
        const Map<LocalOrdinal,GlobalOrdinal,Node> &cmap = *getColMap();
        for (size_t i=0; i< static_cast<size_t>(indices.size()); ++i) {
          if (cmap.isNodeGlobalElement(indices[i])) {
            finds_is_temporary.push_back(indices[i]);
            fvals_is_temporary.push_back(values[i]);
          }
        }
        findices = finds_is_temporary();
        fvalues  = fvals_is_temporary();
      }
      // add the new indices and values
      if (findices.size() > 0) {
        if (myGraph_->indicesAreAllocated() == false) {
          allocateValues(CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateGlobal, ST::zero(), GraphNotYetAllocated);
        }
        // 
        ArrayRCP<GlobalOrdinal> rowindsview;
        ArrayRCP<Scalar>        rowvalsview;
        RowInfo sizeInfo = myGraph_->getFullGlobalViewNonConst(myRow, rowindsview);
        rowvalsview = getFullViewNonConst(myRow, sizeInfo);
        const size_t newSize = sizeInfo.numEntries + findices.size();
        if (newSize > sizeInfo.allocSize) {
          TEST_FOR_EXCEPTION(myGraph_->getProfileType() == StaticProfile, std::runtime_error,
              Teuchos::typeName(*this) << "::insertGlobalValues(): new indices exceed statically allocated graph structure.");
          TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
              "::insertGlobalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
          // update allocation only as much as necessary
          myGraph_->updateGlobalAllocation(myRow,newSize);
          updateAllocation(myRow,newSize);
          // get new views; inefficient, but acceptible in this already inefficient case
          sizeInfo = myGraph_->getFullGlobalViewNonConst(myRow, rowindsview);
          rowvalsview = getFullViewNonConst(myRow, sizeInfo);
        }
        // add new indices, values to graph and matrix rows
        myGraph_->insertGlobalIndicesViaView( myRow, findices, rowindsview + sizeInfo.numEntries );
        std::copy(  fvalues.begin(),  fvalues.end(), rowvalsview + sizeInfo.numEntries );
#ifdef HAVE_TPETRA_DEBUG
        {
          RowInfo sizeInfoNew = myGraph_->getRowInfo(myRow);
          TEST_FOR_EXCEPTION(sizeInfoNew.numEntries != newSize, std::logic_error,
              Teuchos::typeName(*this) << "::insertGlobalValues(): Internal logic error. Please contact Tpetra team.");
        }
#endif
        // 
        rowindsview = Teuchos::null;
        rowvalsview = Teuchos::null;
        // checkInternalState();
      }
    }
    else {
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
      typename Teuchos::ArrayView<const Scalar       >::iterator val =  values.begin();
      nonlocals_[globalRow].reserve( nonlocals_[globalRow].size() + indices.size() );
      for (; val != values.end(); ++val, ++ind) {
        nonlocals_[globalRow].push_back(std::make_pair(*ind, *val));
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::replaceGlobalValues(      
                                        GlobalOrdinal globalRow, 
                                        const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                                        const Teuchos::ArrayView<const Scalar>        &values) {
    using Teuchos::ArrayRCP;
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const        Scalar>::iterator val = values.begin();
    LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION(lrow == LOT::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceGlobalValues(): specified global row does not belong to this processor.");
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    // 
    if (isLocallyIndexed() == true) {
      ArrayRCP<const LocalOrdinal> lindrowview;
      RowInfo sizeInfo = graph->getFullLocalView(lrow, lindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(lrow, sizeInfo);
      while (ind != indices.end()) {
        LocalOrdinal lind = getColMap()->getLocalElement(*ind);
        size_t loc = graph->findLocalIndex(lrow,lind,lindrowview);
        if (loc != STINV) {
          valrowview[loc] = (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      lindrowview = Teuchos::null;
    }
    else if (isGloballyIndexed() == true) {
      ArrayRCP<const GlobalOrdinal> gindrowview;
      RowInfo sizeInfo = graph->getFullGlobalView(lrow, gindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(lrow, sizeInfo);
      while (ind != indices.end()) {
        size_t loc = graph->findGlobalIndex(lrow,*ind,gindrowview);
        if (loc != STINV) {
          valrowview[loc] = (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      gindrowview = Teuchos::null;
    }
    //else {
    // graph indices are not allocated, i.e., are non-existant; nothing to do
    //}
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::replaceLocalValues(      
                                        LocalOrdinal localRow,
                                        const Teuchos::ArrayView<const LocalOrdinal> &indices,
                                        const Teuchos::ArrayView<const Scalar>        &values) {
    using Teuchos::ArrayRCP;
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceLocalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const LocalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const       Scalar>::iterator val = values.begin();
    bool isLocalRow = getRowMap()->isNodeLocalElement(localRow);
    TEST_FOR_EXCEPTION(isLocalRow == false, std::runtime_error,
        Teuchos::typeName(*this) << "::replaceLocalValues(): specified local row does not belong to this processor.");
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    // 
    if (isLocallyIndexed() == true) {
      ArrayRCP<const LocalOrdinal> lindrowview;
      RowInfo sizeInfo = graph->getFullLocalView(localRow, lindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(localRow, sizeInfo);
      while (ind != indices.end()) {
        LocalOrdinal lind = *ind;
        size_t loc = graph->findLocalIndex(localRow,lind,lindrowview);
        if (loc != STINV) {
          valrowview[loc] = (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      lindrowview = Teuchos::null;
    }
    else if (isGloballyIndexed() == true) {
      ArrayRCP<const GlobalOrdinal> gindrowview;
      GlobalOrdinal g_row = graph->getRowMap()->getGlobalElement(localRow);
      RowInfo sizeInfo = graph->getFullGlobalView(g_row, gindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(localRow, sizeInfo);
      while (ind != indices.end()) {
        GlobalOrdinal gind = graph->getColMap()->getGlobalElement(*ind);
        size_t loc = graph->findGlobalIndex(g_row,gind,gindrowview);
        if (loc != STINV) {
          valrowview[loc] = (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      gindrowview = Teuchos::null;
    }
    //else {
    // graph indices are not allocated, i.e., are non-existant; nothing to do
    //}
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::sumIntoGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>        &values) {
    using Teuchos::ArrayRCP;
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::sumIntoGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const        Scalar>::iterator val = values.begin();
    LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION(lrow == LOT::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::sumIntoGlobalValues(): specified global row does not belong to this processor.");
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    //
    if (isLocallyIndexed() == true) {
      ArrayRCP<const LocalOrdinal> lindrowview;
      RowInfo sizeInfo = graph->getFullLocalView(lrow, lindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(lrow, sizeInfo);
      while (ind != indices.end()) {
        LocalOrdinal lind = getColMap()->getLocalElement(*ind);
        size_t loc = graph->findLocalIndex(lrow,lind,lindrowview);
        if (loc != STINV) {
          valrowview[loc] += (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      lindrowview = Teuchos::null;
    }
    else if (isGloballyIndexed() == true) {
      ArrayRCP<const GlobalOrdinal> gindrowview;
      RowInfo sizeInfo = graph->getFullGlobalView(lrow, gindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(lrow, sizeInfo);
      while (ind != indices.end()) {
        size_t loc = graph->findGlobalIndex(lrow,*ind,gindrowview);
        if (loc != STINV) {
          valrowview[loc] += (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      gindrowview = Teuchos::null;
    }
    //else {
      // indices are not allocated; nothing to do
    //}
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalRowCopy(
                                LocalOrdinal LocalRow, 
                                const Teuchos::ArrayView<LocalOrdinal> &indices, 
                                const Teuchos::ArrayView<Scalar>       &values,
                                size_t &numEntries) const {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isGloballyIndexed()==true && hasColMap()==false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(): local indices cannot be produced.");
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,...): specified row (==" << LocalRow << ") is not valid on this node.");
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    if (graph->isLocallyIndexed()) {
      ArrayRCP<const LocalOrdinal> indrowview;
      ArrayRCP<const Scalar>       valrowview; 
      RowInfo sizeInfo = graph->getFullLocalView(LocalRow, indrowview);
      valrowview = getFullView(LocalRow, sizeInfo);
      numEntries = sizeInfo.numEntries;
      TEST_FOR_EXCEPTION(static_cast<size_t>(indices.size()) < numEntries || static_cast<size_t>(values.size()) < numEntries, std::runtime_error, 
          Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
      if (numEntries > 0) {
        std::copy( indrowview.begin(), indrowview.begin() + numEntries, indices.begin() );
        std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      }
      valrowview = Teuchos::null;
      indrowview = Teuchos::null;
    }
    else if (graph->isGloballyIndexed()) {
      ArrayRCP<const GlobalOrdinal> indrowview;
      ArrayRCP<const Scalar>        valrowview; 
      RowInfo sizeInfo = graph->getFullGlobalView(LocalRow, indrowview);
      valrowview = getFullView(LocalRow, sizeInfo);
      numEntries = sizeInfo.numEntries;
      TEST_FOR_EXCEPTION(static_cast<size_t>(indices.size()) < numEntries || static_cast<size_t>(values.size()) < numEntries, std::runtime_error, 
          Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
      if (numEntries > 0) {
        std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      }
      for (size_t j=0; j < numEntries; ++j) {
        indices[j] = getColMap()->getLocalElement(indrowview[j]);
      }
      valrowview = Teuchos::null;
      indrowview = Teuchos::null;
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above if indices are allocated
      TEST_FOR_EXCEPTION( graph->indicesAreAllocated() == true, std::logic_error, 
          Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      numEntries = 0;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalRowCopy(
                                GlobalOrdinal globalRow, 
                                const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                const Teuchos::ArrayView<Scalar>        &values,
                                size_t &numEntries) const {
    using Teuchos::ArrayRCP;
    // Only locally owned rows can be queried, otherwise complain
    LocalOrdinal myRow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION(myRow == LOT::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCopy(globalRow,...): globalRow does not belong to this node.");
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    if (graph->isGloballyIndexed()) {
      ArrayRCP<const GlobalOrdinal> indrowview;
      ArrayRCP<const Scalar>        valrowview; 
      RowInfo sizeInfo = graph->getFullGlobalView(myRow, indrowview);
      valrowview = getFullView(myRow, sizeInfo);
      numEntries = sizeInfo.numEntries;
      TEST_FOR_EXCEPTION(static_cast<size_t>(indices.size()) < numEntries || static_cast<size_t>(values.size()) < numEntries, std::runtime_error, 
          Teuchos::typeName(*this) << "::getGlobalRowCopy(GlobalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
      if (numEntries > 0) {
        std::copy( indrowview.begin(), indrowview.begin() + numEntries, indices.begin() );
        std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      }
      valrowview = Teuchos::null;
      indrowview = Teuchos::null;
    }
    else if (graph->isLocallyIndexed()) {
      ArrayRCP<const LocalOrdinal> indrowview;
      ArrayRCP<const Scalar>       valrowview; 
      RowInfo sizeInfo = graph->getFullLocalView(myRow, indrowview);
      valrowview = getFullView(myRow, sizeInfo);
      numEntries = sizeInfo.numEntries;
      TEST_FOR_EXCEPTION(static_cast<size_t>(indices.size()) < numEntries || static_cast<size_t>(values.size()) < numEntries, std::runtime_error, 
          Teuchos::typeName(*this) << "::getGlobalRowCopy(GlobalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
      if (numEntries > 0) {
        std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      }
      for (size_t j=0; j < numEntries; ++j) {
        indices[j] = getColMap()->getGlobalElement(indrowview[j]);
      }
      valrowview = Teuchos::null;
      indrowview = Teuchos::null;
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above if indices are allocated
      TEST_FOR_EXCEPTION( graph->indicesAreAllocated() == true, std::logic_error, 
          Teuchos::typeName(*this) << "::getGlobalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      numEntries = 0;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalRowView(
                                GlobalOrdinal GlobalRow, 
                                Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
                                Teuchos::ArrayRCP<const Scalar>        &values) const {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(): global indices do not exist; call getLocalRowView().");
    LocalOrdinal lrow = getRowMap()->getLocalElement(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == LOT::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    RowInfo sizeInfo = graph->getFullGlobalView(lrow,indices);
    values = getFullView(lrow,sizeInfo);
    if (sizeInfo.numEntries != sizeInfo.allocSize) {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( indices == Teuchos::null && values != Teuchos::null, std::logic_error, 
          Teuchos::typeName(*this) << "::getGlobalRowView(): Internal logic error. Please contact Tpetra team.");
#endif
      if (indices != Teuchos::null) {
        indices = indices.persistingView(0,sizeInfo.numEntries);
        if (values != Teuchos::null) {
          values  =  values.persistingView(0,sizeInfo.numEntries);
        }
      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalRowView(
                                LocalOrdinal LocalRow, 
                                Teuchos::ArrayRCP<const LocalOrdinal> &indices,
                                Teuchos::ArrayRCP<const Scalar>        &values) const {
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(): local indices do not exist; call getGlobalRowView().");
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    RowInfo sizeInfo = graph->getFullLocalView(LocalRow,indices);
    values = getFullView(LocalRow,sizeInfo);
    if (sizeInfo.numEntries != sizeInfo.allocSize) {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( indices == Teuchos::null && values != Teuchos::null, std::logic_error, 
          Teuchos::typeName(*this) << "::getLocalRowView(): Internal logic error. Please contact Tpetra team.");
#endif
      if (indices != Teuchos::null) {
        indices = indices.persistingView(0,sizeInfo.numEntries);
        if (values != Teuchos::null) {
          values  =  values.persistingView(0,sizeInfo.numEntries);
        }
      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::scale(const Scalar &alpha) {
    // scale all values in the matrix
    // it is easiest to scale all allocated values, instead of scaling only the ones with valid entries
    // however, if there are no valid entries, we can short-circuit
    // furthermore, if the values aren't allocated, we can short-circuit (unallocated values are zero, scaling to zero)
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    const size_t     nlrs = graph->getNodeNumRows(),
                 numAlloc = graph->getNodeAllocationSize(),
               numEntries = graph->getNodeNumEntries();
    if (graph->indicesAreAllocated() == false || numAlloc == 0 || numEntries == 0) {
      // do nothing
    }
    else {
      // do the scale in parallel
      Teuchos::RCP<Node> node = lclMatrix_.getNode();
      Kokkos::ReadyBufferHelper<Node> rbh(node);
      Kokkos::SingleScaleOp<Scalar> wdp;
      wdp.alpha = alpha;
      if (graph->getProfileType() == StaticProfile) {
        rbh.begin();
        wdp.x = rbh.addNonConstBuffer(values1D_);
        rbh.end();
        node->template parallel_for<Kokkos::SingleScaleOp<Scalar> >(0,numAlloc,wdp);
      }
      else if (graph->getProfileType() == DynamicProfile) {
        for (size_t row=0; row < nlrs; ++row) {
          if (values2D_[row] != Teuchos::null) {
            rbh.begin();
            wdp.x = rbh.addNonConstBuffer(values2D_[row]);
            rbh.end();
            node->template parallel_for<Kokkos::SingleScaleOp<Scalar> >(0,values2D_[row].size(),wdp);
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::setAllToScalar(const Scalar &alpha) {
    // set all values in the matrix
    // it is easiest to set all allocated values, instead of setting only the ones with valid entries
    // however, if there are no entries, we can short-circuit
    // this method is equivalent replacing all valid entries with alpha.
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    const size_t     nlrs = graph->getNodeNumRows(),
                 numAlloc = graph->getNodeAllocationSize(),
               numEntries = graph->getNodeNumEntries();
    if (graph->indicesAreAllocated() == false || numAlloc == 0 || numEntries == 0) {
      // do nothing
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( graph->indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::setAllToScalar(): internal logic error. Please contact Tpetra team.");
#endif
      // set the values in parallel
      Teuchos::RCP<Node> node = lclMatrix_.getNode();
      Kokkos::ReadyBufferHelper<Node> rbh(node);
      Kokkos::InitOp<Scalar> wdp;
      wdp.alpha = alpha;
      if (graph->getProfileType() == StaticProfile) {
        rbh.begin();
        wdp.x = rbh.addNonConstBuffer(values1D_);
        rbh.end();
        node->template parallel_for<Kokkos::InitOp<Scalar> >(0,numAlloc,wdp);
      }
      else if (graph->getProfileType() == DynamicProfile) {
        for (size_t row=0; row < nlrs; ++row) {
          if (values2D_[row] != Teuchos::null) {
            rbh.begin();
            wdp.x = rbh.addNonConstBuffer(values2D_[row]);
            rbh.end();
            node->template parallel_for<Kokkos::InitOp<Scalar> >(0,values2D_[row].size(),wdp);
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &dvec) const {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << ": cannot call getLocalDiagCopy() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(dvec.getMap()->isSameAs(*getRowMap()) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalDiagCopy(dvec): dvec must have the same map as the CrsMatrix.");
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
#ifdef HAVE_TPETRA_DEBUG
    size_t numDiagFound = 0;
#endif
    const size_t nlrs = getNodeNumRows();
    Teuchos::ArrayRCP<Scalar> vecView = dvec.get1dViewNonConst();
    for (size_t r=0; r < nlrs; ++r) {
      vecView[r] = Teuchos::ScalarTraits<Scalar>::zero();
      GlobalOrdinal rgid = getRowMap()->getGlobalElement(r);
      if (getColMap()->isNodeGlobalElement(rgid)) {
        LocalOrdinal rlid = getColMap()->getLocalElement(rgid);
        Teuchos::ArrayRCP<const LocalOrdinal> inds;
        RowInfo sizeInfo = graph->getFullLocalView(r,inds);
        if (sizeInfo.numEntries > 0) {
          const size_t j = graph->findLocalIndex(r, rlid, inds);
          if (j != STINV) {
            Teuchos::ArrayRCP<const Scalar> vals;
            vals = getFullView(r,sizeInfo);
            vecView[r] = vals[j];
            vals = Teuchos::null;
#ifdef HAVE_TPETRA_DEBUG
            ++numDiagFound;
#endif
          }
        }
        inds = Teuchos::null;
      }
    }
    vecView = Teuchos::null;
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(numDiagFound != getNodeNumDiags(), std::logic_error, 
        "CrsMatrix::getLocalDiagCopy(): logic error. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::globalAssemble() {
    using Teuchos::arcp;
    using Teuchos::OrdinalTraits;
    using Teuchos::Array;
    using Teuchos::SerialDenseMatrix;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    using Teuchos::outArg;
    using std::pair;
    using std::make_pair;
    using Teuchos::tuple;
    typedef typename std::map<GlobalOrdinal,Teuchos::Array<pair<GlobalOrdinal,Scalar> > >::const_iterator NLITER;
    typedef typename Teuchos::Array<pair<GlobalOrdinal,Scalar> >::const_iterator NLRITER;
    const int numImages = getComm()->getSize();
    const int myImageID = getComm()->getRank();
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
        TEST_FOR_EXCEPTION(gblerror, std::runtime_error,
            Teuchos::typeName(*this) << "::globalAssemble(): non-local entries correspond to invalid rows.");
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
        IdsAndRows.push_back(make_pair<int,GlobalOrdinal>(*id,*nlr));
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
        TEST_FOR_EXCEPTION(sendIDs[numSends] != id, std::logic_error, Teuchos::typeName(*this) << "::globalAssemble(): internal logic error. Contact Tpetra team.");
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
    TEST_FOR_EXCEPTION(Teuchos::as<typename Array<int>::size_type>(numSends) != sendIDs.size(), std::logic_error, Teuchos::typeName(*this) << "::globalAssemble(): internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    nonlocals_.clear();

    ////////////////////////////////////////////////////////////////////////////////////// 
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    ////////////////////////////////////////////////////////////////////////////////////// 
    // perform non-blocking sends: send sizes to our recipients
    Array<Teuchos::RCP<Teuchos::CommRequest> > sendRequests;
    for (size_t s=0; s < numSends ; ++s) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back( Teuchos::isend<int,size_t>(*getComm(),Teuchos::rcpFromRef(sendSizes[s]),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<Teuchos::RCP<Teuchos::CommRequest> > recvRequests;
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
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > sendBuffers(numSends,Teuchos::null);
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
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > recvBuffers(numRecvs,Teuchos::null);
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
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillComplete(OptimizeOption os) {
    fillComplete(getRowMap(),getRowMap(),os);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillComplete(
                                            const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, 
                                            const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, 
                                            OptimizeOption os) {

    if (getCrsGraph()->indicesAreAllocated() == false) {
      // must allocate global, because we have no column map
      allocateValues( CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateGlobal, Teuchos::ScalarTraits<Scalar>::zero(), GraphNotYetAllocated );
    }

    if (getComm()->getSize() > 1) {
      globalAssemble();
    }
    else {
      TEST_FOR_EXCEPTION(nonlocals_.size() > 0, std::runtime_error,
          Teuchos::typeName(*this) << "::fillComplete(): cannot have non-local entries on a serial run. Invalid entry was submitted to the CrsMatrix.");
    }

    // if the matrix was constructed with an un-optimized graph, it must not have been optimized in the interim.
    // this would likely have wrecked the fragile dance between graph and matrix.
    if (isStaticGraph()) {
      if (constructedWithOptimizedGraph_ == false && staticGraph_->isStorageOptimized() == true) {
        TPETRA_ABUSE_WARNING(true,std::runtime_error,
            "::fillComplete(): optimizeStorage() has been called on graph since matrix was constructed.");
        return;
      }
    }
    
    // if we're not allowed to change a static graph, then we can't call optimizeStorage() on it.
    // then we can't call fillComplete() on the matrix with DoOptimizeStorage
    // throw a warning. this is an unfortunate late evaluation; however, we couldn't know when we received the graph that the user would try to 
    // optimize the storage later on.
    if (os == DoOptimizeStorage && isStaticGraph() && staticGraph_->isStorageOptimized() == false) {
      TPETRA_ABUSE_WARNING(true,std::runtime_error,
          "::fillComplete(): requested optimized storage, but static graph does not have optimized storage. Ignoring request to optimize storage.");
      os = DoNotOptimizeStorage;
    }

    if (isStaticGraph() == false) {
      myGraph_->makeIndicesLocal(domainMap,rangeMap);
      sortEntries();
      mergeRedundantEntries();
      // can't optimize graph storage before optimizing our storage; therefore, do not optimize storage yet.
      myGraph_->fillComplete(domainMap,rangeMap,DoNotOptimizeStorage);
    }

    fillComplete_ = true;

    if (os == DoOptimizeStorage && isStorageOptimized() == false) {
      // this will also:
      // * call optimizeStorage() on the graph as well, if isStaticGraph() == false
      //   this will call fillLocalGraph() on the graph
      // * call fillLocalMatrix()
      //   this will initialize the local mat-vec and (if appropriate) the local solve
      optimizeStorage();
    }
    else { 
      // local graph already filled.
      // fill the local matrix.
      fillLocalMatrix();
    }
  
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::sortEntries() {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(myGraph_->isGloballyIndexed() == true, std::logic_error,
        Teuchos::typeName(*this) << "::sortEntries(): sortEntries() must be called after indices are transformed to local.\n"
                                 << "Likely internal logic error. Please contact Tpetra team.");
    if (myGraph_->isSorted()) return;
    if (myGraph_->indicesAreAllocated() && myGraph_->getNodeAllocationSize() > 0) {
      const size_t nlrs = getNodeNumRows();
      for (size_t r=0; r < nlrs; ++r) {
        // TODO: This is slightly inefficient, because it may query pbuf_rowOffsets_ repeatadly. 
        //       However, it is very simple code. Consider rewriting it.
        Teuchos::ArrayRCP<LocalOrdinal> inds;
        Teuchos::ArrayRCP<Scalar>       vals;
        RowInfo sizeInfo = myGraph_->getFullLocalViewNonConst(r, inds);
        vals = getFullViewNonConst(r, sizeInfo);
        sort2(inds.begin(), inds.begin() + sizeInfo.numEntries, vals);
      }
    }
    myGraph_->setSorted(true);  // we just sorted them
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::mergeRedundantEntries() {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(myGraph_->isSorted() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::mergeRedundantEntries() cannot be called before indices are sorted.\n"
                                 << "Likely interal logic error. Please contact Tpetra team.");
    if (myGraph_->notRedundant()) return;
    if (myGraph_->indicesAreAllocated() && myGraph_->getNodeAllocationSize() > 0) {
      for (size_t r=0; r<getNodeNumRows(); ++r) {
        // TODO: This is slightly inefficient, because it may query pbuf_rowOffsets_ repeatadly. 
        //       However, it is very simple code. Consider rewriting it.
        Teuchos::ArrayRCP<LocalOrdinal> inds;
        Teuchos::ArrayRCP<Scalar>       vals;
        RowInfo sizeInfo = myGraph_->getFullLocalViewNonConst(r, inds);
        vals = getFullViewNonConst(r, sizeInfo);
        if (sizeInfo.numEntries > 0) {
          size_t curEntry = 0;
          Scalar curValue = vals[curEntry];
          for (size_t k=1; k < sizeInfo.numEntries; ++k) {
            if (inds[k] == inds[k-1]) {
              curValue += vals[k];
            }
            else {
              vals[curEntry++] = curValue;
              curValue = vals[k];
            }
          }
          vals[curEntry] = curValue;
        }
        vals = Teuchos::null;
        inds = Teuchos::null;
      }
      myGraph_->removeRedundantIndices();
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::optimizeStorage() {
    // optimizeStorage will perform two functions:
    // 1) create a single allocation of memory
    // 2) pack data in that allocation
    // if getProfileType() == StaticProfile, then 1) has already been done
#ifdef HAVE_TPETRA_DEBUG
    std::string err = Teuchos::typeName(*this) + "::optimizeStorage(): Internal logic error. Please contact Tpetra team.";
#endif
    if (isStorageOptimized() == true) return;
    // if storage is not optimized and we have a static graph, we will not be able to optimize the storage. throw an exception.
    TEST_FOR_EXCEPTION( isStaticGraph(), std::logic_error, 
        Teuchos::typeName(*this) << "::optimizeStorage(): Cannot optimize storage with a static graph. Possible internal logic error. Please contact Tpetra team.");
    TEST_FOR_EXCEPTION(isFillComplete() == false || myGraph_->isSorted() == false || myGraph_->notRedundant() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): fillComplete() must be called before optimizeStorage().");
    // 
    Teuchos::RCP<Node> node = lclMatrix_.getNode();
    const size_t       nlrs = getNodeNumRows(),
                 numEntries = myGraph_->getNodeNumEntries();
    if (nlrs > 0) {
      if (myGraph_->getProfileType() == DynamicProfile) {
        // allocate a single memory block
        if (numEntries > 0) {
          // numEntries > 0  =>  allocSize > 0  =>  values2D_ != Teuchos::null
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION( values2D_ == Teuchos::null, std::logic_error, err);
#endif
          values1D_ = node->template allocBuffer<Scalar>(numEntries);
          size_t sofar = 0;
          for (size_t row=0; row<nlrs; ++row) {
            RowInfo sizeInfo = myGraph_->getRowInfo(row);
            if (sizeInfo.numEntries > 0) {
              KOKKOS_NODE_TRACE("CrsMatrix::optimizeStorage()")
              node->template copyBuffers<Scalar>( sizeInfo.numEntries, values2D_[row], values1D_ + sofar );
              values2D_[row] = Teuchos::null;
            }
            sofar += sizeInfo.numEntries;
          }
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(sofar != numEntries, std::logic_error, err);
#endif
        }
        // done with 2D allocation
        values2D_ = Teuchos::null;
      }
      else {
        // storage is already allocated; just need to pack
        if (numEntries > 0) {
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION( values1D_ == Teuchos::null, std::logic_error, err);
#endif
          size_t sofar = 0;
          for (size_t row=0; row<nlrs; ++row) {
            RowInfo sizeInfo = myGraph_->getRowInfo(row);
            if (sizeInfo.numEntries > 0 && sizeInfo.offset1D != sofar) {
              KOKKOS_NODE_TRACE("CrsMatrix::optimizeStorage()")
              node->template copyBuffers<Scalar>( sizeInfo.numEntries, 
                                                  values1D_ + sizeInfo.offset1D,
                                                  values1D_ + sofar );
            }
            sofar += sizeInfo.numEntries;
          }
          values1D_ = values1D_.persistingView(0,sofar);
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(sofar != numEntries, std::logic_error, err);
#endif
        }
      }
    }

    // now we can optimize the graph storage; this process is symmetric w.r.t. the process undertaken in the matrix
    myGraph_->optimizeStorage();

    // if empty, then release buffers
    // this mirrors similar code in CrsGraph::optimizeStorage, and is necessary so that 
    // both the local matrix and local graph are marked empty
    if (myGraph_->getNodeAllocationSize() == 0) {
      values1D_ = Teuchos::null;
    }

    // local graph was filled during myGraph_->optimizeStorage()
    fillLocalMatrix(); 
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(
                                        const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                                        MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                                        Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
    TEST_FOR_EXCEPTION( isFillComplete() == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::apply(): cannot call apply() until fillComplete() has been called.");
    sameScalarMultiplyOp_->apply(X,Y,mode,alpha,beta);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  template <class DomainScalar, class RangeScalar>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::multiply(
                                        const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                                              MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                                              Teuchos::ETransp mode, RangeScalar alpha, RangeScalar beta) const {
    typedef Teuchos::ScalarTraits<RangeScalar> RST;
    const Kokkos::MultiVector<DomainScalar,Node> *lclX = &X.getLocalMV();
    Kokkos::MultiVector<RangeScalar,Node>        *lclY = &Y.getLocalMVNonConst();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << ": cannot call multiply() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(X,Y): X and Y must have the same number of vectors.");
    TEST_FOR_EXCEPTION(X.isConstantStride() == false || Y.isConstantStride() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(X,Y): X and Y must be constant stride.");
    TEST_FOR_EXCEPTION(lclX==lclY, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(X,Y): X and Y cannot share data.");
#endif
    //
    // Call the matvec
    if (beta == RST::zero()) {
      // Y = alpha*op(M)*X with overwrite semantics
      lclMatVec_.template multiply<DomainScalar,RangeScalar>(mode, alpha, *lclX, *lclY);
    }
    else {
      // Y = alpha*op(M) + beta*Y
      lclMatVec_.template multiply<DomainScalar,RangeScalar>(mode, alpha, *lclX, beta, *lclY);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  template <class DomainScalar, class RangeScalar>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::solve(
                                    const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>  &Y, 
                                          MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                                          Teuchos::ETransp mode) const {
    const Kokkos::MultiVector<RangeScalar,Node> *lclY = &Y.getLocalMV();
    Kokkos::MultiVector<DomainScalar,Node>      *lclX = &X.getLocalMVNonConst();
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << ": cannot call solve() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::solve(X,Y): X and Y must have the same number of vectors.");
    TEST_FOR_EXCEPTION(X.isConstantStride() == false || Y.isConstantStride() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(X,Y): X and Y must be constant stride.");
    TEST_FOR_EXCEPTION(isUpperTriangular() == false && isLowerTriangular() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): can only solve() triangular matrices.");
    TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
        Teuchos::typeName(*this) << "::solve() does not currently support transposed solve for complex scalar types.");
#endif
    //
    // Call the solve
    Teuchos::EDiag diag = ( getNodeNumDiags() < getNodeNumRows() ? Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG );
    if (mode == Teuchos::NO_TRANS) {
      if (isUpperTriangular()) {
        lclMatSolve_.template solve<DomainScalar,RangeScalar>(Teuchos::NO_TRANS, Teuchos::UPPER_TRI, diag, *lclY, *lclX);
      }
      else {
        lclMatSolve_.template solve<DomainScalar,RangeScalar>(Teuchos::NO_TRANS, Teuchos::LOWER_TRI, diag, *lclY, *lclX);
      }
    }
    else {
      if (isUpperTriangular()) {
        lclMatSolve_.template solve<DomainScalar,RangeScalar>(Teuchos::CONJ_TRANS, Teuchos::UPPER_TRI, diag, *lclY, *lclX);
      }
      else {
        lclMatSolve_.template solve<DomainScalar,RangeScalar>(Teuchos::CONJ_TRANS, Teuchos::LOWER_TRI, diag, *lclY, *lclX);
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  std::string CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::description() const {
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
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
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
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,11) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    // 
    // for medium and higher, print constituent objects at specified verbLevel
    Teuchos::RCP< const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph = getCrsGraph();
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
        if (getColMap() != Teuchos::null) {
          if (getColMap() == getRowMap()) {
            if (myImageID == 0) out << "\nColumn map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nColumn map: " << std::endl;
            getColMap()->describe(out,vl);
          }
        }
        if (getDomainMap() != Teuchos::null) {
          if (getDomainMap() == getRowMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else if (getDomainMap() == getColMap()) {
            if (myImageID == 0) out << "\nDomain map is row map.";
          }
          else {
            if (myImageID == 0) out << "\nDomain map: " << std::endl;
            getDomainMap()->describe(out,vl);
          }
        }
        if (getRangeMap() != Teuchos::null) {
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
            if (graph->indicesAreAllocated() == false) {
              out << "Node not allocated" << std::endl;
            }
            else {
              out << "Node number of allocated entries = " << graph->getNodeAllocationSize() << std::endl;
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
                  Teuchos::ArrayRCP<const GlobalOrdinal> rowinds;
                  Teuchos::ArrayRCP<const Scalar> rowvals;
                  getGlobalRowView(gid,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << " (" << rowinds[j]
                        << ", " << rowvals[j]
                        << ") ";
                  }
                }
                else if (isLocallyIndexed()) {
                  Teuchos::ArrayRCP<const LocalOrdinal> rowinds;
                  Teuchos::ArrayRCP<const Scalar> rowvals;
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
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::checkSizes(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node> & source)
  {
    // It's not clear what kind of compatibility checks on sizes can be performed here.
    // Epetra_CrsGraph doesn't check any sizes for compatibility.

    // right now, we'll only support import/exporting between CrsMatrix<Scalar>
    // if the source dist object isn't CrsMatrix or some offspring, flag this operation as incompatible.
    try  {
      const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> & A = dynamic_cast<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> &>(source);
      (void)A;
    }
    catch (...) {
      return false;
    }
    return true;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::copyAndPermute(
                          const DistObject<char, LocalOrdinal,GlobalOrdinal,Node> & source,
                          size_t numSameIDs,
                          const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    // this should succeed, because we already tested compatibility in checkSizes()
    const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> & src_mat = dynamic_cast<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> &>(source);
    TEST_FOR_EXCEPTION(permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::copyAndPermute: permuteToLIDs and permuteFromLIDs must have the same size.");
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
        TEST_FOR_EXCEPTION(row_length != check_row_length, std::logic_error,
            Teuchos::typeName(*this) << "::copyAndPermute(): Internal logic error. Please contact Tpetra team.");
#endif
        insertGlobalValues( gid, row_indices(), row_values() );
      }
      else {
        ArrayRCP<const GlobalOrdinal> row_inds; 
        ArrayRCP<const Scalar>        row_vals; 
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
        TEST_FOR_EXCEPTION(row_length != check_row_length, std::logic_error,
            Teuchos::typeName(*this) << "::copyAndPermute(): Internal logic error. Please contact Tpetra team.");
#endif
        insertGlobalValues( mygid, row_indices(), row_values() );
      }
      else {
        ArrayRCP<const GlobalOrdinal> row_inds;
        ArrayRCP<const Scalar>        row_vals;
        src_mat.getGlobalRowView( srcgid, row_inds, row_vals);
        insertGlobalValues( mygid, row_inds(), row_vals());
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::packAndPrepare(
                          const DistObject<char, LocalOrdinal,GlobalOrdinal,Node> & source,
                          const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                          Teuchos::Array<char> &exports,
                          const Teuchos::ArrayView<size_t> & numPacketsPerLID,
                          size_t& constantNumPackets,
                          Distributor &distor)
  {
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;

    TEST_FOR_EXCEPTION(exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::packAndPrepare: exportLIDs and numPacketsPerLID must have the same size.");
    // this should succeed, because we already tested compatibility in checkSizes() and performed this cast in packAndPrepare()
    const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> & src_mat = dynamic_cast<const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> &>(source);
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
          avInds = Teuchos::av_reinterpret_cast<GlobalOrdinal>(avIndsC);
          avVals = Teuchos::av_reinterpret_cast<Scalar       >(avValsC);
          // copy
          std::copy( row_inds.begin(), row_inds.begin()+rowSize, avInds.begin());
          std::copy( row_vals.begin(), row_vals.begin()+rowSize, avVals.begin());
          curOffsetInBytes += SizeOfOrdValPair * rowSize;
        }
      }
      else {
        ArrayRCP<const GlobalOrdinal> row_inds;
        ArrayRCP<const Scalar>        row_vals;
        for (size_t i=0; i<(size_t)exportLIDs.size(); ++i) {
          // get view
          const GlobalOrdinal GID = src_mat.getMap()->getGlobalElement(exportLIDs[i]);
          src_mat.getGlobalRowView(GID, row_inds, row_vals);
          const size_t rowSize = (size_t)row_inds.size();
          // get export views
          avIndsC = exports(curOffsetInBytes,rowSize*sizeof(GlobalOrdinal));
          avValsC = exports(curOffsetInBytes+rowSize*sizeof(GlobalOrdinal),rowSize*sizeof(Scalar));
          avInds = Teuchos::av_reinterpret_cast<GlobalOrdinal>(avIndsC);
          avVals = Teuchos::av_reinterpret_cast<Scalar       >(avValsC);
          // copy
          std::copy( row_inds.begin(), row_inds.end(), avInds.begin());
          std::copy( row_vals.begin(), row_vals.end(), avVals.begin());
          curOffsetInBytes += SizeOfOrdValPair * rowSize;
        }
      }
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(curOffsetInBytes != totalNumBytes, std::logic_error,
          Teuchos::typeName(*this) << "::packAndPrepare(): Internal logic error. Please contact Tpetra team.");
#endif
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::unpackAndCombine(
                            const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                            const Teuchos::ArrayView<const char> &imports,
                            const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                            size_t constantNumPackets,
                            Distributor & /* distor */,
                            CombineMode /* CM */)
  {
    using Teuchos::ArrayView;
    // We are not checking the value of the CombineMode input-argument.
    // Any incoming column-indices are inserted into the target graph. In this context, CombineMode values
    // of ADD vs INSERT are equivalent. What is the meaning of REPLACE for CrsGraph? If a duplicate column-index
    // is inserted, it will be compressed out when fillComplete is called.
    // NOTE: I have added a note to the Tpetra todo list to revisit this discussion. CGB, 6/18/2010

    TEST_FOR_EXCEPTION(importLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::unpackAndCombine: importLIDs and numPacketsPerLID must have the same size.");

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
        // get import views
        avIndsC = imports(curOffsetInBytes,rowSize*sizeof(GlobalOrdinal));
        avValsC = imports(curOffsetInBytes+rowSize*sizeof(GlobalOrdinal),rowSize*sizeof(Scalar));
        avInds = Teuchos::av_reinterpret_cast<const GlobalOrdinal>(avIndsC);
        avVals = Teuchos::av_reinterpret_cast<const Scalar       >(avValsC);
        // do insert
        insertGlobalValues(myGID, avInds(), avVals());
        curOffsetInBytes += rowSize * SizeOfOrdValPair;
      }
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(curOffsetInBytes != totalNumBytes, std::logic_error,
          Teuchos::typeName(*this) << "::packAndPrepare(): Internal logic error. Please contact Tpetra team.");
#endif
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::createViews() const {
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::createViewsNonConst(Kokkos::ReadWriteOption rwo) {
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::releaseViews() const {
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
