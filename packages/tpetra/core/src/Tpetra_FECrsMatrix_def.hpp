// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_FECRSMATRIX_DEF_HPP
#define TPETRA_FECRSMATRIX_DEF_HPP
#include <sstream>
#include "Tpetra_CrsMatrix.hpp"

namespace Tpetra {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
FECrsMatrix(const Teuchos::RCP<const fe_crs_graph_type>& graph,
            const Teuchos::RCP<Teuchos::ParameterList>& params) :
  // We want the OWNED_PLUS_SHARED graph here
  // NOTE: The casts below are terrible, but necesssary
  crs_matrix_type( graph->inactiveCrsGraph_.is_null() ? Teuchos::rcp_const_cast<crs_graph_type>(Teuchos::rcp_dynamic_cast<const crs_graph_type>(graph)) : graph->inactiveCrsGraph_,params),
  feGraph_(graph)

{
  const char tfecfFuncName[] = "FECrsMatrix(RCP<const FECrsGraph>[, RCP<ParameterList>]): ";

  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (graph.is_null (), std::runtime_error, "Input graph is null.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
    (!graph->isFillComplete (), std::runtime_error, "Input graph is not "
     "fill complete. You must call fillComplete on the graph before using "
     "it to construct a FECrsMatrix.  Note that calling resumeFill on the "
     "graph makes it not fill complete, even if you had previously called "
     "fillComplete.  In that case, you must call fillComplete on the graph "
     "again.");
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC
     ( *graph->activeCrsGraph_!= fe_crs_graph_type::FE_ACTIVE_OWNED,std::runtime_error,
      "Input graph must be in FE_ACTIVE_OWNED mode when this constructor is called.");

  bool start_owned = false;
  if (! params.is_null ()) {
    if (params->isParameter ("start owned")) {
      start_owned = params->get<bool>("start owned", start_owned);
    }
  }
  if(start_owned) {
    activeCrsMatrix_  = Teuchos::rcp(new FEWhichActive(FE_ACTIVE_OWNED));
  } else {
    activeCrsMatrix_  = Teuchos::rcp(new FEWhichActive(FE_ACTIVE_OWNED_PLUS_SHARED));
  }

  // Make an "inactive" matrix, if we need to
  if(!graph->inactiveCrsGraph_.is_null() ) {
    // We are *requiring* memory aliasing here, so we'll grab the first chunk of the Owned+Shared matrix's values array to make the
    // guy for the Owned matrix.
    inactiveCrsMatrix_ = Teuchos::rcp(new crs_matrix_type(*this,graph));
  }

  fillState_ = Teuchos::rcp(new FillState(FillState::closed));
}



template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doOwnedPlusSharedToOwned(const CombineMode CM) {
  if(!inactiveCrsMatrix_.is_null() && *activeCrsMatrix_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
    // Do a self-export in "restricted mode"
    this->doExport(*this,*feGraph_->ownedRowsImporter_,CM,true);
    inactiveCrsMatrix_->fillComplete();
  }
  crs_matrix_type::fillComplete();
}//end doOverlapToLocal


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doOwnedToOwnedPlusShared(const CombineMode /* CM */) {
  // This should be a no-op for all of our purposes
}//end doLocalToOverlap

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::switchActiveCrsMatrix() {
  if(*activeCrsMatrix_ == FE_ACTIVE_OWNED_PLUS_SHARED)
    *activeCrsMatrix_ = FE_ACTIVE_OWNED;
  else
    *activeCrsMatrix_ = FE_ACTIVE_OWNED_PLUS_SHARED;

  if(inactiveCrsMatrix_.is_null()) return;

  this->swap(*inactiveCrsMatrix_);

}//end switchActiveCrsMatrix


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::endFill() {
  if(*activeCrsMatrix_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
    doOwnedPlusSharedToOwned(Tpetra::ADD);
    switchActiveCrsMatrix();
  }
  else
    throw std::runtime_error("FECrsMatrix: Local CrsMatrix already active.  Cannot endFill()");
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::beginFill()  {
  // Note: This does not throw an error since the on construction, the FECRS is in overlap mode.  Ergo, calling beginFill(),
  // like one should expect to do in a rational universe, should not cause an error.
  if(*activeCrsMatrix_ == FE_ACTIVE_OWNED) {
    this->resumeFill();
    switchActiveCrsMatrix();
  }
  this->resumeFill();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::beginAssembly() {
  const char tfecfFuncName[] = "FECrsMatrix::beginAssembly: ";
  if (*fillState_ != FillState::closed)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot begin assembly, matrix is not in a closed state "
           << "but is currently open for "
           << (*fillState_ == FillState::open ? "assembly" : "modification");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  *fillState_ = FillState::open;
  this->beginFill();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::endAssembly() {
  const char tfecfFuncName[] = "FECrsMatrix::endAssembly: ";
  if (*fillState_ != FillState::open)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot end assembly, matrix is not open for assembly "
           << "but is currently "
           << (*fillState_ == FillState::closed ? "closed" : "open for modification");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  *fillState_ = FillState::closed;
  this->endFill();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::beginModify() {
  const char tfecfFuncName[] = "FECrsMatrix::beginModify: ";
  if (*fillState_ != FillState::closed)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot begin modifying, matrix is not in a closed state "
           << "but is currently open for "
           << (*fillState_ == FillState::open ? "assembly" : "modification");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  *fillState_ = FillState::modify;
  this->resumeFill();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::endModify() {
  const char tfecfFuncName[] = "FECrsMatrix::endModify: ";
  if (*fillState_ != FillState::modify)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot end modifying, matrix is not open to modify but is currently "
           << (*fillState_ == FillState::open ? "open for assembly" : "closed");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  *fillState_ = FillState::closed;
  this->fillComplete();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValuesImpl(
  impl_scalar_type rowVals[],
  const crs_graph_type& graph,
  const RowInfo& rowInfo,
  const GlobalOrdinal inds[],
  const impl_scalar_type newVals[],
  const LocalOrdinal numElts)
{
  const char tfecfFuncName[] = "FECrsMatrix::replaceGlobalValues: ";
  if (*fillState_ != FillState::open)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot replace global values, matrix is not open for assembly "
           << "but is currently "
           << (*fillState_ == FillState::modify ? "open for modification" : "closed");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  return crs_matrix_type::replaceGlobalValuesImpl(rowVals, graph, rowInfo, inds, newVals, numElts);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValuesImpl(
  impl_scalar_type rowVals[],
  const crs_graph_type& graph,
  const RowInfo& rowInfo,
  const LocalOrdinal inds[],
  const impl_scalar_type newVals[],
  const LocalOrdinal numElts)
{
  const char tfecfFuncName[] = "FECrsMatrix::replaceLocalValues: ";
  if (*fillState_ != FillState::open && *fillState_ != FillState::modify)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot replace local values, matrix is not open to fill/modify. "
           << "The matrix is currently closed";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  return crs_matrix_type::replaceLocalValuesImpl(rowVals, graph, rowInfo, inds, newVals, numElts);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoGlobalValuesImpl(
  impl_scalar_type rowVals[],
  const crs_graph_type& graph,
  const RowInfo& rowInfo,
  const GlobalOrdinal inds[],
  const impl_scalar_type newVals[],
  const LocalOrdinal numElts,
  const bool atomic)
{
  const char tfecfFuncName[] = "FECrsMatrix::sumIntoGlobalValues: ";
  if (*fillState_ != FillState::open)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot sum in to global values, matrix is not open for assembly. "
           << "The matrix is currently "
           << (*fillState_ == FillState::modify ? "open for modification" : "closed");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  return crs_matrix_type::sumIntoGlobalValuesImpl(
    rowVals, graph, rowInfo, inds, newVals, numElts, atomic
  );
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoLocalValuesImpl(
  impl_scalar_type rowVals[],
  const crs_graph_type& graph,
  const RowInfo& rowInfo,
  const LocalOrdinal inds[],
  const impl_scalar_type newVals[],
  const LocalOrdinal numElts,
  const bool atomic)
{
  const char tfecfFuncName[] = "FECrsMatrix::sumIntoLocalValues: ";
  if (*fillState_ != FillState::open)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot sum in to local values, matrix is not open for assembly. "
           << "The matrix is currently "
           << (*fillState_ == FillState::modify ? "open for modification" : "closed");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  return crs_matrix_type::sumIntoLocalValuesImpl(
    rowVals, graph, rowInfo, inds, newVals, numElts, atomic
  );
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertGlobalValuesImpl(
  crs_graph_type& graph,
  RowInfo& rowInfo,
  const GlobalOrdinal gblColInds[],
  const impl_scalar_type vals[],
  const size_t numInputEnt)
{
  const char tfecfFuncName[] = "FECrsMatrix::insertGlobalValues: ";
  if (*fillState_ != FillState::open)
  {
    std::ostringstream errmsg;
    errmsg << "Cannot insert global values, matrix is not open for assembly. "
           << "The matrix is currently "
           << (*fillState_ == FillState::modify ? "open for modification" : "closed");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(true, std::logic_error, errmsg.str());
  }
  return crs_matrix_type::insertGlobalValuesImpl(graph, rowInfo, gblColInds, vals, numInputEnt);
}

}  // end namespace Tpetra


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_FECRSMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  template class FECrsMatrix<SCALAR, LO, GO, NODE>;



#endif // TPETRA_FECRSMATRIX_DEF
