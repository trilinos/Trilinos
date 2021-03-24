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

#ifndef TPETRA_FECRSMATRIX_DEF_HPP
#define TPETRA_FECRSMATRIX_DEF_HPP

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

  activeCrsMatrix_     = Teuchos::rcp(new FEWhichActive(FE_ACTIVE_OWNED_PLUS_SHARED));

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
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::closed,
    std::runtime_error,
    "Cannot beginAssembly, matrix is not in a closed state"
  );
  *fillState_ = FillState::open;
  this->beginFill();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::endAssembly() {
  const char tfecfFuncName[] = "FECrsMatrix::endAssembly: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::runtime_error,
    "Cannot endAssembly, matrix is not open to fill."
  );
  *fillState_ = FillState::closed;
  this->endFill();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::beginModify() {
  const char tfecfFuncName[] = "FECrsMatrix::beginModify: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::closed,
    std::runtime_error,
    "Cannot beginModify, matrix is not in a closed state"
  );
  *fillState_ = FillState::modify;
  this->resumeFill();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::endModify() {
  const char tfecfFuncName[] = "FECrsMatrix::endModify: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::modify,
    std::runtime_error,
    "Cannot endModify, matrix is not open to modify."
  );
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
  const LocalOrdinal numElts) const
{
  const char tfecfFuncName[] = "FECrsMatrix::replaceGlobalValues: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::runtime_error,
    "Cannot replace global values, matrix is not open to fill."
  );
  return CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValuesImpl(
    rowVals, graph, rowInfo, inds, newVals, numElts
  );
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
FECrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValuesImpl(
  impl_scalar_type rowVals[],
  const crs_graph_type& graph,
  const RowInfo& rowInfo,
  const LocalOrdinal inds[],
  const impl_scalar_type newVals[],
  const LocalOrdinal numElts) const
{
  const char tfecfFuncName[] = "FECrsMatrix::replaceLocalValues: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open && *fillState_ != FillState::modify,
    std::runtime_error,
    "Cannot replace local values, matrix is not open to fill/modify."
  );
  return CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValuesImpl(
    rowVals, graph, rowInfo, inds, newVals, numElts
  );
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
  const bool atomic) const
{
  const char tfecfFuncName[] = "FECrsMatrix::sumIntoGlobalValues: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::runtime_error,
    "Cannot sum in to global values, matrix is not open to fill."
  );
  return CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoGlobalValuesImpl(
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
  const bool atomic) const
{
  const char tfecfFuncName[] = "FECrsMatrix::sumIntoLocalValues: ";
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::runtime_error,
    "Cannot sum in to local values, matrix is not open to fill."
  );
  return CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoLocalValuesImpl(
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
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    *fillState_ != FillState::open,
    std::runtime_error,
    "Cannot sum in to local values, matrix is not open to fill."
  );
  return CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::insertGlobalValuesImpl(
    graph, rowInfo, gblColInds, vals, numInputEnt
  );
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
