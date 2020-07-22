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
  typedef typename local_matrix_type::values_type values_type;

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
    values_type myvals = this->getLocalMatrix().values;

    size_t numOwnedVals = graph->getLocalGraph().entries.extent(0); // OwnedVals
    inactiveCrsMatrix_ = Teuchos::rcp(new crs_matrix_type(graph,Kokkos::subview(myvals,Kokkos::pair<size_t,size_t>(0,numOwnedVals))));
  }
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




}  // end namespace Tpetra


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_FECRSMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  template class FECrsMatrix<SCALAR, LO, GO, NODE>;



#endif // TPETRA_FECRSMATRIX_DEF
