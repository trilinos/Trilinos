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

#ifndef TPETRA_FECRSGRAPH_DEF_HPP
#define TPETRA_FECRSGRAPH_DEF_HPP

#include "Tpetra_CrsGraph.hpp"



namespace Tpetra {

template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph(const Teuchos::RCP<const map_type> & ownedRowMap,
           const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap, 
           const size_t maxNumEntriesPerRow,
           const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
           const Teuchos::RCP<const map_type> & domainMap,
           const Teuchos::RCP<const map_type> & rangeMap,
           const Teuchos::RCP<Teuchos::ParameterList>& params): 
  crs_graph_type(ownedPlusSharedRowMap, maxNumEntriesPerRow, StaticProfile, params),
  importer_(ownedPlusSharedToOwnedimporter),
  domainMap_(domainMap),
  rangeMap_(rangeMap)
{
  activeCrsGraph_     = Teuchos::rcp(new FEWhichActive(FE_ACTIVE_OWNED_PLUS_SHARED));
  // FIXME: Finish
  // NOTE: Getting the memory aliasing correct here will be rather tricky
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
FECrsGraph (const Teuchos::RCP<const map_type> & ownedRowMap,
            const Teuchos::RCP<const map_type> & ownedPlusSharedRowMap, 
            const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
            const Teuchos::RCP<const import_type> & ownedPlusSharedToOwnedimporter,
            const Teuchos::RCP<const map_type> & domainMap,
            const Teuchos::RCP<const map_type> & rangeMap,
            const Teuchos::RCP<Teuchos::ParameterList>& params):
  crs_graph_type( ownedPlusSharedRowMap, numEntPerRow, StaticProfile, params),
  importer_(ownedPlusSharedToOwnedimporter),
  domainMap_(domainMap),
  rangeMap_(rangeMap)
{
  activeCrsGraph_     = Teuchos::rcp(new FEWhichActive(FE_ACTIVE_OWNED_PLUS_SHARED));
  // FIXME: Finish
  // NOTE: Getting the memory aliasing correct here will be rather tricky
}




template<class LocalOrdinal, class GlobalOrdinal, class Node>
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>&
FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::
operator=(const FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>& rhs)
{
  return *this;
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doOwnedPlusSharedToOwned(const CombineMode CM) {
  if(!inactiveCrsGraph_.is_null() && *activeCrsGraph_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
    // FIXME: Insert Tim's Magic Here
    //    inactiveCrsGraph_->doExport(*this,*importer_,CM);
  }
}//end doOverlapToLocal


template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::doOwnedToOwnedPlusShared(const CombineMode CM) {
  // This should be a no-op for all of our purposes
}//end doLocalToOverlap

template<class LocalOrdinal, class GlobalOrdinal, class Node>
void FECrsGraph<LocalOrdinal, GlobalOrdinal, Node>::switchActiveCrsGraph() {
  if(*activeCrsGraph_ == FE_ACTIVE_OWNED_PLUS_SHARED)
    *activeCrsGraph_ = FE_ACTIVE_OWNED;
  else
    *activeCrsGraph_ = FE_ACTIVE_OWNED_PLUS_SHARED;

  if(inactiveCrsGraph_.is_null()) return;

  // FIXME: Use CrsGraph's swap routine here once it is implemented!!!
  // 
  //  this->swap(*inactiveCrsGraph_);

}//end switchActiveCrsGraph


}  // end namespace Tpetra


//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
#define TPETRA_FECRSGRAPH_GRAPH_INSTANT(LO,GO,NODE) \
  template class FECrsGraph<LO, GO, NODE>;



#endif // TPETRA_FECRSGRAPH_DEF
