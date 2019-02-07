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

#ifndef TPETRA_FEMULTIVECTOR_DEF_HPP
#define TPETRA_FEMULTIVECTOR_DEF_HPP

/// \file Tpetra_MultiVector_def.hpp
/// \brief Definition of the Tpetra::MultiVector class

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Details_Behavior.hpp"


namespace Tpetra {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
FEMultiVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > & map,
              const Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >& importer,
              const size_t numVecs,
              const bool zeroOut):
  base_type(importer.is_null()? map:importer->getTargetMap(),numVecs,zeroOut),
  importer_(importer) {
  const char tfecfFuncName[] = "FEMultiVector::FEMultiVector(): ";

  activeMultiVector_ = Teuchos::rcp(new FEWhichActive(FE_ACTIVE_OWNED_PLUS_SHARED));

  // Sanity check the importer
  if(!importer_.is_null() && !importer_->getSourceMap()->isSameAs(*map)) {
    throw std::runtime_error("FEMultiVector: 'map' must match 'importer->getSourceMap()' if importer is provided");
  }

  if(!importer_.is_null()) {
    // The locallyFitted check is debug mode only since it is more expensive
    const bool debug = ::Tpetra::Details::Behavior::debug ();
    if(debug) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( !importer->getTargetMap()->isLocallyFitted(*importer->getSourceMap()),
                                             std::runtime_error,"importer->getTargetMap() must be locally fitted to importer->getSourceMap()");
      
    }

    // Memory aliasing is required for FEMultiVector
    inactiveMultiVector_ = Teuchos::rcp(new base_type(importer_->getSourceMap(),Kokkos::subview(this->view_,Kokkos::pair<size_t,size_t>(0,map->getNodeNumElements()),Kokkos::ALL)));
  }

}// end constructor




template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceMap (const Teuchos::RCP<const map_type>& /* newMap */) {
  throw std::runtime_error("Tpetra::FEMultiVector::replaceMap() is not implemented");
}// end replaceMap

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doOwnedPlusSharedToOwned(const CombineMode CM) {
  if(!importer_.is_null() && *activeMultiVector_ == FE_ACTIVE_OWNED_PLUS_SHARED) {
    inactiveMultiVector_->doExport(*this,*importer_,CM);
  }
}//end doTargetToSource


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doOwnedToOwnedPlusShared(const CombineMode CM) {
  if(!importer_.is_null() && *activeMultiVector_ == FE_ACTIVE_OWNED) {
    inactiveMultiVector_->doImport(*this,*importer_,CM);
  }
}//end doTargetToSource

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FEMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::switchActiveMultiVector() {
  if(*activeMultiVector_ == FE_ACTIVE_OWNED_PLUS_SHARED)
    *activeMultiVector_ = FE_ACTIVE_OWNED;
  else
    *activeMultiVector_ = FE_ACTIVE_OWNED_PLUS_SHARED;

  if(importer_.is_null()) return;

  // Use MultiVector's swap routine here
  this->swap(*inactiveMultiVector_);

}//end switchActiveMultiVector

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_FEMULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  template class FEMultiVector< SCALAR , LO , GO , NODE >;

#endif // TPETRA_FEMULTIVECTOR_DEF_HPP


