// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_GATHER_TANGENT_BLOCKED_EPETRA_IMPL_HPP
#define PANZER_GATHER_TANGENT_BLOCKED_EPETRA_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Epetra_Map.h"

template<typename EvalT,typename TRAITS,typename LO,typename GO>
panzer::GatherTangent_BlockedEpetra<EvalT, TRAITS,LO,GO>::
GatherTangent_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                            const Teuchos::ParameterList& p)
  : indexers_(indexers)
  , useTimeDerivativeSolutionVector_(false)
  , globalDataKey_("Tangent Gather Container")
{
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  indexerNames_ = p.get< Teuchos::RCP< std::vector<std::string> > >("Indexer Names");

  Teuchos::RCP<panzer::PureBasis> basis =
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis");

  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] =
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  if (p.isType<bool>("Use Time Derivative Solution Vector"))
    useTimeDerivativeSolutionVector_ = p.get<bool>("Use Time Derivative Solution Vector");

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName("Gather Tangent");
}

// **********************************************************************
template<typename EvalT,typename TRAITS,typename LO,typename GO>
void panzer::GatherTangent_BlockedEpetra<EvalT, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  indexerIds_.resize(gatherFields_.size());
  subFieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = (*indexerNames_)[fd];

    indexerIds_[fd]  = getFieldBlock(fieldName,indexers_);
    subFieldIds_[fd] = indexers_[indexerIds_[fd]]->getFieldNum(fieldName);

    TEUCHOS_ASSERT(indexerIds_[fd]>=0);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template<typename EvalT,typename TRAITS,typename LO,typename GO>
void panzer::GatherTangent_BlockedEpetra<EvalT, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  using Thyra::ProductVectorBase;

  typedef BlockedEpetraLinearObjContainer BLOC;

  // try to extract linear object container
  if (d.gedc.containsDataObject(globalDataKey_)) {
    Teuchos::RCP<const BLOC> blockedContainer = Teuchos::rcp_dynamic_cast<const BLOC>(d.gedc.getDataObject(globalDataKey_),true);

    if (useTimeDerivativeSolutionVector_)
      x_ = Teuchos::rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_dxdt());
    else
      x_ = Teuchos::rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_x());
  }
}

// **********************************************************************
template<typename EvalT,typename TRAITS,typename LO,typename GO>
void panzer::GatherTangent_BlockedEpetra<EvalT, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;

   // If x_ was not initialized, then no global evaluation data
   // container was set, in which case this evaluator becomes a no-op
   if (x_ == Teuchos::null)
     return;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // loop over the fields to be gathered
   Teuchos::ArrayRCP<const double> local_x;
   for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {

      int indexerId   = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // grab local data for inputing
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_->getNonconstVectorBlock(indexerId))->getLocalData(ptrFromRef(local_x));

      auto subRowIndexer = indexers_[indexerId];
      const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

      // gather operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         LO cellLocalId = localCellIds[worksetCellIndex];

         const std::vector<int> & LIDs = subRowIndexer->getElementLIDs(cellLocalId);

         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];

            // TEUCHOS_ASSERT(indexerId==GIDs[offset].first);
            // TEUCHOS_ASSERT(lid<local_x.size() && lid>=0);

            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = local_x[lid];
         }
      }
   }
}

// **********************************************************************

#endif
