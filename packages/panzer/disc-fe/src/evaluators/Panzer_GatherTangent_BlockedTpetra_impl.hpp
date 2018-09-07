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

#ifndef PANZER_GATHER_TANGENT_BLOCKED_TPETRA_IMPL_HPP
#define PANZER_GATHER_TANGENT_BLOCKED_TPETRA_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Tpetra_Map.hpp"

template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherTangent_BlockedTpetra<EvalT, TRAITS,S,LO,GO,NodeT>::
GatherTangent_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
  : gidIndexer_(indexer)
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
template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherTangent_BlockedTpetra<EvalT, TRAITS,S,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = (*indexerNames_)[fd];
    fieldIds_[fd] = gidIndexer_->getFieldNum(fieldName);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherTangent_BlockedTpetra<EvalT, TRAITS,S,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // try to extract linear object container
  if (d.gedc->containsDataObject(globalDataKey_)) {
    blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_),true);
  }
}

// **********************************************************************
template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherTangent_BlockedTpetra<EvalT, TRAITS,S,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

   // If blockedContainer_ was not initialized, then no global evaluation data
   // container was set, in which case this evaluator becomes a no-op
   if (blockedContainer_ == Teuchos::null)
     return;

   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setShowProcRank(true);
   out.setOutputToRootOnly(-1);

   std::vector<std::pair<int,GO> > GIDs;
   std::vector<LO> LIDs;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<ProductVectorBase<double> > x;
   if (useTimeDerivativeSolutionVector_)
     x = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_dxdt());
   else
     x = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_x());

   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      LO cellLocalId = localCellIds[worksetCellIndex];

      gidIndexer_->getElementGIDs(cellLocalId,GIDs,blockId);

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++) {
         // used for doing local ID lookups
         RCP<const MapType> x_map = blockedContainer_->getMapForBlock(GIDs[i].first);

         LIDs[i] = x_map->getLocalElement(GIDs[i].second);
      }

      // loop over the fields to be gathered
      Teuchos::ArrayRCP<const double> local_x;
      for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         int indexerId = gidIndexer_->getFieldBlock(fieldNum);

         // grab local data for inputing
         RCP<SpmdVectorBase<double> > block_x = rcp_dynamic_cast<SpmdVectorBase<double> >(x->getNonconstVectorBlock(indexerId));
         block_x->getLocalData(ptrFromRef(local_x));

         const std::vector<int> & elmtOffset = gidIndexer_->getGIDFieldOffsets(blockId,fieldNum);

         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];

            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = local_x[lid];
         }
      }
   }
}

// **********************************************************************

#endif
