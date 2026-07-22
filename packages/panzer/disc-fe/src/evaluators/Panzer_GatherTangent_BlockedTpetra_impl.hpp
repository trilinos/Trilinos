// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GATHER_TANGENT_BLOCKED_TPETRA_IMPL_HPP
#define PANZER_GATHER_TANGENT_BLOCKED_TPETRA_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Tpetra_Map.hpp"

template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherTangent_BlockedTpetra<EvalT, TRAITS,S,LO,GO,NodeT>::
GatherTangent_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager> & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
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
    // If blockedContainer_ is null, the evalaution is a no-op. In this
    // case we need to preserve zero initial value. Do this by not
    // sharing.
    this->addUnsharedField(gatherFields_[fd].fieldTag().clone());
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
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  const Workset & workset_0 = (*d.worksets_)[0];
  const std::string blockId = this->wda(workset_0).block_id;

  fieldIds_.resize(gatherFields_.size());
  fieldOffsets_.resize(gatherFields_.size());
  fieldGlobalIndexers_.resize(gatherFields_.size());
  productVectorBlockIndex_.resize(gatherFields_.size());
  int maxElementBlockGIDCount = -1;
  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = (*indexerNames_)[fd];
    const int globalFieldNum = globalIndexer_->getFieldNum(fieldName); // Field number in the aggregate BlockDOFManager
    productVectorBlockIndex_[fd] = globalIndexer_->getFieldBlock(globalFieldNum);
    fieldGlobalIndexers_[fd] = globalIndexer_->getFieldDOFManagers()[productVectorBlockIndex_[fd]];
    fieldIds_[fd] = fieldGlobalIndexers_[fd]->getFieldNum(fieldName); // Field number in the sub-global-indexer

    const std::vector<int>& offsets = fieldGlobalIndexers_[fd]->getGIDFieldOffsets(blockId,fieldIds_[fd]);
    fieldOffsets_[fd] = PHX::View<int*>("GatherSolution_BlockedTpetra(Residual):fieldOffsets",offsets.size());
    auto hostFieldOffsets = Kokkos::create_mirror_view(fieldOffsets_[fd]);
    for(std::size_t i=0; i < offsets.size(); ++i)
      hostFieldOffsets(i) = offsets[i];
    Kokkos::deep_copy(fieldOffsets_[fd],hostFieldOffsets);

    maxElementBlockGIDCount = std::max(fieldGlobalIndexers_[fd]->getElementBlockGIDCount(blockId),maxElementBlockGIDCount);
  }

  // We will use one workset lid view for all fields, but has to be
  // sized big enough to hold the largest elementBlockGIDCount in the
  // ProductVector.
  worksetLIDs_ = PHX::View<LO**>("GatherSolution_BlockedTpetra(Residual):worksetLIDs",
                                                gatherFields_[0].extent(0),
                                                maxElementBlockGIDCount);

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherTangent_BlockedTpetra<EvalT, TRAITS,S,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // try to extract linear object container
  if (d.gedc->containsDataObject(globalDataKey_)) {
    Teuchos::RCP<GlobalEvaluationData> ged = d.gedc->getDataObject(globalDataKey_);
    Teuchos::RCP<LOCPair_GlobalEvaluationData> loc_pair =
      Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(ged);

    if(loc_pair!=Teuchos::null) {
      Teuchos::RCP<LinearObjContainer> loc = loc_pair->getGhostedLOC();
      blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(loc,true);
    }

    if(blockedContainer_==Teuchos::null) {
      blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(ged,true);
    }
  }
}

// **********************************************************************
template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherTangent_BlockedTpetra<EvalT, TRAITS,S,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorBase;
  using Thyra::ProductVectorBase;

  // If blockedContainer_ was not initialized, then no global evaluation data
  // container was set, in which case this evaluator becomes a no-op
  if (blockedContainer_ == Teuchos::null) return;
  
  const PHX::View<const int*>& localCellIds = this->wda(workset).cell_local_ids_k;
  
  RCP<ProductVectorBase<ScalarT>> thyraBlockSolution;
  if (useTimeDerivativeSolutionVector_)
    thyraBlockSolution = rcp_dynamic_cast<ProductVectorBase<ScalarT>>(blockedContainer_->get_dxdt(),true);
  else
    thyraBlockSolution = rcp_dynamic_cast<ProductVectorBase<ScalarT>>(blockedContainer_->get_x(),true);
  
  // Loop over gathered fields
  int currentWorksetLIDSubBlock = -1;
  for (std::size_t fieldIndex = 0; fieldIndex < gatherFields_.size(); fieldIndex++) {
    // workset LIDs only change for different sub blocks 
    if (productVectorBlockIndex_[fieldIndex] != currentWorksetLIDSubBlock) {
      const std::string blockId = this->wda(workset).block_id;
      const int num_dofs = fieldGlobalIndexers_[fieldIndex]->getElementBlockGIDCount(blockId);
      fieldGlobalIndexers_[fieldIndex]->getElementLIDs(localCellIds,worksetLIDs_,num_dofs); 
      currentWorksetLIDSubBlock = productVectorBlockIndex_[fieldIndex];
    }

    const auto& tpetraSolution = *((rcp_dynamic_cast<Thyra::TpetraVector<ScalarT,LO,GO,NodeT>>(thyraBlockSolution->getNonconstVectorBlock(productVectorBlockIndex_[fieldIndex]),true))->getTpetraVector());
    const auto& kokkosSolution = tpetraSolution.getLocalViewDevice(Tpetra::Access::ReadOnly);

    // Class data fields for lambda capture
    const auto& fieldOffsets = fieldOffsets_[fieldIndex];
    const auto& worksetLIDs = worksetLIDs_;
    const auto& fieldValues = gatherFields_[fieldIndex];

    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device>(0,workset.num_cells), KOKKOS_LAMBDA (const int& cell) {       
      for(int basis=0; basis < static_cast<int>(fieldOffsets.size()); ++basis) {
        const int lid = worksetLIDs(cell,fieldOffsets(basis));
        fieldValues(cell,basis) = kokkosSolution(lid,0);        
      }
    });
  }

}

// **********************************************************************

#endif
