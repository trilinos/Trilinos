// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GATHER_TANGENT_TPETRA_IMPL_HPP
#define PANZER_GATHER_TANGENT_TPETRA_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"
#include "Panzer_DOFManager.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_Map.hpp"

template<typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::GatherTangent_Tpetra<EvalT, TRAITS,LO,GO,NodeT>::
GatherTangent_Tpetra(
  const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
  , useTimeDerivativeSolutionVector_(false)
  , globalDataKey_("Tangent Gather Container")
{
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  indexerNames_ = p.get< Teuchos::RCP< std::vector<std::string> > >("Indexer Names");

  // this is beging to fix the issues with incorrect use of const
  Teuchos::RCP<const panzer::PureBasis> basis;
  if(p.isType< Teuchos::RCP<panzer::PureBasis> >("Basis"))
    basis = p.get< Teuchos::RCP<panzer::PureBasis> >("Basis");
  else
    basis = p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis");

  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] =
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
    // If tpetraContainer_ is null, the evalaution is a no-op. In this
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
template<typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherTangent_Tpetra<EvalT, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  fieldIds_.resize(gatherFields_.size());

  gatherFieldsVoV_.initialize("GatherSolution_Teptra<Tangent>",gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    const std::string& fieldName = (*indexerNames_)[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);
    gatherFieldsVoV_.addView(gatherFields_[fd].get_static_view(),fd);
  }

  gatherFieldsVoV_.syncHostToDevice();

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template<typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherTangent_Tpetra<EvalT, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

  // try to extract linear object container
  if (d.gedc->containsDataObject(globalDataKey_)) {
    RCP<GlobalEvaluationData> ged = d.gedc->getDataObject(globalDataKey_);
    RCP<LOCPair_GlobalEvaluationData> loc_pair =
      rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(ged);

    if(loc_pair!=Teuchos::null) {
      Teuchos::RCP<LinearObjContainer> loc = loc_pair->getGhostedLOC();
      tpetraContainer_ = rcp_dynamic_cast<LOC>(loc,true);
    }

    if(tpetraContainer_==Teuchos::null) {
      tpetraContainer_ = rcp_dynamic_cast<LOC>(ged,true);
    }
  }
}

// **********************************************************************
template<typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherTangent_Tpetra<EvalT, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
  // If tpetraContainer_ was not initialized, then no global evaluation data
  // container was set, in which case this evaluator becomes a no-op
  if (tpetraContainer_ == Teuchos::null)
    return;

  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;
  // for convenience pull out some objects from workset
  std::string blockId = this->wda(workset).block_id;

  Teuchos::RCP<typename LOC::VectorType> x;
  if (useTimeDerivativeSolutionVector_)
    x = tpetraContainer_->get_dxdt();
  else
    x = tpetraContainer_->get_x();

  auto cellLocalIdsKokkos = this->wda(workset).getLocalCellIDs();
  auto lids = globalIndexer_->getLIDs();
  auto vov = Teuchos::rcp_dynamic_cast<const panzer::DOFManager>(globalIndexer_,true)->getGIDFieldOffsetsKokkos(blockId,fieldIds_);
  auto gidFieldOffsets = vov.getViewDevice();
  auto gatherFieldsDevice = gatherFieldsVoV_.getViewDevice();
  auto x_view = x->getLocalViewDevice(Tpetra::Access::ReadWrite);
  Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{cellLocalIdsKokkos.extent(0),gidFieldOffsets.extent(0)});
  Kokkos::parallel_for("GatherSolutionTpetra<Tangent>",policy,KOKKOS_LAMBDA(const int worksetCellIndex, const int fieldIndex) {
    for(std::size_t basis=0;basis<gidFieldOffsets(fieldIndex).extent(0);basis++) {
      int offset = gidFieldOffsets(fieldIndex)(basis);
      LO lid = lids(cellLocalIdsKokkos(worksetCellIndex),offset);
      auto& gf_ref = (gatherFieldsDevice[fieldIndex])(worksetCellIndex,basis);
      gf_ref = x_view(lid,0);
    }
  });
}

// **********************************************************************

#endif
