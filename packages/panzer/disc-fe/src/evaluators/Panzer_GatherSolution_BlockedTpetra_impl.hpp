// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GATHER_SOLUTION_BLOCKED_EPETRA_IMPL_HPP
#define PANZER_GATHER_SOLUTION_BLOCKED_EPETRA_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_GatherSolution_Input.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Tpetra_Map.hpp"

template <typename EvalT,typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_BlockedTpetra<EvalT, TRAITS,S,LO,GO,NodeT>::
GatherSolution_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager> & indexer,
  const Teuchos::ParameterList& p)
{
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  Teuchos::RCP<panzer::PureBasis> basis =
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis");

  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    PHX::MDField<ScalarT,Cell,NODE> field = PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(field.fieldTag());
  }

  this->setName("Gather Solution");
}

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_BlockedTpetra<panzer::Traits::Residual, TRAITS,S,LO,GO,NodeT>::
GatherSolution_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager> & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
  , has_tangent_fields_(false)
{
  typedef std::vector< std::vector<std::string> > vvstring;

  GatherSolution_Input input;
  input.setParameterList(p);

  const std::vector<std::string> & names      = input.getDofNames();
  Teuchos::RCP<const panzer::PureBasis> basis = input.getBasis();
  const vvstring & tangent_field_names        = input.getTangentNames();

  indexerNames_                    = input.getIndexerNames();
  useTimeDerivativeSolutionVector_ = input.useTimeDerivativeSolutionVector();
  globalDataKey_                   = input.getGlobalDataKey();

  // allocate fields
  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] =
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  // Setup dependent tangent fields if requested
  if (tangent_field_names.size()>0) {
    TEUCHOS_ASSERT(gatherFields_.size() == tangent_field_names.size());

    has_tangent_fields_ = true;
    tangentFields_.resize(gatherFields_.size());
    for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
      tangentFields_[fd].resize(tangent_field_names[fd].size());
      for (std::size_t i=0; i<tangent_field_names[fd].size(); ++i) {
        tangentFields_[fd][i] =
          PHX::MDField<const ScalarT,Cell,NODE>(tangent_field_names[fd][i],basis->functional);
        this->addDependentField(tangentFields_[fd][i]);
      }
    }
  }

  // figure out what the first active name is
  std::string firstName = "<none>";
  if(names.size()>0)
    firstName = names[0];

  std::string n = "GatherSolution (BlockedTpetra): "+firstName+" (Residual)";
  this->setName(n);
}

// **********************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Residual, TRAITS,S,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  const Workset & workset_0 = (*d.worksets_)[0];
  const std::string blockId = this->wda(workset_0).block_id;

  fieldIds_.resize(gatherFields_.size());
  fieldOffsets_.resize(gatherFields_.size());
  fieldGlobalIndexers_.resize(gatherFields_.size());
  productVectorBlockIndex_.resize(gatherFields_.size());
  int maxElementBlockGIDCount = -1;
  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = indexerNames_[fd];
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

  indexerNames_.clear();  // Don't need this anymore
}

// **********************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Residual, TRAITS,S,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_),true);
}

// **********************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Residual, TRAITS,S,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorBase;
  using Thyra::ProductVectorBase;
  
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
// Specialization: Tangent
// **********************************************************************

template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_BlockedTpetra<panzer::Traits::Tangent, TRAITS,S,LO,GO,NodeT>::
GatherSolution_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager> & indexer,
  const Teuchos::ParameterList& p)
  : gidIndexer_(indexer)
  , has_tangent_fields_(false)
{
  typedef std::vector< std::vector<std::string> > vvstring;

  GatherSolution_Input input;
  input.setParameterList(p);

  const std::vector<std::string> & names      = input.getDofNames();
  Teuchos::RCP<const panzer::PureBasis> basis = input.getBasis();
  const vvstring & tangent_field_names        = input.getTangentNames();

  indexerNames_                    = input.getIndexerNames();
  useTimeDerivativeSolutionVector_ = input.useTimeDerivativeSolutionVector();
  globalDataKey_                   = input.getGlobalDataKey();

  // allocate fields
  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] =
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  // Setup dependent tangent fields if requested
  if (tangent_field_names.size()>0) {
    TEUCHOS_ASSERT(gatherFields_.size() == tangent_field_names.size());

    has_tangent_fields_ = true;
    tangentFields_.resize(gatherFields_.size());
    for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
      tangentFields_[fd].resize(tangent_field_names[fd].size());
      for (std::size_t i=0; i<tangent_field_names[fd].size(); ++i) {
        tangentFields_[fd][i] =
          PHX::MDField<const ScalarT,Cell,NODE>(tangent_field_names[fd][i],basis->functional);
        this->addDependentField(tangentFields_[fd][i]);
      }
    }
  }

  // figure out what the first active name is
  std::string firstName = "<none>";
  if(names.size()>0)
    firstName = names[0];

  std::string n = "GatherSolution (BlockedTpetra): "+firstName+" (Tangent)";
  this->setName(n);
}

// **********************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Tangent, TRAITS,S,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = indexerNames_[fd];
    fieldIds_[fd] = gidIndexer_->getFieldNum(fieldName);
  }

  indexerNames_.clear();  // Don't need this anymore
}

// **********************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Tangent, TRAITS,S,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_),true);
}

// **********************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Tangent, TRAITS,S,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

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

      gidIndexer_->getElementGIDsPair(cellLocalId,GIDs,blockId);

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

            if (!has_tangent_fields_)
              (gatherFields_[fieldIndex])(worksetCellIndex,basis) = local_x[lid];
            else {
              (gatherFields_[fieldIndex])(worksetCellIndex,basis).val() = local_x[lid];
              for (std::size_t i=0; i<tangentFields_[fieldIndex].size(); ++i)
                (gatherFields_[fieldIndex])(worksetCellIndex,basis).fastAccessDx(i) =
                  tangentFields_[fieldIndex][i](worksetCellIndex,basis).val();
            }
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,S,LO,GO,NodeT>::
GatherSolution_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager> & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
{
  GatherSolution_Input input;
  input.setParameterList(p);

  const std::vector<std::string> & names      = input.getDofNames();
  Teuchos::RCP<const panzer::PureBasis> basis = input.getBasis();

  indexerNames_                    = input.getIndexerNames();
  useTimeDerivativeSolutionVector_ = input.useTimeDerivativeSolutionVector();
  globalDataKey_                   = input.getGlobalDataKey();

  disableSensitivities_            = !input.firstSensitivitiesAvailable();

  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    PHX::MDField<ScalarT,Cell,NODE> f(names[fd],basis->functional);
    gatherFields_[fd] = f;
    this->addEvaluatedField(gatherFields_[fd]);
  }

  // figure out what the first active name is
  std::string firstName = "<none>";
  if(names.size()>0)
    firstName = names[0];

  // print out convenience
  if(disableSensitivities_) {
    std::string n = "GatherSolution (BlockedTpetra, No Sensitivities): "+firstName+" (Jacobian)";
    this->setName(n);
  }
  else {
    std::string n = "GatherSolution (BlockedTpetra): "+firstName+" ("+PHX::print<EvalT>()+") ";
    this->setName(n);
  }
}

// **********************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,S,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  const Workset & workset_0 = (*d.worksets_)[0];
  const std::string blockId = this->wda(workset_0).block_id;

  fieldIds_.resize(gatherFields_.size());
  fieldOffsets_.resize(gatherFields_.size());
  productVectorBlockIndex_.resize(gatherFields_.size());
  int maxElementBlockGIDCount = -1;
  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {

    const std::string fieldName = indexerNames_[fd];
    const int globalFieldNum = globalIndexer_->getFieldNum(fieldName); // Field number in the aggregate BlockDOFManager
    productVectorBlockIndex_[fd] = globalIndexer_->getFieldBlock(globalFieldNum);
    const auto& subGlobalIndexer = globalIndexer_->getFieldDOFManagers()[productVectorBlockIndex_[fd]];
    fieldIds_[fd] = subGlobalIndexer->getFieldNum(fieldName); // Field number in the sub-global-indexer

    const std::vector<int>& offsets = subGlobalIndexer->getGIDFieldOffsets(blockId,fieldIds_[fd]);
    fieldOffsets_[fd] = PHX::View<int*>("GatherSolution_BlockedTpetra(Jacobian):fieldOffsets",offsets.size());
    auto hostOffsets = Kokkos::create_mirror_view(fieldOffsets_[fd]);
    for (std::size_t i=0; i < offsets.size(); ++i)
      hostOffsets(i) = offsets[i];
    Kokkos::deep_copy(fieldOffsets_[fd], hostOffsets);
    maxElementBlockGIDCount = std::max(subGlobalIndexer->getElementBlockGIDCount(blockId),maxElementBlockGIDCount);
  }

  // We will use one workset lid view for all fields, but has to be
  // sized big enough to hold the largest elementBlockGIDCount in the
  // ProductVector.
  worksetLIDs_ = PHX::View<LO**>("ScatterResidual_BlockedTpetra(Residual):worksetLIDs",
                                                gatherFields_[0].extent(0),
                                                maxElementBlockGIDCount);

  // Compute the block offsets
  const auto& blockGlobalIndexers = globalIndexer_->getFieldDOFManagers();
  const int numBlocks = static_cast<int>(globalIndexer_->getFieldDOFManagers().size());
  blockOffsets_ = PHX::View<LO*>("GatherSolution_BlockedTpetra(Jacobian):blockOffsets_",
                                                numBlocks+1); // Number of blocks, plus a sentinel
  const auto hostBlockOffsets = Kokkos::create_mirror_view(blockOffsets_);
  for (int blk=0;blk<numBlocks;++blk) {
    int blockOffset = globalIndexer_->getBlockGIDOffset(blockId,blk);
    hostBlockOffsets(blk) = blockOffset;
  }
  hostBlockOffsets(numBlocks) = hostBlockOffsets(numBlocks-1) + blockGlobalIndexers[blockGlobalIndexers.size()-1]->getElementBlockGIDCount(blockId);
  Kokkos::deep_copy(blockOffsets_,hostBlockOffsets);

  indexerNames_.clear();  // Don't need this anymore
}

template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,S,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_),true);
}

// **********************************************************************
template <typename TRAITS,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,S,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorBase;
  using Thyra::ProductVectorBase;

  const auto& localCellIds = this->wda(workset).cell_local_ids_k;
  
  RealType seedValue = RealType(0.0);
  RCP<ProductVectorBase<double>> blockedSolution;
  if (useTimeDerivativeSolutionVector_) {
    blockedSolution = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_dxdt());
    seedValue = workset.alpha;
  }
  else {
    blockedSolution = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_x());
    seedValue = workset.beta;
  }

  // turn off sensitivies: this may be faster if we don't expand the term
  // but I suspect not because anywhere it is used the full complement of
  // sensitivies will be needed anyway.
  if(disableSensitivities_)
    seedValue = 0.0;

  // Loop over fields to gather
  int currentWorksetLIDSubBlock = -1;
  for (std::size_t fieldIndex = 0; fieldIndex < gatherFields_.size(); fieldIndex++) {
    // workset LIDs only change if in different sub blocks 
    if (productVectorBlockIndex_[fieldIndex] != currentWorksetLIDSubBlock) {
      const auto& blockIndexer = globalIndexer_->getFieldDOFManagers()[productVectorBlockIndex_[fieldIndex]];
      const std::string blockId = this->wda(workset).block_id;
      const int num_dofs = globalIndexer_->getFieldDOFManagers()[productVectorBlockIndex_[fieldIndex]]->getElementBlockGIDCount(blockId);
      blockIndexer->getElementLIDs(localCellIds,worksetLIDs_,num_dofs); 
      currentWorksetLIDSubBlock = productVectorBlockIndex_[fieldIndex];
    }

    const int blockRowIndex = productVectorBlockIndex_[fieldIndex];
    const auto& subblockSolution = *((rcp_dynamic_cast<Thyra::TpetraVector<RealType,LO,GO,NodeT>>(blockedSolution->getNonconstVectorBlock(blockRowIndex),true))->getTpetraVector());
    const auto kokkosSolution = subblockSolution.getLocalViewDevice(Tpetra::Access::ReadOnly);

    // Class data fields for lambda capture
    const PHX::View<const int*> fieldOffsets = fieldOffsets_[fieldIndex];
    const PHX::View<const LO**> worksetLIDs = worksetLIDs_;
    const PHX::View<ScalarT**> fieldValues = gatherFields_[fieldIndex].get_static_view();        
    const PHX::View<const LO*> blockOffsets = blockOffsets_;
    auto blockOffsets_h = Kokkos::create_mirror_view(blockOffsets);
    Kokkos::deep_copy(blockOffsets_h, blockOffsets);
    const int blockStart = blockOffsets_h(blockRowIndex);

    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device>(0,workset.num_cells), KOKKOS_LAMBDA (const int& cell) {  
      for (int basis=0; basis < static_cast<int>(fieldOffsets.size()); ++basis) {
        const int rowLID = worksetLIDs(cell,fieldOffsets(basis));
	fieldValues(cell,basis).zero();
        fieldValues(cell,basis).val() = kokkosSolution(rowLID,0);
        fieldValues(cell,basis).fastAccessDx(blockStart+fieldOffsets(basis)) = seedValue;
      }
    });

  }
}

// **********************************************************************

#endif
