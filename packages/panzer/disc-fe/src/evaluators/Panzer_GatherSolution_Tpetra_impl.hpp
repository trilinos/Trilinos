// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GATHER_SOLUTION_TPETRA_IMPL_HPP
#define PANZER_GATHER_SOLUTION_TPETRA_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_TpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GatherSolution_Input.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"
#include "Panzer_DOFManager.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_Map.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
GatherSolution_Tpetra(
  const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
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

  std::string n = "GatherSolution (Tpetra): "+firstName+" (Residual)";
  this->setName(n);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  fieldIds_.resize(gatherFields_.size());

  const Workset & workset_0 = (*d.worksets_)[0];
  std::string blockId = this->wda(workset_0).block_id;
  scratch_offsets_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    const std::string& fieldName = indexerNames_[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    int fieldNum = fieldIds_[fd];
    const std::vector<int> & offsets = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
    scratch_offsets_[fd] = PHX::View<int*>("offsets",offsets.size());
    Kokkos::deep_copy(scratch_offsets_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(offsets.data(), offsets.size()));
  }

  scratch_lids_ = PHX::View<LO**>("lids",gatherFields_[0].extent(0),
                                                 globalIndexer_->getElementBlockGIDCount(blockId));

  indexerNames_.clear();  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   // extract linear object container
   tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject(globalDataKey_));

   if(tpetraContainer_==Teuchos::null) {
      // extract linear object container
      Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
      tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
   }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> x;
   if (useTimeDerivativeSolutionVector_)
     x = tpetraContainer_->get_dxdt();
   else
     x = tpetraContainer_->get_x();

   auto x_data = x->getLocalViewDevice(Tpetra::Access::ReadOnly);

   globalIndexer_->getElementLIDs(this->wda(workset).cell_local_ids_k,scratch_lids_);

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // gather operation for each cell in workset

   auto lids = scratch_lids_;
   for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {
     auto offsets = scratch_offsets_[fieldIndex];
     auto gather_field = gatherFields_[fieldIndex].get_static_view();

     Kokkos::parallel_for(localCellIds.size(), KOKKOS_LAMBDA (std::size_t worksetCellIndex) {
       // loop over basis functions and fill the fields
       for(std::size_t basis=0;basis<offsets.extent(0);basis++) {
         int offset = offsets(basis);
         LO lid    = lids(worksetCellIndex,offset);

         // set the value and seed the FAD object
         gather_field(worksetCellIndex,basis) = x_data(lid,0);
       }
     });
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
GatherSolution_Tpetra(
  const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
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
    // Don't allow for sharing so that we can avoid zeroing out the
    // off-diagonal values of the FAD derivative array.
    this->addUnsharedField(gatherFields_[fd].fieldTag().clone());
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
          PHX::MDField<const RealT,Cell,NODE>(tangent_field_names[fd][i],basis->functional);
        this->addDependentField(tangentFields_[fd][i]);
      }
    }
  }

  // figure out what the first active name is
  std::string firstName = "<none>";
  if(names.size()>0)
    firstName = names[0];

  std::string n = "GatherSolution (Tpetra): "+firstName+" (Tangent)";
  this->setName(n);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  fieldIds_.resize(gatherFields_.size());

  // Original implementation of tangentFields used vector of
  // vectors. The inner vectors could have different sizes for each
  // [fd]. With UVM removal, we need to use a rank 2 view of views. So
  // we need an extra vector to carry around the inner vector sizes.
  tangentInnerVectorSizes_ = PHX::View<size_t*>("tangentInnerVectorSizes_",gatherFields_.size());
  auto tangentInnerVectorSizes_host = Kokkos::create_mirror_view(tangentInnerVectorSizes_);
  size_t inner_vector_max_size = 0;
  for (std::size_t fd = 0; fd < tangentFields_.size(); ++fd) {
    inner_vector_max_size = std::max(inner_vector_max_size,tangentFields_[fd].size());
    tangentInnerVectorSizes_host(fd) = tangentFields_[fd].size();
  }
  Kokkos::deep_copy(tangentInnerVectorSizes_,tangentInnerVectorSizes_host);

  gatherFieldsVoV_.initialize("GatherSolution_Teptra<Tangent>::gatherFieldsVoV_",gatherFields_.size());
  tangentFieldsVoV_.initialize("GatherSolution_Teptra<Tangent>::tangentFieldsVoV_",gatherFields_.size(),inner_vector_max_size);

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    const std::string& fieldName = indexerNames_[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);
    gatherFieldsVoV_.addView(gatherFields_[fd].get_static_view(),fd);

    if (has_tangent_fields_) {
      for (std::size_t i=0; i<tangentFields_[fd].size(); ++i) {
        tangentFieldsVoV_.addView(tangentFields_[fd][i].get_static_view(),fd,i);
      }
    }
  }

  gatherFieldsVoV_.syncHostToDevice();
  tangentFieldsVoV_.syncHostToDevice();

  indexerNames_.clear();  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   // extract linear object container
   tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc->getDataObject(globalDataKey_));

   if(tpetraContainer_==Teuchos::null) {
      // extract linear object container
      Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
      tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
   }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;

   Teuchos::RCP<typename LOC::VectorType> x;
   if (useTimeDerivativeSolutionVector_)
     x = tpetraContainer_->get_dxdt();
   else
     x = tpetraContainer_->get_x();

   typedef typename PHX::MDField<ScalarT,Cell,NODE>::array_type::reference_type reference_type;
   auto cellLocalIdsKokkos = this->wda(workset).getLocalCellIDs();
   auto lids = globalIndexer_->getLIDs();
   auto gidFieldOffsetsVoV = Teuchos::rcp_dynamic_cast<const panzer::DOFManager>(globalIndexer_,true)->getGIDFieldOffsetsKokkos(blockId,fieldIds_);
   auto gidFieldOffsets = gidFieldOffsetsVoV.getViewDevice();
   auto gatherFieldsDevice = gatherFieldsVoV_.getViewDevice();
   auto x_view = x->getLocalViewDevice(Tpetra::Access::ReadOnly);
   auto tangentInnerVectorSizes = this->tangentInnerVectorSizes_;

   if (has_tangent_fields_) {
     auto tangentFieldsDevice = tangentFieldsVoV_.getViewDevice();
     Kokkos::parallel_for("GatherSolutionTpetra<Tangent>",cellLocalIdsKokkos.extent(0),KOKKOS_LAMBDA(const int worksetCellIndex) {
       for (size_t fieldIndex = 0; fieldIndex < gidFieldOffsets.extent(0); ++fieldIndex) { 
         for(size_t basis=0;basis<gidFieldOffsets(fieldIndex).extent(0);basis++) {
           int offset = gidFieldOffsets(fieldIndex)(basis);
           LO lid = lids(cellLocalIdsKokkos(worksetCellIndex),offset);
           auto gf_ref = (gatherFieldsDevice[fieldIndex])(worksetCellIndex,basis);
           gf_ref.val() = x_view(lid,0);
           for (std::size_t i=0; i<tangentInnerVectorSizes(fieldIndex); ++i) {
             gf_ref.fastAccessDx(i) = tangentFieldsDevice(fieldIndex,i)(worksetCellIndex,basis);
           }
         }
       }
     });
   }
   else {
     Kokkos::parallel_for("GatherSolutionTpetra<Tangent>",cellLocalIdsKokkos.extent(0),KOKKOS_LAMBDA(const int worksetCellIndex) {
       for (size_t fieldIndex = 0; fieldIndex < gidFieldOffsets.extent(0); ++fieldIndex) { 
         for(size_t basis=0;basis<gidFieldOffsets(fieldIndex).extent(0);basis++) {
           int offset = gidFieldOffsets(fieldIndex)(basis);
           LO lid = lids(cellLocalIdsKokkos(worksetCellIndex),offset);
           reference_type gf_ref = (gatherFieldsDevice[fieldIndex])(worksetCellIndex,basis);
           gf_ref.val() = x_view(lid,0);
         }
       }
     });
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
GatherSolution_Tpetra(
  const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
{
  // typedef std::vector< std::vector<std::string> > vvstring;

  GatherSolution_Input input;
  input.setParameterList(p);

  const std::vector<std::string> & names      = input.getDofNames();
  Teuchos::RCP<const panzer::PureBasis> basis = input.getBasis();
  //const vvstring & tangent_field_names        = input.getTangentNames();

  indexerNames_                    = input.getIndexerNames();
  useTimeDerivativeSolutionVector_ = input.useTimeDerivativeSolutionVector();
  globalDataKey_                   = input.getGlobalDataKey();

  gatherSeedIndex_                 = input.getGatherSeedIndex();
  sensitivitiesName_               = input.getSensitivitiesName();
  disableSensitivities_            = !input.firstSensitivitiesAvailable();

  gatherFields_.resize(names.size());
  scratch_offsets_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    PHX::MDField<ScalarT,Cell,NODE> f(names[fd],basis->functional);
    gatherFields_[fd] = f;
    this->addEvaluatedField(gatherFields_[fd]);
    // Don't allow for sharing so that we can avoid zeroing out the
    // off-diagonal values of the FAD derivative array.
    this->addUnsharedField(gatherFields_[fd].fieldTag().clone());
  }

  // figure out what the first active name is
  std::string firstName = "<none>";
  if(names.size()>0)
    firstName = names[0];

  // print out convenience
  if(disableSensitivities_) {
    std::string n = "GatherSolution (Tpetra, No Sensitivities): "+firstName+" (Jacobian)";
    this->setName(n);
  }
  else {
    std::string n = "GatherSolution (Tpetra): "+firstName+" (Jacobian) ";
    this->setName(n);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  fieldIds_.resize(gatherFields_.size());

  const Workset & workset_0 = (*d.worksets_)[0];
  std::string blockId = this->wda(workset_0).block_id;

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = indexerNames_[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    int fieldNum = fieldIds_[fd];
    const std::vector<int> & offsets = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
    scratch_offsets_[fd] = PHX::View<int*>("offsets",offsets.size());
    Kokkos::deep_copy(scratch_offsets_[fd], Kokkos::View<const int*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(offsets.data(), offsets.size()));
  }

  scratch_lids_ = PHX::View<LO**>("lids",gatherFields_[0].extent(0),
                                                 globalIndexer_->getElementBlockGIDCount(blockId));

  indexerNames_.clear();  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;
  typedef TpetraVector_ReadOnly_GlobalEvaluationData<double,LO,GO,NodeT> RO_GED;

  // manage sensitivities
  ////////////////////////////////////////////////////////////
  if(!disableSensitivities_) {
    if(d.first_sensitivities_name==sensitivitiesName_)
      applySensitivities_ = true;
    else
      applySensitivities_ = false;
  }
  else
    applySensitivities_ = false;

  ////////////////////////////////////////////////////////////

  RCP<GlobalEvaluationData> ged;

  // first try refactored ReadOnly container
  std::string post = useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X";
  if(d.gedc->containsDataObject(globalDataKey_+post)) {
    ged = d.gedc->getDataObject(globalDataKey_+post);

    RCP<RO_GED> ro_ged = rcp_dynamic_cast<RO_GED>(ged,true);

    x_vector = ro_ged->getGhostedVector_Tpetra();

    return;
  }

  ged = d.gedc->getDataObject(globalDataKey_);

  // try to extract linear object container
  {
    RCP<LOC> tpetraContainer = rcp_dynamic_cast<LOC>(ged);
    RCP<LOCPair_GlobalEvaluationData> loc_pair = rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(ged);

    if(loc_pair!=Teuchos::null) {
      Teuchos::RCP<LinearObjContainer> loc = loc_pair->getGhostedLOC();
      // extract linear object container
      tpetraContainer = rcp_dynamic_cast<LOC>(loc);
    }

    if(tpetraContainer!=Teuchos::null) {
      if (useTimeDerivativeSolutionVector_)
        x_vector = tpetraContainer->get_dxdt();
      else
        x_vector = tpetraContainer->get_x();

      return; // epetraContainer was found
    }
  }

  // try to extract an EpetraVector_ReadOnly object (this is the last resort!, it throws if not found)
  {
    RCP<RO_GED> ro_ged = rcp_dynamic_cast<RO_GED>(ged,true);

    x_vector = ro_ged->getGhostedVector_Tpetra();
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;

   double seed_value = 0.0;
   if (useTimeDerivativeSolutionVector_) {
     seed_value = workset.alpha;
   }
   else if (gatherSeedIndex_<0) {
     seed_value = workset.beta;
   }
   else if(!useTimeDerivativeSolutionVector_) {
     seed_value = workset.gather_seeds[gatherSeedIndex_];
   }
   else {
     TEUCHOS_ASSERT(false);
   }

   // turn off sensitivies: this may be faster if we don't expand the term
   // but I suspect not because anywhere it is used the full complement of
   // sensitivies will be needed anyway.
   if(!applySensitivities_)
      seed_value = 0.0;

   // Interface worksets handle DOFs from two element blocks.  The
   // derivative offset for the other element block must be shifted by
   // the derivative side of my element block.
   functor_data.dos = 0;
   if (this->wda.getDetailsIndex() == 1)
   {
     // Get the DOF count for my element block.
     functor_data.dos = globalIndexer_->getElementBlockGIDCount(workset.details(0).block_id);
   }

   // switch to a faster assembly
   bool use_seed = true;
   if(seed_value==0.0)
     use_seed = false;

   globalIndexer_->getElementLIDs(this->wda(workset).cell_local_ids_k,scratch_lids_);

   // now setup the fuctor_data, and run the parallel_for loop
   //////////////////////////////////////////////////////////////////////////////////

   functor_data.x_data = x_vector->getLocalViewDevice(Tpetra::Access::ReadOnly);
   functor_data.seed_value = seed_value;
   functor_data.lids = scratch_lids_;

   // loop over the fields to be gathered
   for(std::size_t fieldIndex=0;
       fieldIndex<gatherFields_.size();fieldIndex++) {

     // setup functor data
     functor_data.offsets = scratch_offsets_[fieldIndex];
     functor_data.field   = gatherFields_[fieldIndex];

     if(use_seed)
       Kokkos::parallel_for(workset.num_cells,*this);
     else
       Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device,NoSeed>(0,workset.num_cells),*this);
   }
   functor_data.x_data = Kokkos::View<const double**, Kokkos::LayoutLeft,PHX::Device>();
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
KOKKOS_INLINE_FUNCTION
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
operator()(const int worksetCellIndex) const
{
  // loop over basis functions and fill the fields
  for(std::size_t basis=0;basis<functor_data.offsets.extent(0);basis++) {
    int offset = functor_data.offsets(basis);
    LO lid    = functor_data.lids(worksetCellIndex,offset);

    // set the value and seed the FAD object
    if (functor_data.dos == 0)
      functor_data.field(worksetCellIndex,basis).val() = functor_data.x_data(lid,0);
    else // Interface conditions need to zero out derivative array
      functor_data.field(worksetCellIndex,basis) = ScalarT(functor_data.x_data(lid,0));

    functor_data.field(worksetCellIndex,basis).fastAccessDx(functor_data.dos + offset) = functor_data.seed_value;
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
KOKKOS_INLINE_FUNCTION
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
operator()(const NoSeed,const int worksetCellIndex) const
{
  // loop over basis functions and fill the fields
  for(std::size_t basis=0;basis<functor_data.offsets.extent(0);basis++) {
    int offset = functor_data.offsets(basis);
    LO lid    = functor_data.lids(worksetCellIndex,offset);

    // set the value and seed the FAD object
    functor_data.field(worksetCellIndex,basis).val() = functor_data.x_data(lid,0);
  }
}

// **********************************************************************

#endif
