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

#ifndef PANZER_GATHER_SOLUTION_BLOCKED_EPETRA_IMPL_HPP
#define PANZER_GATHER_SOLUTION_BLOCKED_EPETRA_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GatherSolution_Input.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Epetra_Map.h"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                             const Teuchos::ParameterList& p)
  : indexers_(indexers)
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

  std::string n = "GatherSolution (BlockedEpetra): "+firstName+" ()";
  this->setName(n);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  indexerIds_.resize(gatherFields_.size());
  subFieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = indexerNames_[fd];

    indexerIds_[fd]  = getFieldBlock(fieldName,indexers_);
    subFieldIds_[fd] = indexers_[indexerIds_[fd]]->getFieldNum(fieldName);

    TEUCHOS_ASSERT(indexerIds_[fd]>=0);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  if (has_tangent_fields_) {
    for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd)
      for (std::size_t i=0; i<tangentFields_[fd].size(); ++i)
        this->utils.setFieldData(tangentFields_[fd][i],fm);
  }

  indexerNames_.clear();  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef BlockedEpetraLinearObjContainer BLOC;

  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  RCP<GlobalEvaluationData> ged;

  // first try refactored ReadOnly container
  std::string post = useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X";
  if(d.gedc.containsDataObject(globalDataKey_+post)) {
    ged = d.gedc.getDataObject(globalDataKey_+post);

    RCP<BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<BlockedVector_ReadOnly_GlobalEvaluationData>(ged,true);

    x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(ro_ged->getGhostedVector());

    return;
  }
  else {
    ged = d.gedc.getDataObject(globalDataKey_);

    // extract linear object container
    RCP<const BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<const BlockedVector_ReadOnly_GlobalEvaluationData>(ged);
    RCP<const BlockedEpetraLinearObjContainer> blockedContainer = rcp_dynamic_cast<const BLOC>(ged);

    if(ro_ged!=Teuchos::null) {
      RCP<BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<BlockedVector_ReadOnly_GlobalEvaluationData>(ged,true);
  
      x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(ro_ged->getGhostedVector());
    }
    else if(blockedContainer!=Teuchos::null) {
      if (useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_dxdt());
      else
       x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_x());
    }
  }

  // post condition
  TEUCHOS_ASSERT(x_!=Teuchos::null); // someone has to find the x_ vector
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

   typedef BlockedEpetraLinearObjContainer BLOC;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // loop over the fields to be gathered
   Teuchos::ArrayRCP<const double> local_x;
   for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {

      PHX::MDField<ScalarT,Cell,NODE> & field = gatherFields_[fieldIndex];

      int indexerId   = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // grab local data for inputing
      Teuchos::ArrayRCP<const double> local_x;
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

            field(worksetCellIndex,basis) = local_x[lid];
         }
      }
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                             const Teuchos::ParameterList& p)
  : indexers_(indexers)
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

  std::string n = "GatherSolution Tangent (BlockedEpetra): "+firstName+" ()";
  this->setName(n);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  indexerIds_.resize(gatherFields_.size());
  subFieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = indexerNames_[fd];

    indexerIds_[fd]  = getFieldBlock(fieldName,indexers_);
    subFieldIds_[fd] = indexers_[indexerIds_[fd]]->getFieldNum(fieldName);

    TEUCHOS_ASSERT(indexerIds_[fd]>=0);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  if (has_tangent_fields_) {
    for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd)
      for (std::size_t i=0; i<tangentFields_[fd].size(); ++i)
        this->utils.setFieldData(tangentFields_[fd][i],fm);
  }

  indexerNames_.clear();  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef BlockedEpetraLinearObjContainer BLOC;

  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  RCP<GlobalEvaluationData> ged;

  // first try refactored ReadOnly container
  std::string post = useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X";
  if(d.gedc.containsDataObject(globalDataKey_+post)) {
    ged = d.gedc.getDataObject(globalDataKey_+post);

    RCP<BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<BlockedVector_ReadOnly_GlobalEvaluationData>(ged,true);

    x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(ro_ged->getGhostedVector());

    return;
  }
  else {
    ged = d.gedc.getDataObject(globalDataKey_);

    // extract linear object container
    RCP<const BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<const BlockedVector_ReadOnly_GlobalEvaluationData>(ged);
    RCP<const BlockedEpetraLinearObjContainer> blockedContainer = rcp_dynamic_cast<const BLOC>(ged);

    if(ro_ged!=Teuchos::null) {
      RCP<BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<BlockedVector_ReadOnly_GlobalEvaluationData>(ged,true);
  
      x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(ro_ged->getGhostedVector());
    }
    else if(blockedContainer!=Teuchos::null) {
      if (useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_dxdt());
      else
       x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_x());
    }
  }

  // post condition
  TEUCHOS_ASSERT(x_!=Teuchos::null); // someone has to find the x_ vector
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

   typedef BlockedEpetraLinearObjContainer BLOC;

   std::vector<std::pair<int,GO> > GIDs;
   std::vector<int> LIDs;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // loop over the fields to be gathered
   Teuchos::ArrayRCP<const double> local_x;
   for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {

      PHX::MDField<ScalarT,Cell,NODE> & field = gatherFields_[fieldIndex];

      int indexerId   = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // grab local data for inputing
      Teuchos::ArrayRCP<const double> local_x;
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

            if (!has_tangent_fields_)
              field(worksetCellIndex,basis) = local_x[lid];
            else {
              field(worksetCellIndex,basis).val() = local_x[lid];
              for (std::size_t i=0; i<tangentFields_[fieldIndex].size(); ++i)
                field(worksetCellIndex,basis).fastAccessDx(i) =
                  tangentFields_[fieldIndex][i](worksetCellIndex,basis).val();
            }
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                             const Teuchos::ParameterList& p)
  : indexers_(indexers)
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

  gatherSeedIndex_                 = input.getGatherSeedIndex();
  sensitivitiesName_               = input.getSensitivitiesName();
  disableSensitivities_            = !input.firstSensitivitiesAvailable();

  // allocate fields
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
    std::string n = "GatherSolution (BlockedEpetra, No Sensitivities): "+firstName+" ("+PHX::typeAsString<EvalT>()+")";
    this->setName(n);
  }
  else {
    std::string n = "GatherSolution (BlockedEpetra): "+firstName+" ("+PHX::typeAsString<EvalT>()+") ";
    this->setName(n);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  indexerIds_.resize(gatherFields_.size());
  subFieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = indexerNames_[fd];

    indexerIds_[fd]  = getFieldBlock(fieldName,indexers_);
    subFieldIds_[fd] = indexers_[indexerIds_[fd]]->getFieldNum(fieldName);

    TEUCHOS_ASSERT(indexerIds_[fd]>=0);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  indexerNames_.clear();  // Don't need this anymore
}

template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef BlockedEpetraLinearObjContainer BLOC;

  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

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
  if(d.gedc.containsDataObject(globalDataKey_+post)) {
    ged = d.gedc.getDataObject(globalDataKey_+post);

    RCP<BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<BlockedVector_ReadOnly_GlobalEvaluationData>(ged,true);

    x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(ro_ged->getGhostedVector());

    return;
  }
  else {
    ged = d.gedc.getDataObject(globalDataKey_);

    // extract linear object container
    RCP<const BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<const BlockedVector_ReadOnly_GlobalEvaluationData>(ged);
    RCP<const BlockedEpetraLinearObjContainer> blockedContainer = rcp_dynamic_cast<const BLOC>(ged);

    if(ro_ged!=Teuchos::null) {
      RCP<BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<BlockedVector_ReadOnly_GlobalEvaluationData>(ged,true);
  
      x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(ro_ged->getGhostedVector());
    }
    else if(blockedContainer!=Teuchos::null) {
      if (useTimeDerivativeSolutionVector_)
        x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_dxdt());
      else
       x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_x());
    }
  }

  // post condition
  TEUCHOS_ASSERT(x_!=Teuchos::null); // someone has to find the x_ vector
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   double seed_value = 0.0;
   if (useTimeDerivativeSolutionVector_) {
     seed_value = workset.alpha;
   }
   else {
     seed_value = workset.beta;
   }

   // turn off sensitivies: this may be faster if we don't expand the term
   // but I suspect not because anywhere it is used the full complement of
   // sensitivies will be needed anyway.
   if(disableSensitivities_) {
      seed_value = 0.0;
   }

   std::vector<int> blockOffsets;
   computeBlockOffsets(blockId,indexers_,blockOffsets);

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   int numDerivs = blockOffsets[blockOffsets.size()-1];

   // loop over the fields to be gathered
   for(std::size_t fieldIndex=0;
       fieldIndex<gatherFields_.size();fieldIndex++) {

      PHX::MDField<ScalarT,Cell,NODE> & field = gatherFields_[fieldIndex];

      int indexerId   = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // grab local data for inputing
      Teuchos::ArrayRCP<const double> local_x;
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_->getNonconstVectorBlock(indexerId))->getLocalData(ptrFromRef(local_x));

      auto subRowIndexer = indexers_[indexerId];
      const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

      int startBlkOffset = blockOffsets[indexerId];

      // gather operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         LO cellLocalId = localCellIds[worksetCellIndex];

         const std::vector<int> & LIDs = subRowIndexer->getElementLIDs(cellLocalId);

         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
           int offset = elmtOffset[basis];
           int lid = LIDs[offset];

           // set the value and seed the FAD object
           field(worksetCellIndex,basis) = ScalarT(numDerivs, local_x[lid]);
           field(worksetCellIndex,basis).fastAccessDx(startBlkOffset+offset) = seed_value;
         }
      }
   }
}

// **********************************************************************

#endif
