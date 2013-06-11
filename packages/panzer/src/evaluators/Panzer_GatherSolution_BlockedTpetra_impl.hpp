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
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

#include "Tpetra_Map.hpp"

template <typename EvalT,typename Traits,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_BlockedTpetra<EvalT, Traits,S,LO,GO,NodeT>::
GatherSolution_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
{ 
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  Teuchos::RCP<panzer::PureBasis> basis = 
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis");

  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    PHX::MDField<ScalarT,Cell,NODE> field = PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(field);
  }

  this->setName("Gather Solution");
}

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_BlockedTpetra<panzer::Traits::Residual, Traits,S,LO,GO,NodeT>::
GatherSolution_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
  : gidIndexer_(indexer)
  , useTimeDerivativeSolutionVector_(false)
  , globalDataKey_("Solution Gather Container")
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

  this->setName("Gather Solution");
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Residual, Traits,S,LO,GO,NodeT>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = (*indexerNames_)[fd];
    fieldIds_[fd] = gidIndexer_->getFieldNum(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Residual, Traits,S,LO,GO,NodeT>::
preEvaluate(typename Traits::PreEvalData d)
{
   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.getDataObject(globalDataKey_),true);
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Residual, Traits,S,LO,GO,NodeT>::
evaluateFields(typename Traits::EvalData workset)
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
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

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
// Specialization: Jacobian
// **********************************************************************

template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_BlockedTpetra<panzer::Traits::Jacobian, Traits,S,LO,GO,NodeT>::
GatherSolution_BlockedTpetra(
  const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
  : gidIndexer_(indexer)
  , useTimeDerivativeSolutionVector_(false)
  , disableSensitivities_(false)
  , globalDataKey_("Solution Gather Container")
{ 
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  indexerNames_ = p.get< Teuchos::RCP< std::vector<std::string> > >("Indexer Names");

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    PHX::MDField<ScalarT,Cell,NODE> f(names[fd],dl);
    gatherFields_[fd] = f;
    this->addEvaluatedField(gatherFields_[fd]);
  }

  if (p.isType<bool>("Use Time Derivative Solution Vector"))
    useTimeDerivativeSolutionVector_ = p.get<bool>("Use Time Derivative Solution Vector");

  if (p.isType<bool>("Disable Sensitivities"))
    disableSensitivities_ = p.get<bool>("Disable Sensitivities");

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  // print out convenience
  if(disableSensitivities_)
     this->setName("Gather Solution (No Sensitivities)");
  else
     this->setName("Gather Solution");
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Jacobian, Traits,S,LO,GO,NodeT>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    //std::string fieldName = gatherFields_[fd].fieldTag().name();
    const std::string& fieldName = (*indexerNames_)[fd];
    fieldIds_[fd] = gidIndexer_->getFieldNum(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Jacobian, Traits,S,LO,GO,NodeT>::
preEvaluate(typename Traits::PreEvalData d)
{
   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.getDataObject(globalDataKey_),true);
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_BlockedTpetra<panzer::Traits::Jacobian, Traits,S,LO,GO,NodeT>::
evaluateFields(typename Traits::EvalData workset)
{ 
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

   std::vector<std::pair<int,GO> > GIDs;
   std::vector<LO> LIDs;

   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   double seed_value = 0.0;
   Teuchos::RCP<ProductVectorBase<double> > x;
   if (useTimeDerivativeSolutionVector_) {
     x = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_dxdt());
     seed_value = workset.alpha;
   }
   else {
     x = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_x()); 
     seed_value = workset.beta;
   }

   // turn off sensitivies: this may be faster if we don't expand the term
   // but I suspect not because anywhere it is used the full complement of
   // sensitivies will be needed anyway.
   if(disableSensitivities_)
      seed_value = 0.0;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

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
      for(std::size_t fieldIndex=0;
          fieldIndex<gatherFields_.size();fieldIndex++) {
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

            // set the value and seed the FAD object
            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = ScalarT(GIDs.size(), local_x[lid]);
            (gatherFields_[fieldIndex])(worksetCellIndex,basis).fastAccessDx(offset) = seed_value;
         }
      }
   }
}

// **********************************************************************

#endif
