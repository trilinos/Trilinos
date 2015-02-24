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

#ifndef PANZER_GATHER_SOLUTION_TPETRA_IMPL_HPP
#define PANZER_GATHER_SOLUTION_TPETRA_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_Map.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
GatherSolution_Tpetra(
  const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
  , useTimeDerivativeSolutionVector_(false)
  , globalDataKey_("Solution Gather Container")
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
  }

  if (p.isType<bool>("Use Time Derivative Solution Vector"))
    useTimeDerivativeSolutionVector_ = p.get<bool>("Use Time Derivative Solution Vector");

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName("Gather Solution");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d, 
		      PHX::FieldManager<TRAITS>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    const std::string& fieldName = (*indexerNames_)[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   // extract linear object container
   tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc.getDataObject(globalDataKey_));

   if(tpetraContainer_==Teuchos::null) {
      // extract linear object container
      Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc.getDataObject(globalDataKey_),true)->getGhostedLOC();
      tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
   }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   std::vector<LO> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> x;
   if (useTimeDerivativeSolutionVector_)
     x = tpetraContainer_->get_dxdt();
   else
     x = tpetraContainer_->get_x(); 

   Teuchos::ArrayRCP<const double> x_array = x->get1dView();
 
   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!
 
   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];
 
      LIDs = globalIndexer_->getElementLIDs(cellLocalId); 
 
      // loop over the fields to be gathered
      for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
 
         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            LO lid = LIDs[offset];
            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = x_array[lid];
         }
      }
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
GatherSolution_Tpetra(
  const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
  , useTimeDerivativeSolutionVector_(false)
  , globalDataKey_("Solution Gather Container")
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
  }

  if (p.isType<bool>("Use Time Derivative Solution Vector"))
    useTimeDerivativeSolutionVector_ = p.get<bool>("Use Time Derivative Solution Vector");

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName("Gather Solution");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d, 
		      PHX::FieldManager<TRAITS>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    const std::string& fieldName = (*indexerNames_)[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   // extract linear object container
   tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc.getDataObject(globalDataKey_));

   if(tpetraContainer_==Teuchos::null) {
      // extract linear object container
      Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc.getDataObject(globalDataKey_),true)->getGhostedLOC();
      tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
   }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Tangent, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   std::vector<LO> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<typename LOC::VectorType> x;
   if (useTimeDerivativeSolutionVector_)
     x = tpetraContainer_->get_dxdt();
   else
     x = tpetraContainer_->get_x(); 

   Teuchos::ArrayRCP<const double> x_array = x->get1dView();
 
   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!
 
   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];
 
      LIDs = globalIndexer_->getElementLIDs(cellLocalId); 
 
      // loop over the fields to be gathered
      for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
 
         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            LO lid = LIDs[offset];
            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = x_array[lid];
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
GatherSolution_Tpetra(
  const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
  , useTimeDerivativeSolutionVector_(false)
  , globalDataKey_("Solution Gather Container")
  , gatherSeedIndex_(-1)
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
    PHX::MDField<ScalarT,Cell,NODE> f(names[fd],basis->functional);
    gatherFields_[fd] = f;
    this->addEvaluatedField(gatherFields_[fd]);
  }

  if (p.isType<bool>("Use Time Derivative Solution Vector"))
    useTimeDerivativeSolutionVector_ = p.get<bool>("Use Time Derivative Solution Vector");

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  if (p.isType<int>("Gather Seed Index")) {
     gatherSeedIndex_ = p.get<int>("Gather Seed Index");
  }

  this->setName("Gather Solution");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d, 
		      PHX::FieldManager<TRAITS>& fm)
{
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_->size());

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = (*indexerNames_)[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   // extract linear object container
   tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(d.gedc.getDataObject(globalDataKey_));

   if(tpetraContainer_==Teuchos::null) {
      // extract linear object container
      Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc.getDataObject(globalDataKey_),true)->getGhostedLOC();
      tpetraContainer_ = Teuchos::rcp_dynamic_cast<LOC>(loc);
   }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   typedef TpetraLinearObjContainer<double,LO,GO,NodeT> LOC;

   std::vector<LO> LIDs;

   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;

   Teuchos::RCP<typename LOC::VectorType> x;
   double seed_value = 0.0;
   if (useTimeDerivativeSolutionVector_) {
     x = tpetraContainer_->get_dxdt();
     seed_value = workset.alpha;
   }
   else if (gatherSeedIndex_<0) {
     x = tpetraContainer_->get_x();
     seed_value = workset.beta;
   }
   else if(!useTimeDerivativeSolutionVector_) {
     x = tpetraContainer_->get_x();
     seed_value = workset.gather_seeds[gatherSeedIndex_];
   }
   else {
     TEUCHOS_ASSERT(false);
   }

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

/*
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
   Teuchos::ArrayRCP<const double> x_array = x->get1dView();

   // loop over the fields to be gathered
   for(std::size_t fieldIndex=0;
       fieldIndex<gatherFields_.size();fieldIndex++) {
     int fieldNum = fieldIds_[fieldIndex];
     const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

     PHX::MDField<ScalarT,Cell,NODE> field = gatherFields_[fieldIndex];

     // gather operation for each cell in workset
     for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
       std::size_t cellLocalId = localCellIds[worksetCellIndex];

       LIDs = globalIndexer_->getElementLIDs(cellLocalId); 

       // loop over basis functions and fill the fields
       for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
         int offset = elmtOffset[basis];
         LO lid = LIDs[offset];

         // set the value and seed the FAD object
         field(worksetCellIndex,basis) = ScalarT(LIDs.size(), x_array[lid]);
         field(worksetCellIndex,basis).fastAccessDx(offset) = seed_value;
       }
     }
   }
*/
   fctr_localCellIds = workset.cell_local_ids;
   fctr_x_array = x->get1dView();
   fctr_seed_value = seed_value;

   // loop over the fields to be gathered
   for(std::size_t fieldIndex=0;
       fieldIndex<gatherFields_.size();fieldIndex++) {
     int fieldNum = fieldIds_[fieldIndex];

     // setup functor data
     fctr_elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
     fctr_field = gatherFields_[fieldIndex];

     Kokkos::parallel_for(workset.num_cells,*this);
   }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO,typename NodeT>
KOKKOS_INLINE_FUNCTION
void panzer::GatherSolution_Tpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
operator()(const int worksetCellIndex) const
{
  std::vector<LO> LIDs;

  // gather operation for each cell in workset
  std::size_t cellLocalId = fctr_localCellIds[worksetCellIndex];

  LIDs = globalIndexer_->getElementLIDs(cellLocalId); 

  // loop over basis functions and fill the fields
  for(std::size_t basis=0;basis<fctr_elmtOffset.size();basis++) {
    int offset = fctr_elmtOffset[basis];
    LO lid = LIDs[offset];

    // set the value and seed the FAD object
    fctr_field(worksetCellIndex,basis) = ScalarT(LIDs.size(), fctr_x_array[lid]);
    fctr_field(worksetCellIndex,basis).fastAccessDx(offset) = fctr_seed_value;
  }
}

// **********************************************************************

#endif
