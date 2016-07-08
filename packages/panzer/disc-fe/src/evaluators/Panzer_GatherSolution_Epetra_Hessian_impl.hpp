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

#ifndef __Panzer_GatherSolution_Epetra_Hessian_impl_hpp__
#define __Panzer_GatherSolution_Epetra_Hessian_impl_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_ParameterList_GlobalEvaluationData.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace panzer {

// **********************************************************************
// Specialization: Hessian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::GatherSolution_Epetra<panzer::Traits::Hessian, TRAITS,LO,GO>::
GatherSolution_Epetra(
  const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
  : globalIndexer_(indexer)
{
  typedef std::vector< std::vector<std::string> > vvstring;

  GatherSolution_Input input;
  input.setParameterList(p);

  const std::vector<std::string> & names      = input.getDofNames();
  Teuchos::RCP<const panzer::PureBasis> basis = input.getBasis();

  indexerNames_                    = input.getIndexerNames();
  useTimeDerivativeSolutionVector_ = input.useTimeDerivativeSolutionVector();
  globalDataKey_                   = input.getGlobalDataKey();

  gatherSeedIndex_                 = input.getGatherSeedIndex();
  sensitivitiesName_               = input.getSensitivitiesName();
  firstSensitivitiesAvailable_     = input.firstSensitivitiesAvailable();

  secondSensitivitiesAvailable_    = input.secondSensitivitiesAvailable();
  sensitivities2ndPrefix_          = input.getSecondSensitivityDataKeyPrefix();

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
  if(!firstSensitivitiesAvailable_) {
    std::string n = "GatherSolution (Epetra, No First Sensitivities): "+firstName+" (Hessian)";
    this->setName(n);
  }
  else {
    std::string n = "GatherSolution (Epetra): "+firstName+" ("+PHX::typeAsString<EvalT>()+") ";
    this->setName(n);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_Epetra<panzer::Traits::Hessian, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& fm)
{
  // globalIndexer_ = d.globalIndexer_;
  TEUCHOS_ASSERT(gatherFields_.size() == indexerNames_.size());

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = indexerNames_[fd];
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // this is the error return code, raise the alarm
    if(fieldIds_[fd]==-1) {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "GatherSolution_Epetra<Hessian>: Could not find field \"" + fieldName + "\" in the global indexer. ");
         // wouldn't it be nice to print more information???
    }

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }

  indexerNames_.clear();  // Don't need this anymore
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_Epetra<panzer::Traits::Hessian, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // manage sensitivities
  ////////////////////////////////////////////////////////////
  if(firstSensitivitiesAvailable_) {
    if(d.first_sensitivities_name==sensitivitiesName_)
      firstApplySensitivities_ = true;
    else
      firstApplySensitivities_ = false;
  }
  else
    firstApplySensitivities_ = false;

  if(secondSensitivitiesAvailable_) {
    if(d.second_sensitivities_name==sensitivitiesName_)
      secondApplySensitivities_ = true;
    else
      secondApplySensitivities_ = false;
  }
  else
    secondApplySensitivities_ = false;

  ////////////////////////////////////////////////////////////

  RCP<GlobalEvaluationData> ged;

  // first try refactored ReadOnly container
  std::string post = useTimeDerivativeSolutionVector_ ? " - Xdot" : " - X";
  if(d.gedc.containsDataObject(globalDataKey_+post)) {
    ged = d.gedc.getDataObject(globalDataKey_+post);

    RCP<EpetraVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<EpetraVector_ReadOnly_GlobalEvaluationData>(ged,true);

    x_ = ro_ged->getGhostedVector_Epetra();
  }
  else if(d.gedc.containsDataObject(globalDataKey_)) {

    ged = d.gedc.getDataObject(globalDataKey_);
  
    // try to extract linear object container
    RCP<EpetraVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<EpetraVector_ReadOnly_GlobalEvaluationData>(ged);
    RCP<EpetraLinearObjContainer> epetraContainer = rcp_dynamic_cast<EpetraLinearObjContainer>(ged);

    // handle LOCPair case
    {
      RCP<LOCPair_GlobalEvaluationData> loc_pair = rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(ged);

      if(loc_pair!=Teuchos::null) {
        Teuchos::RCP<LinearObjContainer> loc = loc_pair->getGhostedLOC();
        // extract linear object container
        epetraContainer = rcp_dynamic_cast<EpetraLinearObjContainer>(loc);
      }
    }

    if(ro_ged!=Teuchos::null) {
      RCP<EpetraVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<EpetraVector_ReadOnly_GlobalEvaluationData>(ged,true);
  
      x_ = ro_ged->getGhostedVector_Epetra();
    }
    else if(epetraContainer!=Teuchos::null) {
      if (useTimeDerivativeSolutionVector_)
        x_ = epetraContainer->get_dxdt();
      else
        x_ = epetraContainer->get_x();
    }
  } // end "else if" contains(globalDataKey)

  // post condition
  TEUCHOS_ASSERT(x_!=Teuchos::null); // someone has to find the x_ vector

  // don't try to extract the dx, if its not required
  if(!secondApplySensitivities_) {
    dx_ = Teuchos::null; // do set it to null though!
    return;
  }

  // get second derivative perturbation
  ////////////////////////////////////////////////////////////
 
  // now parse the second derivative direction
  if(d.gedc.containsDataObject(sensitivities2ndPrefix_+globalDataKey_)) {
    ged = d.gedc.getDataObject(sensitivities2ndPrefix_+globalDataKey_);

    RCP<EpetraVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<EpetraVector_ReadOnly_GlobalEvaluationData>(ged,true);

    dx_ = ro_ged->getGhostedVector_Epetra();
  }

  // post conditions
  TEUCHOS_TEST_FOR_EXCEPTION(dx_==Teuchos::null,std::logic_error,
                             "Cannot find sensitivity vector associated with \""+sensitivities2ndPrefix_+globalDataKey_+"\" and \""+post+"\""); // someone has to find the dx_ vector
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_Epetra<panzer::Traits::Hessian, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   double seed_value = 0.0;
   if(firstApplySensitivities_) {
     if (useTimeDerivativeSolutionVector_ && gatherSeedIndex_<0) {
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
   }

   Epetra_Vector & x = *x_;

   // turn off sensitivies: this may be faster if we don't expand the term
   // but I suspect not because anywhere it is used the full complement of
   // sensitivies will be needed anyway.
   if(!firstApplySensitivities_)
      seed_value = 0.0;

   // Interface worksets handle DOFs from two element blocks. The derivative
   // offset for the other element block must be shifted by the derivative side
   // of my element block.
   int dos = 0;
   if (this->wda.getDetailsIndex() == 1) {
     // Get the DOF count for my element block.
     dos = globalIndexer_->getElementBlockGIDCount(workset.details(0).block_id);
   }

   // loop over the fields to be gathered
   for(std::size_t fieldIndex=0;
       fieldIndex<gatherFields_.size();fieldIndex++) {
      PHX::MDField<ScalarT,Cell,NODE> & field = gatherFields_[fieldIndex];
      int fieldNum = fieldIds_[fieldIndex];
      const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

      // gather operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

         const std::vector<int> & LIDs = globalIndexer_->getElementLIDs(cellLocalId);

         if(!firstApplySensitivities_) {
           // loop over basis functions and fill the fields
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs[offset];

             // set the value and seed the FAD object
             field(worksetCellIndex,basis) = x[lid];
           }
         }
         else {
           // loop over basis functions and fill the fields
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs[offset];

             // set the value and seed the FAD object
             field(worksetCellIndex,basis).val() = x[lid];
             field(worksetCellIndex,basis).fastAccessDx(dos + offset) = seed_value;
           }
         }

         // this is the direction to use for the second derivative
         if(secondApplySensitivities_) {
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs[offset];
             field(worksetCellIndex,basis).val().fastAccessDx(0) = (*dx_)[lid];
           }
         }
      }
   }
}

// **********************************************************************

} // end namespace panzer

#endif

#endif
