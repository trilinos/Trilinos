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

#ifndef __Panzer_GatherSolution_BlockedEpetra_Hessian_impl_hpp__
#define __Panzer_GatherSolution_BlockedEpetra_Hessian_impl_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_ParameterList_GlobalEvaluationData.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"

namespace panzer {

// **********************************************************************
// Specialization: Hessian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS,LO,GO>::
GatherSolution_BlockedEpetra(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,int> > > & indexers,
                             const Teuchos::ParameterList& p)
  : indexers_(indexers)
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
    std::string n = "GatherSolution (BlockedEpetra, No First Sensitivities): "+firstName+" (Hessian)";
    this->setName(n);
  }
  else {
    std::string n = "GatherSolution (BlockedEpetra): "+firstName+" ("+PHX::typeAsString<EvalT>()+") ";
    this->setName(n);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS,LO,GO>::
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

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  typedef BlockedEpetraLinearObjContainer BLOC;

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

    RCP<BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<BlockedVector_ReadOnly_GlobalEvaluationData>(ged,true);

    x_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(ro_ged->getGhostedVector());

    // post condition
    TEUCHOS_TEST_FOR_EXCEPTION(x_==Teuchos::null,std::logic_error,
                               "Gather Hessian: Can't find x vector in GEDC \"" << globalDataKey_ << post << "\""
                               "A cast failed for " << ged << ". Type is " << Teuchos::typeName(*ged) << std::endl); 
  }
  else if(d.gedc.containsDataObject(globalDataKey_)) {
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

    // post condition
    TEUCHOS_TEST_FOR_EXCEPTION(x_==Teuchos::null,std::logic_error,
                               "Gather Hessian: Can't find x vector in GEDC \"" << globalDataKey_ << "\" (" << post << "). "
                               "A cast failed for " << ged << ". Type is " << Teuchos::typeName(*ged)); 
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

    RCP<BlockedVector_ReadOnly_GlobalEvaluationData> ro_ged = rcp_dynamic_cast<BlockedVector_ReadOnly_GlobalEvaluationData>(ged,true);

    dx_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(ro_ged->getGhostedVector());
  }

  // post conditions
  TEUCHOS_TEST_FOR_EXCEPTION(dx_==Teuchos::null,std::logic_error,
                             "Cannot find sensitivity vector associated with \""+sensitivities2ndPrefix_+globalDataKey_+"\" and \""+post+"\""); // someone has to find the dx_ vector
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::GatherSolution_BlockedEpetra<panzer::Traits::Hessian, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::ptrFromRef;
   using Teuchos::ArrayRCP;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

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

   // turn off sensitivies: this may be faster if we don't expand the term
   // but I suspect not because anywhere it is used the full complement of
   // sensitivies will be needed anyway.
   if(!firstApplySensitivities_)
      seed_value = 0.0;

   std::vector<int> blockOffsets;
   computeBlockOffsets(blockId,indexers_,blockOffsets);

   // loop over the fields to be gathered
   for(std::size_t fieldIndex=0;
       fieldIndex<gatherFields_.size();fieldIndex++) {

      PHX::MDField<ScalarT,Cell,NODE> & field = gatherFields_[fieldIndex];

      int indexerId   = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // grab local data for inputing
      Teuchos::ArrayRCP<const double> local_x, local_dx;
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_->getNonconstVectorBlock(indexerId))->getLocalData(ptrFromRef(local_x));
      if(secondApplySensitivities_)
        rcp_dynamic_cast<SpmdVectorBase<double> >(dx_->getNonconstVectorBlock(indexerId))->getLocalData(ptrFromRef(local_dx));

      auto subRowIndexer = indexers_[indexerId];
      const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

      int startBlkOffset = blockOffsets[indexerId];

      // gather operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

         const std::vector<int> & LIDs = subRowIndexer->getElementLIDs(cellLocalId);

         if(!firstApplySensitivities_) {
           // loop over basis functions and fill the fields
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs[offset];

             // set the value and seed the FAD object
             field(worksetCellIndex,basis) = local_x[lid];
           }
         }
         else {
           // loop over basis functions and fill the fields
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs[offset];

             // set the value and seed the FAD object
             field(worksetCellIndex,basis).val() = local_x[lid];
             field(worksetCellIndex,basis).fastAccessDx(startBlkOffset+offset) = seed_value;
           }
         }

         // this is the direction to use for the second derivative
         if(secondApplySensitivities_) {
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs[offset];
             field(worksetCellIndex,basis).val().fastAccessDx(0) = local_dx[lid];
           }
         }
      }
   }
}

// **********************************************************************

} // end namespace panzer

#endif

#endif
