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

#ifndef PANZER_GATHER_ORIENTATION_IMPL_HPP
#define PANZER_GATHER_ORIENTATION_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"
#include "Panzer_PureBasis.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename TRAITS,typename LO,typename GO>
panzer::GatherOrientation<EvalT, TRAITS,LO,GO>::
GatherOrientation(
  const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & indexer,
  const Teuchos::ParameterList& p)
{ 
  indexers_.push_back(indexer);

  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  indexerNames_ = p.get< Teuchos::RCP< std::vector<std::string> > >("Indexer Names");

  // this is beging to fix the issues with incorrect use of const
  Teuchos::RCP<const panzer::PureBasis> basis;
  if(p.isType< Teuchos::RCP<panzer::PureBasis> >("Basis"))
    basis = p.get< Teuchos::RCP<panzer::PureBasis> >("Basis");
  else
    basis = p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis");

  gatherFieldOrientations_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFieldOrientations_[fd] = 
      // PHX::MDField<ScalarT,Cell,NODE>(names[fd]+" Orientation",basis->functional);
      PHX::MDField<ScalarT,Cell,NODE>(basis->name()+" Orientation",basis->functional);
    this->addEvaluatedField(gatherFieldOrientations_[fd]);
  }

  this->setName("Gather Orientation");
}

template<typename EvalT,typename TRAITS,typename LO,typename GO>
panzer::GatherOrientation<EvalT, TRAITS,LO,GO>::
GatherOrientation(const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LO,GO> > > & indexers,
                  const Teuchos::ParameterList& p)
  : indexers_(indexers)
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

  gatherFieldOrientations_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFieldOrientations_[fd] = 
      PHX::MDField<ScalarT,Cell,NODE>(basis->name()+" Orientation",basis->functional);
    this->addEvaluatedField(gatherFieldOrientations_[fd]);
  }

  this->setName("Gather Orientation");
}

// **********************************************************************
template<typename EvalT,typename TRAITS,typename LO,typename GO>
void panzer::GatherOrientation<EvalT, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData /* d */, 
		      PHX::FieldManager<TRAITS>& /* fm */)
{
  TEUCHOS_ASSERT(gatherFieldOrientations_.size() == indexerNames_->size());

  indexerIds_.resize(gatherFieldOrientations_.size());
  subFieldIds_.resize(gatherFieldOrientations_.size());

  for (std::size_t fd = 0; fd < gatherFieldOrientations_.size(); ++fd) {
    // get field ID from DOF manager
    const std::string& fieldName = (*indexerNames_)[fd];

    indexerIds_[fd]  = getFieldBlock(fieldName,indexers_);
    subFieldIds_[fd] = indexers_[indexerIds_[fd]]->getFieldNum(fieldName);
  }

  indexerNames_ = Teuchos::null;  // Don't need this anymore
}

// **********************************************************************
template<typename EvalT,typename TRAITS,typename LO,typename GO>
void panzer::GatherOrientation<EvalT, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   std::vector<double> orientation;
 
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // loop over the fields to be gathered
   for (std::size_t fieldIndex=0; fieldIndex<gatherFieldOrientations_.size();fieldIndex++) {

      int indexerId   = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      auto subRowIndexer = indexers_[indexerId];
      const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);
 
      // gather operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];
 
         subRowIndexer->getElementOrientation(cellLocalId,orientation); 

         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            (gatherFieldOrientations_[fieldIndex])(worksetCellIndex,basis) = orientation[offset];
            (gatherFieldOrientations_[fieldIndex])(worksetCellIndex,basis) = std::sqrt(-1.0);
         }
      }
   }
}

#endif
