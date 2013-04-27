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

#ifndef PANZER_SCATTER_RESIDUAL_BLOCKEDEPETRA_IMPL_HPP
#define PANZER_SCATTER_RESIDUAL_BLOCKEDEPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Teuchos_FancyOStream.hpp"

template <typename EvalT,typename Traits,typename S,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_BlockedTpetra<EvalT, Traits,S,LO,GO,NodeT>::
ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
                              const Teuchos::ParameterList& p)
{ 
  std::string scatterName = p.get<std::string>("Scatter Name");
  Teuchos::RCP<PHX::FieldTag> scatterHolder =
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;
  
  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<ScalarT,Cell,NODE> field = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(field);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder);

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Residual, Traits,S,LO,GO,NodeT>::
ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
                              const Teuchos::ParameterList& p)
  : globalIndexer_(indexer) 
  , globalDataKey_("Residual Scatter Container")
{ 
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ = 
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;
  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Residual, Traits,S,LO,GO,NodeT>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  fieldIds_.resize(scatterFields_.size());

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Residual, Traits,S,LO,GO,NodeT>::
preEvaluate(typename Traits::PreEvalData d)
{
   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.getDataObject(globalDataKey_),true);
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Residual, Traits,S,LO,GO,NodeT>::
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
   std::vector<int> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   RCP<const ContainerType> blockedContainer = blockedContainer_;

   RCP<ProductVectorBase<double> > r = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_f(),true);

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      globalIndexer_->getElementGIDs(cellLocalId,GIDs,blockId); 

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++) {
         // used for doing local ID lookups
         RCP<const MapType> r_map = blockedContainer->getMapForBlock(GIDs[i].first);

         LIDs[i] = r_map->getLocalElement(GIDs[i].second);
      }

      // loop over each field to be scattered
      Teuchos::ArrayRCP<double> local_r;
      for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         int indexerId = globalIndexer_->getFieldBlock(fieldNum);

         // grab local data for inputing
         RCP<SpmdVectorBase<double> > block_r = rcp_dynamic_cast<SpmdVectorBase<double> >(r->getNonconstVectorBlock(indexerId));
         block_r->getNonconstLocalData(ptrFromRef(local_r));

         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            local_r[lid] += (scatterFields_[fieldIndex])(worksetCellIndex,basis);
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian, Traits,S,LO,GO,NodeT>::
ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
                              const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
   , globalDataKey_("Residual Scatter Container")
{ 
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ = 
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;
  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian, Traits,S,LO,GO,NodeT>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  fieldIds_.resize(scatterFields_.size());

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian, Traits,S,LO,GO,NodeT>::
preEvaluate(typename Traits::PreEvalData d)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   // extract linear object container
   blockedContainer_ = rcp_dynamic_cast<const ContainerType>(d.getDataObject(globalDataKey_));

   if(blockedContainer_==Teuchos::null) {
     RCP<const LOCPair_GlobalEvaluationData> gdata = rcp_dynamic_cast<const LOCPair_GlobalEvaluationData>(d.getDataObject(globalDataKey_),true);
     blockedContainer_ = rcp_dynamic_cast<const ContainerType>(gdata->getGhostedLOC());
   }
}

// **********************************************************************
template <typename Traits,typename S,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian, Traits,S,LO,GO,NodeT>::
evaluateFields(typename Traits::EvalData workset)
{ 
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;
   using Thyra::BlockedLinearOpBase;

   std::vector<std::pair<int,GO> > GIDs;
   std::vector<LO> LIDs;
   std::vector<double> jacRow;

   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   RCP<const ContainerType> blockedContainer = blockedContainer_;

   RCP<ProductVectorBase<double> > r = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_f());
   Teuchos::RCP<BlockedLinearOpBase<double> > Jac = rcp_dynamic_cast<BlockedLinearOpBase<double> >(blockedContainer->get_A());

   int numFieldBlocks = globalIndexer_->getNumFieldBlocks();
   std::vector<int> blockOffsets(numFieldBlocks+1); // number of fields, plus a sentinnel
   for(int blk=0;blk<numFieldBlocks;blk++) {
      int blockOffset = globalIndexer_->getBlockGIDOffset(blockId,blk);
      blockOffsets[blk] = blockOffset;
   }

   boost::unordered_map<std::pair<int,int>,Teuchos::RCP<CrsMatrixType> > jacTpetraBlocks;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets" may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      globalIndexer_->getElementGIDs(cellLocalId,GIDs,blockId); 

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++) {
         // used for doing local ID lookups
         RCP<const MapType> r_map = blockedContainer->getMapForBlock(GIDs[i].first);

         LIDs[i] = r_map->getLocalElement(GIDs[i].second);
      }

      // loop over each field to be scattered
      Teuchos::ArrayRCP<double> local_r;
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         int blockRowIndex = globalIndexer_->getFieldBlock(fieldNum);

         // grab local data for inputing
         if(r!=Teuchos::null) {
            RCP<SpmdVectorBase<double> > block_r = rcp_dynamic_cast<SpmdVectorBase<double> >(r->getNonconstVectorBlock(blockRowIndex));
            block_r->getNonconstLocalData(ptrFromRef(local_r));
         }

         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
        
         // loop over the basis functions (currently they are nodes)
         for(std::size_t rowBasisNum = 0; rowBasisNum < elmtOffset.size(); rowBasisNum++) {
            const ScalarT & scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,rowBasisNum);
            int rowOffset = elmtOffset[rowBasisNum];
            int r_lid = LIDs[rowOffset];
    
            // Sum residual
            if(local_r!=Teuchos::null)
               local_r[r_lid] += (scatterField.val());

            blockOffsets[numFieldBlocks] = scatterField.size(); // add the sentinel
    
            // loop over the sensitivity indices: all DOFs on a cell
            jacRow.resize(scatterField.size());
            
            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex) {
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);
            }
    
            for(int blockColIndex=0;blockColIndex<numFieldBlocks;blockColIndex++) {
               int start = blockOffsets[blockColIndex];
               int end = blockOffsets[blockColIndex+1];

               if(end-start<=0) 
                  continue;

               // check hash table for jacobian sub block
               std::pair<int,int> blockIndex = std::make_pair(blockRowIndex,blockColIndex);
               Teuchos::RCP<CrsMatrixType> subJac = jacTpetraBlocks[blockIndex];

               // if you didn't find one before, add it to the hash table
               if(subJac==Teuchos::null) {
                  Teuchos::RCP<Thyra::LinearOpBase<double> > tOp = Jac->getNonconstBlock(blockIndex.first,blockIndex.second); 

                  // block operator is null, don't do anything (it is excluded)
                  if(Teuchos::is_null(tOp))
                     continue;

                  Teuchos::RCP<OperatorType> tpetra_Op = rcp_dynamic_cast<ThyraLinearOp>(tOp)->getTpetraOperator();
                  subJac = rcp_dynamic_cast<CrsMatrixType>(tpetra_Op,true);
                  jacTpetraBlocks[blockIndex] = subJac;
               }

               // Sum Jacobian
               subJac->sumIntoLocalValues(r_lid, Teuchos::arrayViewFromVector(LIDs).view(start,end-start),
                                                 Teuchos::arrayViewFromVector(jacRow).view(start,end-start));
            }
         } // end rowBasisNum
      } // end fieldIndex
   }
}

// **********************************************************************

#endif
