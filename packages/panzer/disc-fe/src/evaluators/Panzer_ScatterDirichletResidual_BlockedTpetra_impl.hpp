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

#ifndef PANZER_SCATTER_DIRICHLET_RESIDUAL_BLOCEDEPETRA_IMPL_HPP
#define PANZER_SCATTER_DIRICHLET_RESIDUAL_BLOCEDEPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Teuchos_FancyOStream.hpp"

#include <unordered_map>

template <typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_BlockedTpetra<EvalT,TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & /* indexer */,
                                       const Teuchos::ParameterList& p)
{ 
  std::string scatterName = p.get<std::string>("Scatter Name");
  Teuchos::RCP<PHX::FieldTag> scatterHolder = 
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<const ScalarT,Cell,NODE> scatterField = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterField.fieldTag());
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder);

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
// Specialization: Residual
// **********************************************************************


template <typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
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

  // determine if we are scattering an initial condition
  scatterIC_ = p.isParameter("Scatter Initial Condition") ? p.get<bool>("Scatter Initial Condition") : false;

  Teuchos::RCP<PHX::DataLayout> dl = (!scatterIC_) ?
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional :
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;
  if (!scatterIC_) {
    side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
    local_side_id_ = p.get<int>("Local Side ID");
  }
  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.isParameter("Check Apply BC") ? p.get<bool>("Check Apply BC") : false;
  if (checkApplyBC_) {
    applyBC_.resize(names.size());
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<const bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData /* d */, 
                      PHX::FieldManager<TRAITS>& fm)
{
  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);

    if (checkApplyBC_)
      this->utils.setFieldData(applyBC_[fd],fm);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   // extract dirichlet counter from container
   Teuchos::RCP<const ContainerType> blockContainer 
         = Teuchos::rcp_dynamic_cast<ContainerType>(d.gedc->getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_f(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));

   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(blockedContainer_));
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
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

   RCP<ProductVectorBase<double> > r = (!scatterIC_) ? 
     rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_f(),true) :
     rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_x(),true);

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      globalIndexer_->getElementGIDs(cellLocalId,GIDs); 

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++) {
         // used for doing local ID lookups
         RCP<const MapType> r_map = blockedContainer_->getMapForBlock(GIDs[i].first);

         LIDs[i] = r_map->getLocalElement(GIDs[i].second);
      }

      // loop over each field to be scattered
      Teuchos::ArrayRCP<double> local_r;
      Teuchos::ArrayRCP<double> local_dc;
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         int indexerId = globalIndexer_->getFieldBlock(fieldNum);

         RCP<SpmdVectorBase<double> > dc = rcp_dynamic_cast<SpmdVectorBase<double> >(dirichletCounter_->getNonconstVectorBlock(indexerId));
         dc->getNonconstLocalData(ptrFromRef(local_dc));

         // grab local data for inputing
         RCP<SpmdVectorBase<double> > block_r = rcp_dynamic_cast<SpmdVectorBase<double> >(r->getNonconstVectorBlock(indexerId));
         block_r->getNonconstLocalData(ptrFromRef(local_r));

         if (!scatterIC_) {
           // this call "should" get the right ordering according to the Intrepid2 basis
           const std::pair<std::vector<int>,std::vector<int> > & indicePair 
             = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
           const std::vector<int> & elmtOffset = indicePair.first;
           const std::vector<int> & basisIdMap = indicePair.second;

           // loop over basis functions
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs[offset];
             if(lid<0) // not on this processor!
               continue;

             int basisId = basisIdMap[basis];

             if (checkApplyBC_)
               if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
                 continue;

             local_r[lid] = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);

             // record that you set a dirichlet condition
             local_dc[lid] = 1.0;
           }
         } else {
           // this call "should" get the right ordering according to the Intrepid2 basis
           const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

           // loop over basis functions
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs[offset];
             if(lid<0) // not on this processor!
               continue;

            local_r[lid] = (scatterFields_[fieldIndex])(worksetCellIndex,basis);

            // record that you set a dirichlet condition
            local_dc[lid] = 1.0;
           }
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template <typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
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
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
  local_side_id_ = p.get<int>("Local Side ID");
  
  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.get<bool>("Check Apply BC");
  if (checkApplyBC_) {
    applyBC_.resize(names.size());
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<const bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& fm)
{
  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);

    if (checkApplyBC_)
      this->utils.setFieldData(applyBC_[fd],fm);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
  num_eq = scatterFields_.size();
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   // extract dirichlet counter from container
   Teuchos::RCP<const ContainerType> blockContainer 
         = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_f(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));

   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(blockedContainer_));
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
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
 
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   RCP<ProductVectorBase<double> > r = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_f());
   Teuchos::RCP<BlockedLinearOpBase<double> > Jac = rcp_dynamic_cast<BlockedLinearOpBase<double> >(blockedContainer_->get_A());

   int numFieldBlocks = globalIndexer_->getNumFieldBlocks();
   std::vector<int> blockOffsets(numFieldBlocks+1); // number of fields, plus a sentinnel
   for(int blk=0;blk<numFieldBlocks;blk++) {
      int blockOffset = globalIndexer_->getBlockGIDOffset(blockId,blk);
      blockOffsets[blk] = blockOffset;
   }

   std::unordered_map<std::pair<int,int>,Teuchos::RCP<CrsMatrixType>,panzer::pair_hash> jacTpetraBlocks;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      globalIndexer_->getElementGIDs(cellLocalId,GIDs); 
      blockOffsets[numFieldBlocks] = GIDs.size();

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++) {
         // used for doing local ID lookups
         RCP<const MapType> r_map = blockedContainer_->getMapForBlock(GIDs[i].first);

         LIDs[i] = r_map->getLocalElement(GIDs[i].second);
      }

      // loop over each field to be scattered
      Teuchos::ArrayRCP<double> local_r, local_dc;
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         int blockRowIndex = globalIndexer_->getFieldBlock(fieldNum);

         RCP<SpmdVectorBase<double> > dc = rcp_dynamic_cast<SpmdVectorBase<double> >(dirichletCounter_->getNonconstVectorBlock(blockRowIndex));
         dc->getNonconstLocalData(ptrFromRef(local_dc));

         // grab local data for inputing
         RCP<SpmdVectorBase<double> > block_r = rcp_dynamic_cast<SpmdVectorBase<double> >(r->getNonconstVectorBlock(blockRowIndex));
         block_r->getNonconstLocalData(ptrFromRef(local_r));
   
         // this call "should" get the right ordering according to the Intrepid2 basis
         const std::pair<std::vector<int>,std::vector<int> > & indicePair 
               = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
         const std::vector<int> & elmtOffset = indicePair.first;
         const std::vector<int> & basisIdMap = indicePair.second;
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            if(lid<0) // not on this processor
               continue;

            int basisId = basisIdMap[basis];
            
            if (checkApplyBC_)
              if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
                continue;

            // zero out matrix row
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

               std::size_t sz = subJac->getNumEntriesInLocalRow(lid);
               std::size_t numEntries = 0;
               Teuchos::Array<LO> rowIndices(sz);
               Teuchos::Array<double> rowValues(sz);

               subJac->getLocalRowCopy(lid,rowIndices,rowValues,numEntries);

               for(std::size_t i=0;i<numEntries;i++)
                  rowValues[i] = 0.0;

               subJac->replaceLocalValues(lid,rowIndices,rowValues);
            }
 
            const ScalarT scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);
    
            local_r[lid] = scatterField.val();
            local_dc[lid] = 1.0; // mark row as dirichlet
    
            // loop over the sensitivity indices: all DOFs on a cell
            std::vector<double> jacRow(scatterField.size(),0.0);
    
            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);
            TEUCHOS_ASSERT(jacRow.size()==GIDs.size());
    
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
               subJac->replaceLocalValues(lid, Teuchos::arrayViewFromVector(LIDs).view(start,end-start), 
                                               Teuchos::arrayViewFromVector(jacRow).view(start,end-start));
            }
         }
      }
   }
}

// **********************************************************************

#endif
