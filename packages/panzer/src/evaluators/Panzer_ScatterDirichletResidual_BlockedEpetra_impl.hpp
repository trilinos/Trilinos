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
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Teuchos_FancyOStream.hpp"

#include <boost/unordered_map.hpp>

template<typename EvalT,typename Traits,typename LO,typename GO>
panzer::ScatterDirichletResidual_BlockedEpetra<EvalT, Traits,LO,GO>::
ScatterDirichletResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,int> > & indexer,
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
    PHX::MDField<ScalarT,Cell,NODE> scatterField = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterField);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder);

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
// Specialization: Residual
// **********************************************************************


template<typename Traits,typename LO,typename GO>
panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual, Traits,LO,GO>::
ScatterDirichletResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
                                const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
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
    scatterFields_[eq] = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template<typename Traits,typename LO,typename GO> 
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual, Traits,LO,GO>::
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

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
}

// **********************************************************************
template<typename Traits,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual, Traits,LO,GO>::
preEvaluate(typename Traits::PreEvalData d)
{
   typedef BlockedEpetraLinearObjContainer BLOC;

   // extract dirichlet counter from container
   Teuchos::RCP<BLOC> blockContainer 
         = Teuchos::rcp_dynamic_cast<BLOC>(d.getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_x(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
}

// **********************************************************************
template<typename Traits,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual, Traits,LO,GO>::
evaluateFields(typename Traits::EvalData workset)
{ 
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

   typedef BlockedEpetraLinearObjContainer BLOC;

   Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
   out.setShowProcRank(true);   
   out.setOutputToRootOnly(-1);   

   std::vector<std::pair<int,GO> > GIDs;
   std::vector<int> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   RCP<BLOC> blockedContainer = rcp_dynamic_cast<BLOC>(workset.ghostedLinContainer);
   RCP<ProductVectorBase<double> > r = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_f(),true);

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
         RCP<const Epetra_Map> r_map = blockedContainer->getMapForBlock(GIDs[i].first);

         LIDs[i] = r_map->LID(GIDs[i].second);
      }

      // std::vector<bool> is_owned(GIDs.size(), false);
      // globalIndexer_->ownedIndices(GIDs,is_owned);

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
   
         // this call "should" get the right ordering according to the Intrepid basis
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
            local_r[lid] = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);

            // record that you set a dirichlet condition
            local_dc[lid] = 1.0;
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits,typename LO,typename GO>
panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian, Traits,LO,GO>::
ScatterDirichletResidual_BlockedEpetra(const Teuchos::RCP<const BlockedDOFManager<LO,GO> > & indexer,
                                const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
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
    scatterFields_[eq] = PHX::MDField<ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template<typename Traits,typename LO,typename GO> 
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian, Traits,LO,GO>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  // globalIndexer_ = d.globalIndexer_;

  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
  num_eq = scatterFields_.size();
}

// **********************************************************************
template<typename Traits,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian, Traits,LO,GO>::
preEvaluate(typename Traits::PreEvalData d)
{
   typedef BlockedEpetraLinearObjContainer BLOC;

   // extract dirichlet counter from container
   Teuchos::RCP<BLOC> blockContainer 
         = Teuchos::rcp_dynamic_cast<BLOC>(d.getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_x(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
}

// **********************************************************************
template<typename Traits,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian, Traits,LO,GO>::
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

   typedef BlockedEpetraLinearObjContainer BLOC;

   std::vector<std::pair<int,GO> > GIDs;
   std::vector<int> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   RCP<BLOC> blockedContainer = rcp_dynamic_cast<BLOC>(workset.ghostedLinContainer);
   RCP<ProductVectorBase<double> > r = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_f());
   Teuchos::RCP<BlockedLinearOpBase<double> > Jac = rcp_dynamic_cast<BlockedLinearOpBase<double> >(blockedContainer->get_A());

   int numFieldBlocks = globalIndexer_->getNumFieldBlocks();
   std::vector<int> blockOffsets(numFieldBlocks+1); // number of fields, plus a sentinnel
   for(int blk=0;blk<numFieldBlocks;blk++) {
      int blockOffset = globalIndexer_->getBlockGIDOffset(blockId,blk);
      blockOffsets[blk] = blockOffset;
   }

   boost::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsMatrix> > jacEpetraBlocks;

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
         RCP<const Epetra_Map> r_map = blockedContainer->getMapForBlock(GIDs[i].first);

         LIDs[i] = r_map->LID(GIDs[i].second);
      }

      std::vector<bool> is_owned(GIDs.size(), false);
      globalIndexer_->ownedIndices(GIDs,is_owned);

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
   
         // this call "should" get the right ordering accordint to the Intrepid basis
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

            // zero out matrix row
            for(int blockColIndex=0;blockColIndex<numFieldBlocks;blockColIndex++) {
               int start = blockOffsets[blockColIndex];
               int end = blockOffsets[blockColIndex+1];

               if(end-start<=0) 
                  continue;

               // check hash table for jacobian sub block
               std::pair<int,int> blockIndex = std::make_pair(blockRowIndex,blockColIndex);
               Teuchos::RCP<Epetra_CrsMatrix> subJac = jacEpetraBlocks[blockIndex];

               // if you didn't find one before, add it to the hash table
               if(subJac==Teuchos::null) {
                  Teuchos::RCP<Epetra_Operator> eOp
                     = Thyra::get_Epetra_Operator(*Jac->getNonconstBlock(blockIndex.first,blockIndex.second)); 
                  subJac = rcp_dynamic_cast<Epetra_CrsMatrix>(eOp,true);
                  jacEpetraBlocks[blockIndex] = subJac;
               }

               int numEntries = 0;
               int * rowIndices = 0;
               double * rowValues = 0;

               subJac->ExtractMyRowView(lid,numEntries,rowValues,rowIndices);

               for(int i=0;i<numEntries;i++)
                  rowValues[i] = 0.0;
            }
 
            int basisId = basisIdMap[basis];
            const ScalarT & scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);
    
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
               Teuchos::RCP<Epetra_CrsMatrix> subJac = jacEpetraBlocks[blockIndex];

               // if you didn't find one before, add it to the hash table
               if(subJac==Teuchos::null) {
                  Teuchos::RCP<Epetra_Operator> eOp
                     = Thyra::get_Epetra_Operator(*Jac->getNonconstBlock(blockIndex.first,blockIndex.second)); 
                  subJac = rcp_dynamic_cast<Epetra_CrsMatrix>(eOp,true);
                  jacEpetraBlocks[blockIndex] = subJac;
               }

               // Sum Jacobian
               int err = subJac->ReplaceMyValues(lid, end-start, &jacRow[start],&LIDs[start]);
               TEUCHOS_ASSERT_EQUALITY(err,0);
            }
         }
      }
   }
}

// **********************************************************************

#endif
