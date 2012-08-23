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

#ifndef PANZER_SCATTER_DIRICHLET_RESIDUAL_EPETRA_SG_IMPL_HPP
#define PANZER_SCATTER_DIRICHLET_RESIDUAL_EPETRA_SG_IMPL_HPP

#include "Panzer_config.hpp"
#ifdef HAVE_STOKHOS

#include "Panzer_SGEpetraLinearObjContainer.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// **********************************************************************
// Specialization: SGResidual
// **********************************************************************

template<typename Traits,typename LO,typename GO>
panzer::ScatterDirichletResidual_Epetra<panzer::Traits::SGResidual, Traits,LO,GO>::
ScatterDirichletResidual_Epetra(const Teuchos::RCP<const UniqueGlobalIndexer<LO,GO> > & indexer,
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

  this->setName(scatterName+" Scatter Residual (SGResidual)");
}

// **********************************************************************
template<typename Traits,typename LO,typename GO> 
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::SGResidual, Traits,LO,GO>::
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
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::SGResidual, Traits,LO,GO>::
preEvaluate(typename Traits::PreEvalData d)
{
   // extract dirichlet counter from container
   Teuchos::RCP<EpetraLinearObjContainer> epetraContainer 
         = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = epetraContainer->get_x();
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
}

// **********************************************************************
template<typename Traits,typename LO,typename GO>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::SGResidual, Traits,LO,GO>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<GO> GIDs;
   std::vector<int> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<SGEpetraLinearObjContainer> sgEpetraContainer 
         = Teuchos::rcp_dynamic_cast<SGEpetraLinearObjContainer>(workset.ghostedLinContainer);
   Teuchos::RCP<Epetra_Vector> r_template = (*sgEpetraContainer->begin())->get_f();
   const Epetra_BlockMap & map = r_template->Map();

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
      for(std::size_t i=0;i<GIDs.size();i++)
         LIDs[i] = map.LID(GIDs[i]);

      std::vector<bool> is_owned(GIDs.size(), false);
      globalIndexer_->ownedIndices(GIDs,is_owned);

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
   
         // this call "should" get the right ordering accordint to the Intrepid basis
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
            const ScalarT & field = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);

            // loop over stochastic basis scatter field values to residual vectors
            int stochIndex = 0;
            panzer::SGEpetraLinearObjContainer::iterator itr; 
            for(itr=sgEpetraContainer->begin();itr!=sgEpetraContainer->end();++itr,++stochIndex)
               (*(*itr)->get_f())[lid] = field.coeff(stochIndex);

            // dirichlet condition application
            (*dirichletCounter_)[lid] = 1.0;
         }
      }
   }
}

// **********************************************************************
// Specialization: SGJacobian
// **********************************************************************

template<typename Traits,typename LO,typename GO>
panzer::ScatterDirichletResidual_Epetra<panzer::Traits::SGJacobian, Traits,LO,GO>::
ScatterDirichletResidual_Epetra(const Teuchos::RCP<const UniqueGlobalIndexer<LO,GO> > & indexer,
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

  this->setName(scatterName+" Scatter Residual (SGJacobian)");
}

// **********************************************************************
template<typename Traits,typename LO,typename GO> 
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::SGJacobian, Traits,LO,GO>::
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
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::SGJacobian, Traits,LO,GO>::
preEvaluate(typename Traits::PreEvalData d)
{
   // extract dirichlet counter from container
   Teuchos::RCP<EpetraLinearObjContainer> epetraContainer 
         = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = epetraContainer->get_x();
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
}

// **********************************************************************
template<typename Traits,typename LO,typename GO>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::SGJacobian, Traits,LO,GO>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<GO> GIDs;
   std::vector<int> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   Teuchos::RCP<SGEpetraLinearObjContainer> sgEpetraContainer 
         = Teuchos::rcp_dynamic_cast<SGEpetraLinearObjContainer>(workset.ghostedLinContainer);
   Teuchos::RCP<Epetra_CrsMatrix> Jac_template = (*sgEpetraContainer->begin())->get_A();
   const Epetra_BlockMap & map = Jac_template->RowMap();

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
      for(std::size_t i=0;i<GIDs.size();i++)
         LIDs[i] = map.LID(GIDs[i]);

      std::vector<bool> is_owned(GIDs.size(), false);
      globalIndexer_->ownedIndices(GIDs,is_owned);

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
   
         // this call "should" get the right ordering accordint to the Intrepid basis
         const std::pair<std::vector<int>,std::vector<int> > & indicePair 
               = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
         const std::vector<int> & elmtOffset = indicePair.first;
         const std::vector<int> & basisIdMap = indicePair.second;
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            int gid = GIDs[offset];
            int basisId = basisIdMap[basis];
            if(lid<0) // not on this processor
               continue;

            const ScalarT & scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);

            // loop over stochastic basis scatter field values to residual vectors
            int stochIndex = 0;
            panzer::SGEpetraLinearObjContainer::iterator itr; 
            for(itr=sgEpetraContainer->begin();itr!=sgEpetraContainer->end();++itr,++stochIndex) {
               Teuchos::RCP<Epetra_Vector> r = (*itr)->get_f();
               Teuchos::RCP<Epetra_CrsMatrix> Jac = (*itr)->get_A();

               // zero out matrix row
               {
                  int numEntries = 0;
                  int * rowIndices = 0;
                  double * rowValues = 0;
   
                  Jac->ExtractMyRowView(lid,numEntries,rowValues,rowIndices);
   
                  for(int i=0;i<numEntries;i++)
                     rowValues[i] = 0.0;
               }
    
               (*r)[lid] = scatterField.val().coeff(stochIndex);
       
               // loop over the sensitivity indices: all DOFs on a cell
               std::vector<double> jacRow(scatterField.size(),0.0);
       
               for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
                  jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex).coeff(stochIndex);
       
               int err = Jac->ReplaceGlobalValues(gid, scatterField.size(), &jacRow[0],&GIDs[0]);
               TEUCHOS_ASSERT(err==0); 
            }

            // dirichlet condition application
            (*dirichletCounter_)[lid] = 1.0;
         }
      }
   }
}

// **********************************************************************
#endif

#endif
