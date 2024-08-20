// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ScatterDirichletResidual_Epetra_Hessian_impl_hpp__
#define __Panzer_ScatterDirichletResidual_Epetra_Hessian_impl_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

#include "Panzer_PtrFromStlVector.hpp"

// the includes for this file come in as a result of the includes in the main 
// Epetra scatter dirichlet residual file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
ScatterDirichletResidual_Epetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
ScatterDirichletResidual_Epetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                                const Teuchos::RCP<const panzer::GlobalIndexer> & cIndexer,
                                const Teuchos::ParameterList& p) 
   : globalIndexer_(indexer)
   , colGlobalIndexer_(cIndexer) 
   , globalDataKey_("Residual Scatter Container")
   , preserveDiagonal_(false)
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

  if (p.isType<bool>("Preserve Diagonal"))
     preserveDiagonal_ = p.get<bool>("Preserve Diagonal");

  this->setName(scatterName+" Scatter Residual (Hessian)");
}
  
template<typename TRAITS,typename LO,typename GO>
void
ScatterDirichletResidual_Epetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  fieldIds_.resize(scatterFields_.size());

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].extent(1);
  num_eq = scatterFields_.size();
}

template<typename TRAITS,typename LO,typename GO>
void
ScatterDirichletResidual_Epetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d) 
{
  // extract linear object container
  epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject(globalDataKey_));
 
  if(epetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(loc,true);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<GlobalEvaluationData> dataContainer = d.gedc->getDataObject("Dirichlet Counter");
    Teuchos::RCP<EpetraLinearObjContainer> epetraContainer = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(dataContainer,true);

    dirichletCounter_ = epetraContainer->get_f();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }
}
  
template<typename TRAITS,typename LO,typename GO>
void
ScatterDirichletResidual_Epetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset) 
{
  using panzer::ptrFromStlVector;
  using std::vector;

   // TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
   //                           "ScatterDirichletResidual_Epetra<Hessian> is not yet implemented"); // just in case

  PHX::View<const int*, Kokkos::LayoutRight, PHX::Device> cLIDs, rLIDs;
   std::vector<double> jacRow;

   bool useColumnIndexer = colGlobalIndexer_!=Teuchos::null;
 
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<const EpetraLinearObjContainer> epetraContainer = epetraContainer_;
   TEUCHOS_ASSERT(epetraContainer!=Teuchos::null);
   Teuchos::RCP<Epetra_CrsMatrix> Jac = epetraContainer->get_A();
   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      rLIDs = globalIndexer_->getElementLIDs(cellLocalId); 
      if(useColumnIndexer)
        cLIDs = colGlobalIndexer_->getElementLIDs(cellLocalId); 
      else
        cLIDs = rLIDs;

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
   
         // this call "should" get the right ordering according to the Intrepid2 basis
         const std::pair<std::vector<int>,std::vector<int> > & indicePair 
               = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
         const std::vector<int> & elmtOffset = indicePair.first;
         const std::vector<int> & basisIdMap = indicePair.second;
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int row = rLIDs[offset];
            if(row<0) // not on this processor
               continue;

            int basisId = basisIdMap[basis];

            if (checkApplyBC_)
              if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
                continue;

            // zero out matrix row
            {
               int numEntries = 0;
               int * rowIndices = 0;
               double * rowValues = 0;

               Jac->ExtractMyRowView(row,numEntries,rowValues,rowIndices);

               for(int i=0;i<numEntries;i++) {
                  if(preserveDiagonal_) {
                    if(row!=rowIndices[i])
                      rowValues[i] = 0.0;
                  }
                  else
                    rowValues[i] = 0.0;
               }
            }
 
            // int gid = GIDs[offset];
            const ScalarT scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);

            if(dirichletCounter_!=Teuchos::null) {
              // std::cout << "Writing " << row << " " << dirichletCounter_->MyLength() << std::endl;
              (*dirichletCounter_)[row] = 1.0; // mark row as dirichlet
            }

            // loop over the sensitivity indices: all DOFs on a cell
            jacRow.resize(scatterField.size());
            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
              jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex).fastAccessDx(0);
    
            if(!preserveDiagonal_) {
              int err = Jac->ReplaceMyValues(row, cLIDs.size(), 
                ptrFromStlVector(jacRow), &cLIDs[0]);
              TEUCHOS_ASSERT(err==0); 
            }
         }
      }
   }
}
  
}

// **************************************************************
#endif

#endif
