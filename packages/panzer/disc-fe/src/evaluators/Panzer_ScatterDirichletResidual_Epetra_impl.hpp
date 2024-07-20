// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCATTER_DIRICHLET_RESIDUAL_EPETRA_IMPL_HPP
#define PANZER_SCATTER_DIRICHLET_RESIDUAL_EPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_ParameterList_GlobalEvaluationData.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Teuchos_FancyOStream.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************


template<typename TRAITS,typename LO,typename GO>
panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Residual, TRAITS,LO,GO>::
ScatterDirichletResidual_Epetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                                const Teuchos::RCP<const panzer::GlobalIndexer> & /* cIndexer */,
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
template<typename TRAITS,typename LO,typename GO> 
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Residual, TRAITS,LO,GO>::
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
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Residual, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // extract linear object container
  epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject(globalDataKey_));
 
  if(epetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(loc);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<EpetraLinearObjContainer> epetraContainer 
       = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject("Dirichlet Counter"),true);

    dirichletCounter_ = epetraContainer->get_f();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Residual, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<const EpetraLinearObjContainer> epetraContainer = epetraContainer_;
   Teuchos::RCP<Epetra_Vector> r = (!scatterIC_) ? 
     epetraContainer->get_f() : 
     epetraContainer->get_x();
   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   auto LIDs = globalIndexer_->getLIDs();
   auto LIDs_h = Kokkos::create_mirror_view(LIDs);
   Kokkos::deep_copy(LIDs_h, LIDs);


   // loop over each field to be scattered
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
     int fieldNum = fieldIds_[fieldIndex];
     auto field = PHX::as_view(scatterFields_[fieldIndex]);
     auto field_h = Kokkos::create_mirror_view(field);
     Kokkos::deep_copy(field_h, field);

     for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
       std::size_t cellLocalId = localCellIds[worksetCellIndex];

       if (!scatterIC_) {
	 // this call "should" get the right ordering according to the Intrepid2 basis
	 const std::pair<std::vector<int>,std::vector<int> > & indicePair 
	   = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
	 const std::vector<int> & elmtOffset = indicePair.first;
	 const std::vector<int> & basisIdMap = indicePair.second;
	 
	 // loop over basis functions
	 for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
	   int offset = elmtOffset[basis];
	   int lid = LIDs_h(cellLocalId, offset);
	   if(lid<0) // not on this processor!
	     continue;

	   int basisId = basisIdMap[basis];
	   
	   if (checkApplyBC_)
	     if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
	       continue;

	   (*r)[lid] = field_h(worksetCellIndex,basisId);
            
	   // record that you set a dirichlet condition
	   if(dirichletCounter_!=Teuchos::null)
	     (*dirichletCounter_)[lid] = 1.0; 
	 }
       } else {
	 // this call "should" get the right ordering according to the Intrepid2 basis
	 const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
	 
           // loop over basis functions
	 for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
	   int offset = elmtOffset[basis];
	   int lid = LIDs_h(cellLocalId, offset);
	   if(lid<0) // not on this processor!
	     continue;

	   (*r)[lid] = field_h(worksetCellIndex,basis);
            
	   // record that you set a dirichlet condition
	   if(dirichletCounter_!=Teuchos::null)
	     (*dirichletCounter_)[lid] = 1.0; 
	 }
       }
     }
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************


template<typename TRAITS,typename LO,typename GO>
panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
ScatterDirichletResidual_Epetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                                const Teuchos::RCP<const panzer::GlobalIndexer> & /* cIndexer */,
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

  this->setName(scatterName+" Scatter Tangent");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO> 
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
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
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // extract linear object container
  epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject(globalDataKey_));
 
  if(epetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(loc);

    dirichletCounter_ = Teuchos::null;
  }
  else {
    // extract dirichlet counter from container
    Teuchos::RCP<EpetraLinearObjContainer> epetraContainer 
       = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject("Dirichlet Counter"),true);

    dirichletCounter_ = epetraContainer->get_f();
    TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));
  }

  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // this is the list of parameters and their names that this scatter has to account for
  std::vector<std::string> activeParameters =
    rcp_dynamic_cast<ParameterList_GlobalEvaluationData>(d.gedc->getDataObject("PARAMETER_NAMES"))->getActiveParameters();

  // ETP 02/03/16:  This code needs to be updated to properly handle scatterIC_
  TEUCHOS_ASSERT(!scatterIC_);
  dfdp_vectors_.clear();
  for(std::size_t i=0;i<activeParameters.size();i++) {
    RCP<Epetra_Vector> vec = rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject(activeParameters[i]),true)->get_f();
    dfdp_vectors_.push_back(vec);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<const EpetraLinearObjContainer> epetraContainer = epetraContainer_;
   Teuchos::RCP<Epetra_Vector> r = (!scatterIC_) ? 
     epetraContainer->get_f() : 
     epetraContainer->get_x();
   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // loop over each field to be scattered
   auto LIDs = globalIndexer_->getLIDs();
   auto LIDs_h = Kokkos::create_mirror_view(LIDs);
   Kokkos::deep_copy(LIDs_h, LIDs);
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
     int fieldNum = fieldIds_[fieldIndex];
     auto scatterField_h = Kokkos::create_mirror_view(scatterFields_[fieldIndex].get_static_view());
     Kokkos::deep_copy(scatterField_h, scatterFields_[fieldIndex].get_static_view());

     // scatter operation for each cell in workset
     for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
       std::size_t cellLocalId = localCellIds[worksetCellIndex];

         if (!scatterIC_) {
           // this call "should" get the right ordering according to the Intrepid2 basis
           const std::pair<std::vector<int>,std::vector<int> > & indicePair 
             = globalIndexer_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
           const std::vector<int> & elmtOffset = indicePair.first;
           const std::vector<int> & basisIdMap = indicePair.second;
   
           // loop over basis functions
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs_h(cellLocalId, offset);
             if(lid<0) // not on this processor!
               continue;
             
             int basisId = basisIdMap[basis];

             if (checkApplyBC_)
               if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
                 continue;

             ScalarT value = scatterField_h(worksetCellIndex,basisId);
             // (*r)[lid] = (scatterFields_[fieldIndex])(worksetCellIndex,basisId).val();

             // then scatter the sensitivity vectors
             if(value.size()==0) 
               for(std::size_t d=0;d<dfdp_vectors_.size();d++)
                 (*dfdp_vectors_[d])[lid] = 0.0;
             else
               for(int d=0;d<value.size();d++) {
                 (*dfdp_vectors_[d])[lid] = value.fastAccessDx(d);
               }
            
             // record that you set a dirichlet condition
             if(dirichletCounter_!=Teuchos::null)
               (*dirichletCounter_)[lid] = 1.0;
           }
         } else {
           // this call "should" get the right ordering according to the Intrepid2 basis
           const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

           // loop over basis functions
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs_h(cellLocalId, offset);
             if(lid<0) // not on this processor!
               continue;
             
             ScalarT value = scatterField_h(worksetCellIndex,basis);
             // (*r)[lid] = (scatterFields_[fieldIndex])(worksetCellIndex,basis).val();

             // then scatter the sensitivity vectors
             if(value.size()==0) 
               for(std::size_t d=0;d<dfdp_vectors_.size();d++)
                 (*dfdp_vectors_[d])[lid] = 0.0;
             else
               for(int d=0;d<value.size();d++) {
                 (*dfdp_vectors_[d])[lid] = value.fastAccessDx(d);
               }
            
             // record that you set a dirichlet condition
             if(dirichletCounter_!=Teuchos::null)
               (*dirichletCounter_)[lid] = 1.0;
           }
         }
     }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
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

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO> 
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
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

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
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

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
   Kokkos::View<const int*, Kokkos::LayoutRight, PHX::Device> cLIDs, rLIDs;
   int gidCount(0);
   bool useColumnIndexer = colGlobalIndexer_!=Teuchos::null;
   const Teuchos::RCP<const panzer::GlobalIndexer>&
     colGlobalIndexer = useColumnIndexer ? colGlobalIndexer_ : globalIndexer_;
 
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<const EpetraLinearObjContainer> epetraContainer = epetraContainer_;
   TEUCHOS_ASSERT(epetraContainer!=Teuchos::null);
   Teuchos::RCP<Epetra_Vector> r = epetraContainer->get_f(); 
   Teuchos::RCP<Epetra_CrsMatrix> Jac = epetraContainer->get_A();
   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset

   auto LIDs = globalIndexer_->getLIDs();
   auto colLIDs = colGlobalIndexer->getLIDs();
   auto LIDs_h = Kokkos::create_mirror_view(LIDs);
   auto colLIDs_h = Kokkos::create_mirror_view(colLIDs);
   Kokkos::deep_copy(LIDs_h, LIDs);
   Kokkos::deep_copy(colLIDs_h, colLIDs);

   std::vector<typename decltype(scatterFields_[0].get_static_view())::HostMirror> scatterFields_h;
   for ( std::size_t i=0; i< scatterFields_.size(); ++i) {
     scatterFields_h.push_back(Kokkos::create_mirror_view(scatterFields_[i].get_static_view()));
     Kokkos::deep_copy(scatterFields_h[i], scatterFields_[i].get_static_view());
   }

   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      gidCount = colGlobalIndexer->getElementBlockGIDCount(blockId);

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
            int row = LIDs_h(cellLocalId, offset);
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
            const ScalarT scatterField = (scatterFields_h[fieldIndex])(worksetCellIndex,basisId);

            if(r!=Teuchos::null) 
              (*r)[row] = scatterField.val();
            if(dirichletCounter_!=Teuchos::null) {
              // std::cout << "Writing " << row << " " << dirichletCounter_->MyLength() << std::endl;
              (*dirichletCounter_)[row] = 1.0; // mark row as dirichlet
            }

            // loop over the sensitivity indices: all DOFs on a cell
            std::vector<double> jacRow(scatterField.size(),0.0);
    
            if(!preserveDiagonal_) {
              int err = Jac->ReplaceMyValues(row, gidCount, scatterField.dx(),
					     &colLIDs_h(cellLocalId,0));
              TEUCHOS_ASSERT(err==0); 
            }
         }
      }
   }
}

// **********************************************************************

#endif
