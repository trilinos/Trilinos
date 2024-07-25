// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCATTER_DIRICHLET_RESIDUAL_BLOCEDEPETRA_IMPL_HPP
#define PANZER_SCATTER_DIRICHLET_RESIDUAL_BLOCEDEPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Teuchos_FancyOStream.hpp"

#include <unordered_map>

// **********************************************************************
// Specialization: Residual
// **********************************************************************


template<typename TRAITS,typename LO,typename GO>
panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                       const std::vector<Teuchos::RCP<const GlobalIndexer> > & cIndexers,
                                       const Teuchos::ParameterList& p,
                                       bool /* useDiscreteAdjoint */)
   : rowIndexers_(rIndexers)
   , colIndexers_(cIndexers)
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
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  indexerIds_.resize(scatterFields_.size());
  subFieldIds_.resize(scatterFields_.size());

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;

    indexerIds_[fd]  = getFieldBlock(fieldName,rowIndexers_);
    subFieldIds_[fd] = rowIndexers_[indexerIds_[fd]]->getFieldNum(fieldName);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].extent(1);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef BlockedEpetraLinearObjContainer BLOC;
   typedef BlockedEpetraLinearObjContainer ELOC;

   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;

   // extract dirichlet counter from container
   Teuchos::RCP<BLOC> blockContainer
         = Teuchos::rcp_dynamic_cast<BLOC>(d.gedc->getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_f(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));

   // extract linear object container
   Teuchos::RCP<const BLOC> blockedContainer = Teuchos::rcp_dynamic_cast<const BLOC>(d.gedc->getDataObject(globalDataKey_));
   Teuchos::RCP<const ELOC> epetraContainer  = Teuchos::rcp_dynamic_cast<const ELOC>(d.gedc->getDataObject(globalDataKey_));

   // if its blocked do this
   if(blockedContainer!=Teuchos::null)
     r_ = (!scatterIC_) ?
            rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_f(),true) :
            rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_x(),true);
   else if(epetraContainer!=Teuchos::null) // if its straight up epetra do this
     r_ = (!scatterIC_) ?
            Thyra::castOrCreateNonconstProductVectorBase<double>(epetraContainer->get_f_th()) :
            Thyra::castOrCreateNonconstProductVectorBase<double>(epetraContainer->get_x_th());

   TEUCHOS_ASSERT(r_!=Teuchos::null);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // loop over each field to be scattered
   Teuchos::ArrayRCP<double> local_r;
   Teuchos::ArrayRCP<double> local_dc;
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      int rowIndexer  = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      rcp_dynamic_cast<SpmdVectorBase<double> >(dirichletCounter_->getNonconstVectorBlock(rowIndexer))
                                                                 ->getNonconstLocalData(ptrFromRef(local_dc));

      // grab local data for inputing
      rcp_dynamic_cast<SpmdVectorBase<double> >(r_->getNonconstVectorBlock(rowIndexer))
                                                  ->getNonconstLocalData(ptrFromRef(local_r));

      auto subRowIndexer = rowIndexers_[rowIndexer];
      auto LIDs = subRowIndexer->getLIDs();
      auto LIDs_h = Kokkos::create_mirror_view(LIDs);
      Kokkos::deep_copy(LIDs_h, LIDs);

      auto field = scatterFields_[fieldIndex].get_view();
      auto field_h = Kokkos::create_mirror_view(field);
      Kokkos::deep_copy(field_h, field);

      BCFieldType::array_type::HostMirror applyBC_h;
      if(checkApplyBC_){
        auto applyBC = applyBC_[fieldIndex].get_static_view();
        applyBC_h = Kokkos::create_mirror_view(applyBC);
        Kokkos::deep_copy(applyBC_h, applyBC);
      }

      // scatter operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

         if (!scatterIC_) {
           // this call "should" get the right ordering according to the Intrepid2 basis
           const std::pair<std::vector<int>,std::vector<int> > & indicePair
             = subRowIndexer->getGIDFieldOffsets_closure(blockId,subFieldNum, side_subcell_dim_, local_side_id_);
           const std::vector<int> & elmtOffset = indicePair.first;
           const std::vector<int> & basisIdMap = indicePair.second;

           // loop over basis functions
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs_h(cellLocalId, offset);
             if(lid<0) // not on this processor!
               continue;

             int basisId = basisIdMap[basis];

             if (checkApplyBC_ and !applyBC_h(worksetCellIndex,basisId))
               continue;

             local_r[lid] = field_h(worksetCellIndex,basisId);

             // record that you set a dirichlet condition
             local_dc[lid] = 1.0;
           }
         } else {
           // this call "should" get the right ordering according to the Intrepid2 basis
           const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

           // loop over basis functions
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs_h(cellLocalId, offset);
             if(lid<0) // not on this processor!
               continue;

             local_r[lid] = field_h(worksetCellIndex,basis);

             // record that you set a dirichlet condition
             local_dc[lid] = 1.0;
           }
         }
      }
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************


template<typename TRAITS,typename LO,typename GO>
panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                       const std::vector<Teuchos::RCP<const GlobalIndexer> > & cIndexers,
                                       const Teuchos::ParameterList& p,
                                       bool /* useDiscreteAdjoint */)
   : rowIndexers_(rIndexers)
   , colIndexers_(cIndexers)
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
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  indexerIds_.resize(scatterFields_.size());
  subFieldIds_.resize(scatterFields_.size());

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;

    indexerIds_[fd]  = getFieldBlock(fieldName,rowIndexers_);
    subFieldIds_[fd] = rowIndexers_[indexerIds_[fd]]->getFieldNum(fieldName);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].extent(1);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef BlockedEpetraLinearObjContainer BLOC;
   typedef BlockedEpetraLinearObjContainer ELOC;

   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;

   // extract dirichlet counter from container
   Teuchos::RCP<BLOC> blockContainer
         = Teuchos::rcp_dynamic_cast<BLOC>(d.gedc->getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_f(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));

   // extract linear object container
   Teuchos::RCP<const BLOC> blockedContainer = Teuchos::rcp_dynamic_cast<const BLOC>(d.gedc->getDataObject(globalDataKey_));
   Teuchos::RCP<const ELOC> epetraContainer  = Teuchos::rcp_dynamic_cast<const ELOC>(d.gedc->getDataObject(globalDataKey_));

   // if its blocked do this
   if(blockedContainer!=Teuchos::null)
     r_ = (!scatterIC_) ?
            rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_f(),true) :
            rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer->get_x(),true);
   else if(epetraContainer!=Teuchos::null) // if its straight up epetra do this
     r_ = (!scatterIC_) ?
            Thyra::castOrCreateNonconstProductVectorBase<double>(epetraContainer->get_f_th()) :
            Thyra::castOrCreateNonconstProductVectorBase<double>(epetraContainer->get_x_th());

   TEUCHOS_ASSERT(r_!=Teuchos::null);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   TEUCHOS_ASSERT(false);

   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::VectorBase;
   using Thyra::SpmdVectorBase;
   using Thyra::ProductVectorBase;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // loop over each field to be scattered
   Teuchos::ArrayRCP<double> local_r;
   Teuchos::ArrayRCP<double> local_dc;
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      int rowIndexer  = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      rcp_dynamic_cast<SpmdVectorBase<double> >(dirichletCounter_->getNonconstVectorBlock(rowIndexer))
                                                                 ->getNonconstLocalData(ptrFromRef(local_dc));

      // grab local data for inputing
      rcp_dynamic_cast<SpmdVectorBase<double> >(r_->getNonconstVectorBlock(rowIndexer))
                                                  ->getNonconstLocalData(ptrFromRef(local_r));

      auto subRowIndexer = rowIndexers_[rowIndexer];

      auto LIDs = subRowIndexer->getLIDs();
      auto LIDs_h = Kokkos::create_mirror_view(LIDs);
      Kokkos::deep_copy(LIDs_h, LIDs);

      auto field = scatterFields_[fieldIndex].get_view();
      auto field_h = Kokkos::create_mirror_view(field);
      Kokkos::deep_copy(field_h, field);

      BCFieldType::array_type::HostMirror applyBC_h;
      if(checkApplyBC_){
        auto applyBC = applyBC_[fieldIndex].get_static_view();
        applyBC_h = Kokkos::create_mirror_view(applyBC);
        Kokkos::deep_copy(applyBC_h, applyBC);
      }

      // scatter operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

         if (!scatterIC_) {
           // this call "should" get the right ordering according to the Intrepid2 basis
           const std::pair<std::vector<int>,std::vector<int> > & indicePair
             = subRowIndexer->getGIDFieldOffsets_closure(blockId,subFieldNum, side_subcell_dim_, local_side_id_);
           const std::vector<int> & elmtOffset = indicePair.first;
           const std::vector<int> & basisIdMap = indicePair.second;

           // loop over basis functions
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs_h(cellLocalId, offset);
             if(lid<0) // not on this processor!
               continue;

             int basisId = basisIdMap[basis];

             if (checkApplyBC_ and !applyBC_h(worksetCellIndex,basisId))
               continue;

             local_r[lid] = field_h(worksetCellIndex,basisId).val();

             // record that you set a dirichlet condition
             local_dc[lid] = 1.0;
           }
         } else {
           // this call "should" get the right ordering according to the Intrepid2 basis
           const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

           // loop over basis functions
           for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
             int offset = elmtOffset[basis];
             int lid = LIDs_h(cellLocalId, offset);
             if(lid<0) // not on this processor!
               continue;

             local_r[lid] = field_h(worksetCellIndex,basis).val();

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

template<typename TRAITS,typename LO,typename GO>
panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                                       const std::vector<Teuchos::RCP<const GlobalIndexer> > & cIndexers,
                                       const Teuchos::ParameterList& p,
                                       bool /* useDiscreteAdjoint */)
   : rowIndexers_(rIndexers)
   , colIndexers_(cIndexers)
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

  if(colIndexers_.size()==0)
    colIndexers_ = rowIndexers_;

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  indexerIds_.resize(scatterFields_.size());
  subFieldIds_.resize(scatterFields_.size());

  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;

    indexerIds_[fd]  = getFieldBlock(fieldName,rowIndexers_);
    subFieldIds_[fd] = rowIndexers_[indexerIds_[fd]]->getFieldNum(fieldName);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].extent(1);
  num_eq = scatterFields_.size();
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef BlockedEpetraLinearObjContainer BLOC;

   using Teuchos::rcp_dynamic_cast;

   // extract dirichlet counter from container
   Teuchos::RCP<const BLOC> blockContainer
         = rcp_dynamic_cast<const BLOC>(d.gedc->getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_f(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));

   // extract linear object container
   blockContainer = rcp_dynamic_cast<const BLOC>(d.gedc->getDataObject(globalDataKey_),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(blockContainer));

   r_ = rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_f());
   Jac_ = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(blockContainer->get_A());
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::ArrayRCP;
   using Teuchos::ptrFromRef;
   using Teuchos::rcp_dynamic_cast;

   using Thyra::SpmdVectorBase;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   int numFieldBlocks = Teuchos::as<int>(colIndexers_.size());

   std::vector<int> blockOffsets;
   computeBlockOffsets(blockId,colIndexers_,blockOffsets);

   std::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsMatrix>,panzer::pair_hash> jacEpetraBlocks;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      int rowIndexer  = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // loop over each field to be scattered
      Teuchos::ArrayRCP<double> local_dc;
      rcp_dynamic_cast<SpmdVectorBase<double> >(dirichletCounter_->getNonconstVectorBlock(rowIndexer))
                                                                 ->getNonconstLocalData(ptrFromRef(local_dc));

      // grab local data for inputing
      Teuchos::ArrayRCP<double> local_r;
      if(r_!=Teuchos::null)
        rcp_dynamic_cast<SpmdVectorBase<double> >(r_->getNonconstVectorBlock(rowIndexer))
                                                    ->getNonconstLocalData(ptrFromRef(local_r));

      auto subRowIndexer = rowIndexers_[rowIndexer];
      auto subIndicePair = subRowIndexer->getGIDFieldOffsets_closure(blockId,subFieldNum, side_subcell_dim_, local_side_id_);
      const std::vector<int> & subElmtOffset = subIndicePair.first;
      const std::vector<int> & subBasisIdMap = subIndicePair.second;

      auto rLIDs = subRowIndexer->getLIDs();
      auto rLIDs_h = Kokkos::create_mirror_view(rLIDs);
      Kokkos::deep_copy(rLIDs_h, rLIDs);

      auto field = scatterFields_[fieldIndex].get_view();
      auto field_h = Kokkos::create_mirror_view(field);
      Kokkos::deep_copy(field_h, field);

      BCFieldType::array_type::HostMirror applyBC_h;
      if(checkApplyBC_){
        auto applyBC = applyBC_[fieldIndex].get_static_view();
        applyBC_h = Kokkos::create_mirror_view(applyBC);
        Kokkos::deep_copy(applyBC_h, applyBC);
      }

      // scatter operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

         // loop over basis functions
         for(std::size_t basis=0;basis<subElmtOffset.size();basis++) {
            int offset = subElmtOffset[basis];
            int lid = rLIDs_h(cellLocalId, offset);
            if(lid<0) // not on this processor
               continue;

            int basisId = subBasisIdMap[basis];

            if (checkApplyBC_ and !applyBC_h(worksetCellIndex,basisId))
              continue;

            // zero out matrix row
            for(int colIndexer=0;colIndexer<numFieldBlocks;colIndexer++) {
               int start = blockOffsets[colIndexer];
               int end = blockOffsets[colIndexer+1];

               if(end-start<=0)
                  continue;

               // check hash table for jacobian sub block
               std::pair<int,int> blockIndex = std::make_pair(rowIndexer,colIndexer);
               Teuchos::RCP<Epetra_CrsMatrix> subJac = jacEpetraBlocks[blockIndex];
               // if you didn't find one before, add it to the hash table
               if(subJac==Teuchos::null) {
                  Teuchos::RCP<Thyra::LinearOpBase<double> > tOp = Jac_->getNonconstBlock(blockIndex.first,blockIndex.second);

                  // block operator is null, don't do anything (it is excluded)
                  if(Teuchos::is_null(tOp))
                     continue;

                  Teuchos::RCP<Epetra_Operator> eOp = Thyra::get_Epetra_Operator(*tOp);
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

            const ScalarT scatterField = field_h(worksetCellIndex,basisId);

            if(r_!=Teuchos::null)
              local_r[lid] = scatterField.val();
            local_dc[lid] = 1.0; // mark row as dirichlet

            // loop over the sensitivity indices: all DOFs on a cell
            std::vector<double> jacRow(scatterField.size(),0.0);

            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);

            for(int colIndexer=0;colIndexer<numFieldBlocks;colIndexer++) {
               int start = blockOffsets[colIndexer];
               int end = blockOffsets[colIndexer+1];

               if(end-start<=0)
                  continue;

               auto subColIndexer = colIndexers_[colIndexer];

               auto cLIDs = subColIndexer->getElementLIDs(cellLocalId);
               auto cLIDs_h = Kokkos::create_mirror_view(cLIDs);
               Kokkos::deep_copy(cLIDs_h, cLIDs);

               TEUCHOS_ASSERT(end-start==Teuchos::as<int>(cLIDs.size()));

               // check hash table for jacobian sub block
               std::pair<int,int> blockIndex = std::make_pair(rowIndexer,colIndexer);
               Teuchos::RCP<Epetra_CrsMatrix> subJac = jacEpetraBlocks[blockIndex];

               // if you didn't find one before, add it to the hash table
               if(subJac==Teuchos::null) {
                  Teuchos::RCP<Thyra::LinearOpBase<double> > tOp = Jac_->getNonconstBlock(blockIndex.first,blockIndex.second);

                  // block operator is null, don't do anything (it is excluded)
                  if(Teuchos::is_null(tOp))
                     continue;

                  Teuchos::RCP<Epetra_Operator> eOp = Thyra::get_Epetra_Operator(*tOp);
                  subJac = rcp_dynamic_cast<Epetra_CrsMatrix>(eOp,true);
                  jacEpetraBlocks[blockIndex] = subJac;
               }

               // Sum Jacobian
               int err = subJac->ReplaceMyValues(lid, end-start, &jacRow[start],&cLIDs_h[0]);
               if(err!=0) {
                 std::stringstream ss;
                 ss << "Failed inserting row: " << " (" << lid << "): ";
                 for(int i=0;i<end-start;i++)
                   ss << cLIDs_h[i] << " ";
                 ss << std::endl;
                 ss << "Into block " << rowIndexer << ", " << colIndexer << std::endl;

                 ss << "scatter field = ";
                 scatterFields_[fieldIndex].print(ss);
                 ss << std::endl;

                 TEUCHOS_TEST_FOR_EXCEPTION(err!=0,std::runtime_error,ss.str());
               }

            }
         }
      }
   }
}

// **********************************************************************

#endif
