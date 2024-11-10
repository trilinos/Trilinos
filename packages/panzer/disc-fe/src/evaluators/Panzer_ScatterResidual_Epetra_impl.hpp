// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCATTER_RESIDUAL_EPETRA_IMPL_HPP
#define PANZER_SCATTER_RESIDUAL_EPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_PtrFromStlVector.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_ParameterList_GlobalEvaluationData.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Teuchos_FancyOStream.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::ScatterResidual_Epetra<panzer::Traits::Residual, TRAITS,LO,GO>::
ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                       const Teuchos::RCP<const panzer::GlobalIndexer> & /* cIndexer */,
                       const Teuchos::ParameterList& p,
                       bool useDiscreteAdjoint)
  : globalIndexer_(indexer)
  , globalDataKey_("Residual Scatter Container")
  , useDiscreteAdjoint_(useDiscreteAdjoint)
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
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

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
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Residual, TRAITS,LO,GO>::
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
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Residual, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // extract linear object container
  epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject(globalDataKey_));

  if(epetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(loc);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Residual, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<Epetra_Vector> r = epetraContainer_->get_f();

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   auto LIDs = globalIndexer_->getLIDs();
   auto LIDs_h = Kokkos::create_mirror_view(LIDs);
   Kokkos::deep_copy(LIDs_h, LIDs);
   // loop over each field to be scattered
   for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
     int fieldNum = fieldIds_[fieldIndex];
     auto field = PHX::as_view(scatterFields_[fieldIndex]);
     auto field_h = Kokkos::create_mirror_view(field);
     Kokkos::deep_copy(field_h, field);

     const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
     for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
       std::size_t cellLocalId = localCellIds[worksetCellIndex];

         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs_h(cellLocalId, offset);
            (*r)[lid] += field_h(worksetCellIndex,basis);
         }
      }
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::ScatterResidual_Epetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
ScatterResidual_Epetra(const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
                       const Teuchos::RCP<const panzer::GlobalIndexer> & /* cIndexer */,
                       const Teuchos::ParameterList& p,
                       bool /* useDiscreteAdjoint */)
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
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  this->setName(scatterName+" Scatter Tangent");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
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
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // this is the list of parameters and their names that this scatter has to account for
  std::vector<std::string> activeParameters =
    rcp_dynamic_cast<ParameterList_GlobalEvaluationData>(d.gedc->getDataObject("PARAMETER_NAMES"))->getActiveParameters();

  dfdp_vectors_.clear();
  for(std::size_t i=0;i<activeParameters.size();i++) {
    RCP<Epetra_Vector> vec = rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject(activeParameters[i]),true)->get_f();
    dfdp_vectors_.push_back(vec);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // loop over each field to be scattered
   auto LIDs = globalIndexer_->getLIDs();
   auto LIDs_h = Kokkos::create_mirror_view(LIDs);
   Kokkos::deep_copy(LIDs_h, LIDs);
   for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
     int fieldNum = fieldIds_[fieldIndex];
     const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
     auto scatterField_h = Kokkos::create_mirror_view(scatterFields_[fieldIndex].get_static_view());
     Kokkos::deep_copy(scatterField_h, scatterFields_[fieldIndex].get_static_view());
     // scatter operation for each cell in workset
     for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
       std::size_t cellLocalId = localCellIds[worksetCellIndex];

       auto LIDs = globalIndexer_->getElementLIDs(cellLocalId);
       // loop over basis functions
       for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
	 int offset = elmtOffset[basis];
	 int lid = LIDs_h(cellLocalId, offset);
	 
	 // scatter the sensitivity vectors
	 ScalarT value = scatterField_h(worksetCellIndex,basis);
	 for(int d=0;d<value.size();d++)
	   (*dfdp_vectors_[d])[lid] += value.fastAccessDx(d);
       }
     }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::ScatterResidual_Epetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
ScatterResidual_Epetra(const Teuchos::RCP<const GlobalIndexer> & indexer,
                       const Teuchos::RCP<const panzer::GlobalIndexer> & cIndexer,
                       const Teuchos::ParameterList& p,
                       bool useDiscreteAdjoint)
   : globalIndexer_(indexer)
   , colGlobalIndexer_(cIndexer)
   , globalDataKey_("Residual Scatter Container")
   , useDiscreteAdjoint_(useDiscreteAdjoint)
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
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");
  if (p.isType<bool>("Use Discrete Adjoint"))
     useDiscreteAdjoint = p.get<bool>("Use Discrete Adjoint");

  // discrete adjoint does not work with non-square matrices
  if(useDiscreteAdjoint)
  { TEUCHOS_ASSERT(colGlobalIndexer_==globalIndexer_); }

  this->setName(scatterName+" Scatter Residual Epetra (Jacobian)");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
postRegistrationSetup(typename TRAITS::SetupData /* d */,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  // globalIndexer_ = d.globalIndexer_;

  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // extract linear object container
  epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(d.gedc->getDataObject(globalDataKey_));

  if(epetraContainer_==Teuchos::null) {
    // extract linear object container
    Teuchos::RCP<LinearObjContainer> loc = Teuchos::rcp_dynamic_cast<LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true)->getGhostedLOC();
    epetraContainer_ = Teuchos::rcp_dynamic_cast<EpetraLinearObjContainer>(loc);
  }
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_Epetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   std::vector<double> jacRow;

   bool useColumnIndexer = colGlobalIndexer_!=Teuchos::null;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   Teuchos::RCP<Epetra_Vector> r = epetraContainer_->get_f();
   Teuchos::RCP<Epetra_CrsMatrix> Jac = epetraContainer_->get_A();

   const Teuchos::RCP<const panzer::GlobalIndexer>&
     colGlobalIndexer = useColumnIndexer ? colGlobalIndexer_ : globalIndexer_;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets" may be expensive.  However the
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

      std::vector<int> cLIDs;
      for (int i(0); i < static_cast<int>(colLIDs_h.extent(1)); ++i)
        cLIDs.push_back(colLIDs_h(cellLocalId, i));
      if (Teuchos::nonnull(workset.other)) {
        const std::size_t other_cellLocalId = workset.other->cell_local_ids[worksetCellIndex];
        for (int i(0); i < static_cast<int>(colLIDs_h.extent(1)); ++i)
          cLIDs.push_back(colLIDs_h(other_cellLocalId, i));
      }

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

         // loop over the basis functions (currently they are nodes)
         for(std::size_t rowBasisNum = 0; rowBasisNum < elmtOffset.size(); rowBasisNum++) {
            const ScalarT scatterField = (scatterFields_h[fieldIndex])(worksetCellIndex,rowBasisNum);
            int rowOffset = elmtOffset[rowBasisNum];
            int row = LIDs_h(cellLocalId,rowOffset);

            // Sum residual
            if(r!=Teuchos::null)
                   r->SumIntoMyValue(row,0,scatterField.val());

            // Sum Jacobian
            if(useDiscreteAdjoint_) {
               // loop over the sensitivity indices: all DOFs on a cell
               jacRow.resize(scatterField.size());

               for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
                  jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);

               for(std::size_t c=0;c<cLIDs.size();c++) {
                  int err = Jac->SumIntoMyValues(cLIDs[c], 1, &jacRow[c],&row);
                  TEUCHOS_ASSERT_EQUALITY(err,0);
               }
            }
            else {
               int err = Jac->SumIntoMyValues(
                 row,
                 std::min(cLIDs.size(), static_cast<size_t>(scatterField.size())),
                 scatterField.dx(),
                 panzer::ptrFromStlVector(cLIDs));

               TEUCHOS_ASSERT_EQUALITY(err,0);
            }
         } // end rowBasisNum
      } // end fieldIndex
   }
}

// **********************************************************************

#endif
