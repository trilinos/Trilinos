// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCATTER_RESIDUAL_BLOCKEDEPETRA_IMPL_HPP
#define PANZER_SCATTER_RESIDUAL_BLOCKEDEPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_HashUtils.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Teuchos_FancyOStream.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
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
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
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
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef BlockedEpetraLinearObjContainer BLOC;
   typedef BlockedEpetraLinearObjContainer ELOC;

   // extract linear object container
   Teuchos::RCP<const BLOC> blockedContainer = Teuchos::rcp_dynamic_cast<const BLOC>(d.gedc->getDataObject(globalDataKey_));
   Teuchos::RCP<const ELOC> epetraContainer  = Teuchos::rcp_dynamic_cast<const ELOC>(d.gedc->getDataObject(globalDataKey_));

   // if its blocked do this
   if(blockedContainer!=Teuchos::null)
     r_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_f());
   else if(epetraContainer!=Teuchos::null) // if its straight up epetra do this
     r_ = Thyra::castOrCreateNonconstProductVectorBase<double>(epetraContainer->get_f_th());

   TEUCHOS_ASSERT(r_!=Teuchos::null);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Residual, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
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
   for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      int indexerId   = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // grab local data for inputing
      rcp_dynamic_cast<SpmdVectorBase<double> >(r_->getNonconstVectorBlock(indexerId))->getNonconstLocalData(ptrFromRef(local_r));

      auto subRowIndexer = rowIndexers_[indexerId];
      const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

      auto field = PHX::as_view(scatterFields_[fieldIndex]);
      auto field_h = Kokkos::create_mirror_view(field);
      Kokkos::deep_copy(field_h, field);

      // scatter operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
          std::size_t cellLocalId = localCellIds[worksetCellIndex];

	        auto LIDs = subRowIndexer->getElementLIDs(cellLocalId);
          auto LIDs_h = Kokkos::create_mirror_view(LIDs);
          Kokkos::deep_copy(LIDs_h, LIDs);

         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs_h[offset];
            local_r[lid] += field_h(worksetCellIndex,basis);
         }
      }
   }
}

// **********************************************************************
// Specialization: Tangent
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
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

  this->setName(scatterName+" Scatter Tangent");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
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
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   typedef BlockedEpetraLinearObjContainer BLOC;
   typedef BlockedEpetraLinearObjContainer ELOC;

   // extract linear object container
   Teuchos::RCP<const BLOC> blockedContainer = Teuchos::rcp_dynamic_cast<const BLOC>(d.gedc->getDataObject(globalDataKey_));
   Teuchos::RCP<const ELOC> epetraContainer  = Teuchos::rcp_dynamic_cast<const ELOC>(d.gedc->getDataObject(globalDataKey_));

   // if its blocked do this
   if(blockedContainer!=Teuchos::null)
     r_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_f());
   else if(epetraContainer!=Teuchos::null) // if its straight up epetra do this
     r_ = Thyra::castOrCreateNonconstProductVectorBase<double>(epetraContainer->get_f_th());

   TEUCHOS_ASSERT(r_!=Teuchos::null);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Tangent, TRAITS,LO,GO>::
evaluateFields(typename TRAITS::EvalData workset)
{
   TEUCHOS_ASSERT(false);

   using Teuchos::RCP;
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
   for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      int indexerId   = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // grab local data for inputing
      rcp_dynamic_cast<SpmdVectorBase<double> >(r_->getNonconstVectorBlock(indexerId))->getNonconstLocalData(ptrFromRef(local_r));

      auto subRowIndexer = rowIndexers_[indexerId];
      const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

      // scatter operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

         auto LIDs = subRowIndexer->getElementLIDs(cellLocalId);

         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            local_r[lid] += (scatterFields_[fieldIndex])(worksetCellIndex,basis).val();
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename TRAITS,typename LO,typename GO>
panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer> > & rIndexers,
                              const std::vector<Teuchos::RCP<const GlobalIndexer> > & cIndexers,
                              const Teuchos::ParameterList& p,
                              bool useDiscreteAdjoint)
  : rowIndexers_(rIndexers)
  , colIndexers_(cIndexers)
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
  { TEUCHOS_ASSERT(colIndexers_.size()==0); }

  if(colIndexers_.size()==0)
    colIndexers_ = rowIndexers_;

  this->setName(scatterName+" Scatter Residual BlockedEpetra (Jacobian)");
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
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
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   typedef BlockedEpetraLinearObjContainer BLOC;
   typedef BlockedEpetraLinearObjContainer ELOC;

   // extract linear object container
   RCP<const BLOC> blockedContainer = rcp_dynamic_cast<const BLOC>(d.gedc->getDataObject(globalDataKey_));
   RCP<const ELOC> epetraContainer  = rcp_dynamic_cast<const ELOC>(d.gedc->getDataObject(globalDataKey_));

   // if its blocked do this
   if(blockedContainer!=Teuchos::null) {
     r_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockedContainer->get_f());
     Jac_ = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(blockedContainer->get_A());
   }
   else if(epetraContainer!=Teuchos::null) {
     // if its straight up epetra do this
     if(epetraContainer->get_f_th()!=Teuchos::null)
       r_ = Thyra::castOrCreateNonconstProductVectorBase<double>(epetraContainer->get_f_th());

     // convert it into a blocked operator
     RCP<Thyra::LinearOpBase<double> > J = blockedContainer->get_A_th();
     Jac_ = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(Thyra::nonconstBlock1x1(J));
   }

   TEUCHOS_ASSERT(Jac_!=Teuchos::null);
}

// **********************************************************************
template<typename TRAITS,typename LO,typename GO>
void panzer::ScatterResidual_BlockedEpetra<panzer::Traits::Jacobian, TRAITS,LO,GO>::
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

   std::vector<double> jacRow;

   // for convenience pull out some objects from workset
   std::string blockId = this->wda(workset).block_id;
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   int numFieldBlocks = Teuchos::as<int>(colIndexers_.size());

   std::vector<int> blockOffsets;
   computeBlockOffsets(blockId,colIndexers_,blockOffsets);

   std::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsMatrix>,panzer::pair_hash> jacEpetraBlocks;

   // loop over each field to be scattered
   Teuchos::ArrayRCP<double> local_r;
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      int rowIndexer  = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      // grab local data for inputing
      if(r_!=Teuchos::null)
         rcp_dynamic_cast<SpmdVectorBase<double> >(r_->getNonconstVectorBlock(rowIndexer))->getNonconstLocalData(ptrFromRef(local_r));

      auto subRowIndexer = rowIndexers_[rowIndexer];
      const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

      auto field = scatterFields_[fieldIndex].get_view();
      auto field_h = Kokkos::create_mirror_view(field);
      Kokkos::deep_copy(field_h, field);

      auto rLIDs = subRowIndexer->getLIDs();
      auto rLIDs_h = Kokkos::create_mirror_view(rLIDs);
      Kokkos::deep_copy(rLIDs_h, rLIDs);

      // scatter operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

         // loop over the basis functions (currently they are nodes)
         for(std::size_t rowBasisNum = 0; rowBasisNum < elmtOffset.size(); rowBasisNum++) {
            const ScalarT scatterField = field_h(worksetCellIndex,rowBasisNum);
            int rowOffset = elmtOffset[rowBasisNum];
            int r_lid = rLIDs_h(cellLocalId, rowOffset);

            // Sum residual
            if(local_r!=Teuchos::null)
               local_r[r_lid] += (scatterField.val());

            // loop over the sensitivity indices: all DOFs on a cell
            jacRow.resize(scatterField.size());

            // For Neumann conditions with no dependence on degrees of freedom, there should be no Jacobian contribution
            if(scatterField.size() == 0)
                continue;

            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);

            // scatter the row to each block
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
               std::pair<int,int> blockIndex = useDiscreteAdjoint_ ? std::make_pair(colIndexer,rowIndexer)
                                                                   : std::make_pair(rowIndexer,colIndexer);
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
               if(!useDiscreteAdjoint_) {
                 int err = subJac->SumIntoMyValues(r_lid, end-start, &jacRow[start],&cLIDs_h[0]);
                 if(err!=0) {

                   std::stringstream ss;
                   ss << "Failed inserting row: " << "LID = " << r_lid << "): ";
                   for(int i=start;i<end;i++)
                     ss <<  cLIDs_h[i] << " ";
                   ss << std::endl;
                   ss << "Into block " << rowIndexer << ", " << colIndexer << std::endl;

                   ss << "scatter field = ";
                   scatterFields_[fieldIndex].print(ss);
                   ss << std::endl;

                   TEUCHOS_TEST_FOR_EXCEPTION(err!=0,std::runtime_error,ss.str());
                 }
               }
               else {
                 for(std::size_t c=0;c<cLIDs.size();c++) {
                   int err = subJac->SumIntoMyValues(cLIDs_h[c], 1, &jacRow[start+c],&r_lid);
                   TEUCHOS_ASSERT_EQUALITY(err,0);
                 }
               }
            }
         } // end rowBasisNum
      } // end fieldIndex
   }
}

// **********************************************************************

#endif
