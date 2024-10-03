// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
#ifndef __Panzer_ScatterResidual_BlockedEpetra_Hessian_impl_hpp__
#define __Panzer_ScatterResidual_BlockedEpetra_Hessian_impl_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// Epetra scatter residual file

namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
ScatterResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
ScatterResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > & rIndexers,
                              const std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > & cIndexers,
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

  this->setName(scatterName+" Scatter Residual BlockedEpetra (Hessian)");
}
  
template<typename TRAITS,typename LO,typename GO>
void
ScatterResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
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

template<typename TRAITS,typename LO,typename GO>
void
ScatterResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
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
     Jac_ = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(blockedContainer->get_A());
   }
   else if(epetraContainer!=Teuchos::null) {
     // convert it into a blocked operator
     RCP<Thyra::LinearOpBase<double> > J = blockedContainer->get_A_th();
     Jac_ = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(Thyra::nonconstBlock1x1(J));
   }

   TEUCHOS_ASSERT(Jac_!=Teuchos::null);
}
  
template<typename TRAITS,typename LO,typename GO>
void
ScatterResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
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
   for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
      int rowIndexer  = indexerIds_[fieldIndex];
      int subFieldNum = subFieldIds_[fieldIndex];

      auto subRowIndexer = rowIndexers_[rowIndexer];
      const std::vector<int> & elmtOffset = subRowIndexer->getGIDFieldOffsets(blockId,subFieldNum);

      // scatter operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

	 auto rLIDs = subRowIndexer->getElementLIDs(cellLocalId); 

         // loop over the basis functions (currently they are nodes)
         for(std::size_t rowBasisNum = 0; rowBasisNum < elmtOffset.size(); rowBasisNum++) {
            const ScalarT scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,rowBasisNum);
            int rowOffset = elmtOffset[rowBasisNum];
            int r_lid = rLIDs[rowOffset];

            // loop over the sensitivity indices: all DOFs on a cell
            jacRow.resize(scatterField.size());
  
            // For Neumann conditions with no dependence on degrees of freedom, there should be no Jacobian contribution
            if(scatterField.size() == 0)
                continue;
 
            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex).fastAccessDx(0);
    
            // scatter the row to each block
            for(int colIndexer=0;colIndexer<numFieldBlocks;colIndexer++) {
               int start = blockOffsets[colIndexer];
               int end = blockOffsets[colIndexer+1];

               if(end-start<=0) 
                  continue;

               auto subColIndexer = colIndexers_[colIndexer];
	       auto cLIDs = subColIndexer->getElementLIDs(cellLocalId); 

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
               {
                 int err = subJac->SumIntoMyValues(r_lid, end-start, &jacRow[start],&cLIDs[0]);
                 if(err!=0) {
  
                   std::stringstream ss;
                   ss << "Failed inserting row: " << "LID = " << r_lid << ": ";
                   for(int i=0;i<end-start;i++)
                     ss <<  cLIDs[i] << " ";
                   ss << std::endl;
                   ss << "Into block " << rowIndexer << ", " << colIndexer << std::endl;
  
                   ss << "scatter field = ";
                   scatterFields_[fieldIndex].print(ss);
                   ss << std::endl;

                   ss << "values = ";
                   for(int i=start;i<end;i++)
                     ss <<  jacRow[i] << " ";
                   ss << std::endl;

                   std::cout << ss.str() << std::endl;
                 
                   TEUCHOS_TEST_FOR_EXCEPTION(err!=0,std::runtime_error,ss.str());
                 }
               }
            }
         } // end rowBasisNum
      } // end fieldIndex
   }
}

}

// **************************************************************
#endif

#endif
