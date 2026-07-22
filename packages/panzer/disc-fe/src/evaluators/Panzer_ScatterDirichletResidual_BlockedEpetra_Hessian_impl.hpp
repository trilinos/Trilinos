// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ScatterDirichletResidual_BlockedEpetra_Hessian_impl_hpp__
#define __Panzer_ScatterDirichletResidual_BlockedEpetra_Hessian_impl_hpp__

// only do this if required by the user
#ifdef Panzer_BUILD_HESSIAN_SUPPORT

// the includes for this file come in as a result of the includes in the main 
// blocked Epetra scatter dirichlet residual file


namespace panzer {

// **************************************************************
// Hessian Specialization
// **************************************************************
template<typename TRAITS,typename LO,typename GO>
ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
ScatterDirichletResidual_BlockedEpetra(const std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > & rIndexers,
                                       const std::vector<Teuchos::RCP<const GlobalIndexer<LO,int> > > & cIndexers,
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

  this->setName(scatterName+" Scatter Dirichlet Residual : BlockedEpetra (Hessian)");
}
  
template<typename TRAITS,typename LO,typename GO>
void
ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
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

template<typename TRAITS,typename LO,typename GO>
void
ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
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

   Jac_ = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(blockContainer->get_A());
}
  
template<typename TRAITS,typename LO,typename GO>
void
ScatterDirichletResidual_BlockedEpetra<panzer::Traits::Hessian,TRAITS,LO,GO>::
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

      auto subRowIndexer = rowIndexers_[rowIndexer];
      auto subIndicePair = subRowIndexer->getGIDFieldOffsets_closure(blockId,subFieldNum, side_subcell_dim_, local_side_id_);
      const std::vector<int> & subElmtOffset = subIndicePair.first;
      const std::vector<int> & subBasisIdMap = subIndicePair.second;

      // scatter operation for each cell in workset
      for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
         std::size_t cellLocalId = localCellIds[worksetCellIndex];

	 auto rLIDs = subRowIndexer->getElementLIDs(cellLocalId); 

         // loop over basis functions
         for(std::size_t basis=0;basis<subElmtOffset.size();basis++) {
            int offset = subElmtOffset[basis];
            int lid = rLIDs[offset];
            if(lid<0) // not on this processor
               continue;

            int basisId = subBasisIdMap[basis];

            if (checkApplyBC_)
              if (!applyBC_[fieldIndex](worksetCellIndex,basisId))
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
 
            const ScalarT scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,basisId);

            local_dc[lid] = 1.0; // mark row as dirichlet
    
            // loop over the sensitivity indices: all DOFs on a cell
            std::vector<double> jacRow(scatterField.size(),0.0);
    
            for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
               jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex).fastAccessDx(0);
    
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
               int err = subJac->ReplaceMyValues(lid, end-start, &jacRow[start],&cLIDs[0]);
               if(err!=0) {
                 std::stringstream ss;
                 ss << "Failed inserting row: " << " (" << lid << "): ";
                 for(int i=0;i<end-start;i++)
                   ss << cLIDs[i] << " ";
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

}

// **************************************************************
#endif

#endif
