#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_DOFManager.hpp"
#include "Panzer_Basis.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Teuchos_FancyOStream.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename Traits>
panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Residual, Traits>::
ScatterDirichletResidual_Epetra(const Teuchos::ParameterList& p)
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
    p.get< Teuchos::RCP<panzer::Basis> >("Basis")->functional;

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
template<typename Traits> 
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Residual, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  dofManager_ = d.dofManager_;

  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = dofManager_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
}

// **********************************************************************
template<typename Traits>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<typename Traits::GlobalOrdinal> GIDs;
   std::vector<int> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
   Teuchos::RCP<Epetra_Vector> r = workset.residual_vector;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      dofManager_->getElementGIDs(cellLocalId,GIDs); 

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++)
         LIDs[i] = r->Map().LID(GIDs[i]);

      std::vector<bool> is_owned(GIDs.size(), false);
      for(std::size_t i=0;i<GIDs.size();i++)
         is_owned[i] = dofManager_->getMap()->MyGID(GIDs[i]);

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
   
         // this call "should" get the right ordering accordint to the Intrepid basis
         const std::vector<int> & elmtOffset 
               = dofManager_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
    
            if(is_owned[offset])
               (*r)[lid] = (scatterFields_[fieldIndex])(worksetCellIndex,basis);
            else
               (*r)[lid] = 0.0;
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits>
panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Jacobian, Traits>::
ScatterDirichletResidual_Epetra(const Teuchos::ParameterList& p)
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
    p.get< Teuchos::RCP<panzer::Basis> >("Basis")->functional;

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
template<typename Traits> 
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Jacobian, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  dofManager_ = d.dofManager_;

  fieldIds_.resize(scatterFields_.size());
  // load required field numbers for fast use
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    fieldIds_[fd] = dofManager_->getFieldNum(fieldName);

    // fill field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }

  // get the number of nodes (Should be renamed basis)
  num_nodes = scatterFields_[0].dimension(1);
  num_eq = scatterFields_.size();
}

// **********************************************************************
template<typename Traits>
void panzer::ScatterDirichletResidual_Epetra<panzer::Traits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<typename Traits::GlobalOrdinal> GIDs;
   std::vector<int> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
   Teuchos::RCP<Epetra_Vector> r = workset.residual_vector;
   Teuchos::RCP<Epetra_CrsMatrix> Jac = workset.jacobian_matrix;

   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // scatter operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      dofManager_->getElementGIDs(cellLocalId,GIDs); 

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++)
         LIDs[i] = r->Map().LID(GIDs[i]);

      std::vector<bool> is_owned(GIDs.size(), false);
      for(std::size_t i=0;i<GIDs.size();i++)
         is_owned[i] = dofManager_->getMap()->MyGID(GIDs[i]);

      // loop over each field to be scattered
      for(std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
   
         // this call "should" get the right ordering accordint to the Intrepid basis
         const std::vector<int> & elmtOffset 
               = dofManager_->getGIDFieldOffsets_closure(blockId,fieldNum, side_subcell_dim_, local_side_id_);
   
         // loop over basis functions
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            int gid = GIDs[offset];
            const ScalarT & scatterField = (scatterFields_[fieldIndex])(worksetCellIndex,basis);
    
            if(is_owned[offset])
               (*r)[lid] = scatterField.val();
            else
               (*r)[lid] = 0.0;
    
            // loop over the sensitivity indices: all DOFs on a cell
            std::vector<double> jacRow(scatterField.size(),0.0);
    
            if(is_owned[offset]) {
               for(int sensIndex=0;sensIndex<scatterField.size();++sensIndex)
                  jacRow[sensIndex] = scatterField.fastAccessDx(sensIndex);
            }
            TEUCHOS_ASSERT(jacRow.size()==GIDs.size());
    
            int err = Jac->ReplaceGlobalValues(gid, scatterField.size(), &jacRow[0],&GIDs[0]);
            TEUCHOS_ASSERT(err==0); 
         }
      }
   }
}

// **********************************************************************
