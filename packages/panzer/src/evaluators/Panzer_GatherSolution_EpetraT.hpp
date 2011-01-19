#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_Basis.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename Traits>
Epetra_Vector &
panzer::GatherSolution_Epetra<panzer::Traits::Residual, Traits>::
getEpetraVector(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & in_v) const
{
   Teuchos::RCP<const Epetra_Map> eMap = Teuchos::get_extra_data<Teuchos::RCP<const Epetra_Map> >(in_v,"epetra_map");
   return *Thyra::get_Epetra_MultiVector(*eMap,in_v)->operator()(0);
}

template<typename Traits>
panzer::GatherSolution_Epetra<panzer::Traits::Residual, Traits>::
GatherSolution_Epetra(const Teuchos::ParameterList& p)
{ 
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  Teuchos::RCP<panzer::Basis> basis = 
    p.get< Teuchos::RCP<panzer::Basis> >("Basis");

  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] = 
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  this->setName("Gather Solution");
}

// **********************************************************************
template<typename Traits> 
void panzer::GatherSolution_Epetra<panzer::Traits::Residual, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  globalIndexer_ = d.globalIndexer_;

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    std::string fieldName = gatherFields_[fd].fieldTag().name();
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }
}

// **********************************************************************
template<typename Traits>
void panzer::GatherSolution_Epetra<panzer::Traits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<typename Traits::GlobalOrdinal> GIDs;
   std::vector<int> LIDs;
 
   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
   Teuchos::RCP<Epetra_Vector> x; 
   if(workset.th_solution_vector!=Teuchos::null) {
      Epetra_Vector & temp_X = getEpetraVector(workset.th_solution_vector);
      x = Teuchos::rcpFromRef(temp_X);
   }
   else
      x = workset.solution_vector;
 
   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!
 
   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];
 
      globalIndexer_->getElementGIDs(cellLocalId,GIDs); 
 
      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++)
         LIDs[i] = x->Map().LID(GIDs[i]);
 
      // loop over the fields to be gathered
      for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);
 
         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];
            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = (*x)[lid];
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits>
Epetra_Vector &
panzer::GatherSolution_Epetra<panzer::Traits::Jacobian, Traits>::
getEpetraVector(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & in_v) const
{
   Teuchos::RCP<const Epetra_Map> eMap = Teuchos::get_extra_data<Teuchos::RCP<const Epetra_Map> >(in_v,"epetra_map");
   return *Thyra::get_Epetra_MultiVector(*eMap,in_v)->operator()(0);
}

template<typename Traits>
panzer::GatherSolution_Epetra<panzer::Traits::Jacobian, Traits>::
GatherSolution_Epetra(const Teuchos::ParameterList& p)
{ 
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("DOF Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<panzer::Basis> >("Basis")->functional;

  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    PHX::MDField<ScalarT,Cell,NODE> f(names[fd],dl);
    gatherFields_[fd] = f;
    this->addEvaluatedField(gatherFields_[fd]);
  }

  this->setName("Gather Solution");
}

// **********************************************************************
template<typename Traits> 
void panzer::GatherSolution_Epetra<panzer::Traits::Jacobian, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  globalIndexer_ = d.globalIndexer_;

  fieldIds_.resize(gatherFields_.size());

  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    std::string fieldName = gatherFields_[fd].fieldTag().name();
    fieldIds_[fd] = globalIndexer_->getFieldNum(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }
}

// **********************************************************************
template<typename Traits>
void panzer::GatherSolution_Epetra<panzer::Traits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
   std::vector<typename Traits::GlobalOrdinal> GIDs;
   std::vector<int> LIDs;

   // for convenience pull out some objects from workset
   std::string blockId = workset.block_id;
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
   Teuchos::RCP<Epetra_Vector> x; 
   if(workset.th_solution_vector!=Teuchos::null)
      x = Teuchos::rcpFromRef(getEpetraVector(workset.th_solution_vector));
   else
      x = workset.solution_vector;
 
   // NOTE: A reordering of these loops will likely improve performance
   //       The "getGIDFieldOffsets may be expensive.  However the
   //       "getElementGIDs" can be cheaper. However the lookup for LIDs
   //       may be more expensive!

   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];

      globalIndexer_->getElementGIDs(cellLocalId,GIDs); 

      // caculate the local IDs for this element
      LIDs.resize(GIDs.size());
      for(std::size_t i=0;i<GIDs.size();i++)
        LIDs[i] = x->Map().LID(GIDs[i]);

      // loop over the fields to be gathered
      for(std::size_t fieldIndex=0;
          fieldIndex<gatherFields_.size();fieldIndex++) {
         int fieldNum = fieldIds_[fieldIndex];
         const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
            int offset = elmtOffset[basis];
            int lid = LIDs[offset];

            // set the value and seed the FAD object
            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = Sacado::Fad::DFad<double>(GIDs.size(), (*x)[lid]);
            (gatherFields_[fieldIndex])(worksetCellIndex,basis).fastAccessDx(offset) = 1.0;
         }
      }
   }
}

// **********************************************************************
