#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_Basis.hpp"

#include "Teuchos_FancyOStream.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename Traits>
panzer_stk::GatherFields<panzer::Traits::Residual, Traits>::
  GatherFields(const Teuchos::RCP<const STK_Interface> & mesh,const Teuchos::ParameterList& p)
{ 
  using panzer::Cell;
  using panzer::NODE;
 
  mesh_ = mesh;

  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::Basis> basis = 
    p.get< Teuchos::RCP<panzer::Basis> >("Basis");

  gatherFields_.resize(names.size());
  stkFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] = 
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  this->setName("STK-Gather Fields");
}

// **********************************************************************
template<typename Traits>
void panzer_stk::GatherFields<panzer::Traits::Residual, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    std::string fieldName = gatherFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->get_field<VariableField>(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }
}

// **********************************************************************
template<typename Traits>
void panzer_stk::GatherFields<panzer::Traits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
   const std::vector<stk::mesh::Entity*> & localElements = *mesh_->getElementsOrderedByLID();
 
   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
 
   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];
      stk::mesh::PairIterRelation relations = localElements[cellLocalId]->relations(mesh_->getNodeRank());
 
      // loop over the fields to be gathered
      for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {
         VariableField * field = stkFields_[fieldIndex];

         std::size_t basisCnt = gatherFields_[fieldIndex].dimension(1);

         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<basisCnt;basis++) {
            stk::mesh::Entity * node = relations[basis].entity();
            stk::mesh::EntityArray<VariableField> fieldData(*field,*node);
            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = fieldData(); // from STK
         }
      }
   }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits>
panzer_stk::GatherFields<panzer::Traits::Jacobian, Traits>::
  GatherFields(const Teuchos::RCP<const STK_Interface> & mesh,const Teuchos::ParameterList& p)
{ 
  using panzer::Cell;
  using panzer::NODE;

  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::Basis> basis = 
    p.get< Teuchos::RCP<panzer::Basis> >("Basis");

  gatherFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] = 
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  this->setName("Gather STK Fields");
}

// **********************************************************************
template<typename Traits> 
void panzer_stk::GatherFields<panzer::Traits::Jacobian, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    // get field ID from DOF manager
    std::string fieldName = gatherFields_[fd].fieldTag().name();

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }
}

// **********************************************************************
template<typename Traits>
void panzer_stk::GatherFields<panzer::Traits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
}

// **********************************************************************
