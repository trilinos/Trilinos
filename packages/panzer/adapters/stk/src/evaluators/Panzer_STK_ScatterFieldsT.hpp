#include "Teuchos_TestForException.hpp"

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_Basis.hpp"

#include "Teuchos_FancyOStream.hpp"

namespace panzer_stk {

PHX_EVALUATOR_CTOR(ScatterFields,p) :
   mesh_(p.get<Teuchos::RCP<STK_Interface> >("Mesh"))
{
  using panzer::Cell;
  using panzer::NODE;

  std::string scatterName = p.get<std::string>("Scatter Name");
 
  const std::vector<std::string> & names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::Basis> basis = 
    p.get< Teuchos::RCP<panzer::Basis> >("Basis");

  // build dependent fields
  scatterFields_.resize(names.size());
  stkFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    scatterFields_[fd] = 
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addDependentField(scatterFields_[fd]);
  }

  // setup a dummy field to evaluate
  PHX::Tag<ScalarT> scatterHolder(scatterName,Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
  this->addEvaluatedField(scatterHolder);

  this->setName(scatterName+": STK-Scatter Fields");
}

PHX_POST_REGISTRATION_SETUP(ScatterFields,d,fm)
{
  for (std::size_t fd = 0; fd < scatterFields_.size(); ++fd) {
    std::string fieldName = scatterFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->get_field<VariableField>(fieldName);

    // setup the field data object
    this->utils.setFieldData(scatterFields_[fd],fm);
  }
}

PHX_EVALUATE_FIELDS(ScatterFields,workset)
{
   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
   std::string blockId = workset.block_id;

   for(std::size_t fieldIndex=0; fieldIndex<scatterFields_.size();fieldIndex++)
      mesh_->setSolutionFieldData(scatterFields_[fieldIndex].fieldTag().name(),blockId,
                                  localCellIds,scatterFields_[fieldIndex]);
}

} // end panzer_stk
