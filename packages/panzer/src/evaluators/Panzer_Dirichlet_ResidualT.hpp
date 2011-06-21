
#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(DirichletResidual,p)
{
  std::string residual_name = p.get<std::string>("Residual Name");
  std::string dof_name = p.get<std::string>("DOF Name"); 
  std::string value_name = p.get<std::string>("Value Name");

  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout"); 

  residual = PHX::MDField<ScalarT>(residual_name, data_layout);
  dof = PHX::MDField<ScalarT>(dof_name, data_layout);
  value = PHX::MDField<ScalarT>(value_name, data_layout);
  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
  this->addDependentField(value);
 
  std::string n = "Dirichlet Residual Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DirichletResidual,worksets,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(dof,fm);
  this->utils.setFieldData(value,fm);

  cell_data_size = 
    residual.size() / residual.fieldTag().dataLayout().dimension(0);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DirichletResidual,workset)
{ 
  std::size_t length = workset.num_cells * cell_data_size;

  for (std::size_t i = 0; i < length; ++i)
    residual[i] = dof[i] - value[i];
}

//**********************************************************************

}
