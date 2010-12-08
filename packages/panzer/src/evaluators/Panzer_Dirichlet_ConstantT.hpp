
#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(DirichletConstant,p)
{
  std::string residual_name = p.get<std::string>("Residual Name");
  std::string dof_name = p.get<std::string>("DOF Name");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");  
  value = p.get<double>("Value");

  residual = PHX::MDField<ScalarT>(residual_name, data_layout);
  dof = PHX::MDField<ScalarT>(dof_name, data_layout);
  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
 
  std::string n = "Dirichlet Constant Value Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DirichletConstant,worksets,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(dof,fm);

  cell_data_size = 
    residual.size() / residual.fieldTag().dataLayout().dimension(0);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DirichletConstant,workset)
{ 
  std::size_t length = workset.num_cells * cell_data_size;

  for (std::size_t i = 0; i < length; ++i)
    residual[i] = dof[i] - value;
}

//**********************************************************************

}
