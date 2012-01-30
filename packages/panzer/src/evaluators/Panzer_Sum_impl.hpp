#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Sum,p)
{
  std::string sum_name = p.get<std::string>("Sum Name");
  Teuchos::RCP<std::vector<std::string> > value_names = 
    p.get<Teuchos::RCP<std::vector<std::string> > >("Values Names");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  sum = PHX::MDField<ScalarT>(sum_name, data_layout);
  
  this->addEvaluatedField(sum);
 
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<ScalarT>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
 
  std::string n = "Sum Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Sum,worksets,fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < values.size(); ++i)
    this->utils.setFieldData(values[i],fm);

  cell_data_size = sum.size() / sum.fieldTag().dataLayout().dimension(0);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Sum,workset)
{ 
  std::size_t length = workset.num_cells * cell_data_size;

  for (std::size_t i = 0; i < length; ++i) {
    sum[i] = 0.0;
    for (std::size_t j = 0; j < values.size(); ++j)
      sum[i] += (values[j])[i];
  }

}

//**********************************************************************

}
