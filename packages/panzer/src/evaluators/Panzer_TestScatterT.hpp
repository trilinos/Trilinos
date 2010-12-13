#ifndef PANZER_TEST_SCATTER_T_HPP
#define PANZER_TEST_SCATTER_T_HPP

#include "Panzer_DOFManager.hpp"

template <typename EvalT,typename Traits>
int panzer::TestScatter<EvalT, Traits>::offset = 0;

namespace panzer {

PHX_EVALUATOR_CTOR(TestScatter,p)
{
  std::string test_name     = p.get<std::string>("Test Name");
  std::string test_name_res = p.get<std::string>("Test Name Residual");
  Teuchos::RCP<PHX::DataLayout> dl = p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  value = PHX::MDField<ScalarT,Cell,NODE>(p.get<std::string>("Test Name"), dl);
  scatter_value = PHX::MDField<ScalarT,Cell,NODE>(test_name_res, dl);

  this->addDependentField(value);
  this->addEvaluatedField(scatter_value);

  localOffset = offset;

  if(offset==0) offset = 10000;
  else offset *= 10;

  std::string n = scatter_value.fieldTag().name();
  this->setName(n);
}

PHX_POST_REGISTRATION_SETUP(TestScatter,setupData,fm)
{
  this->utils.setFieldData(scatter_value,fm);
  this->utils.setFieldData(value,fm);

  num_nodes = scatter_value.dimension(1);
}

PHX_EVALUATE_FIELDS(TestScatter,workset)
{ 
  for (int i=0; i < scatter_value.size(); ++i)
    scatter_value[i] = 0.0;

  for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
    ScalarT sum = 0.0;
    for (std::size_t node = 0; node < num_nodes; ++node) 
       sum += value(cell,node);
    sum = sum / double(num_nodes);

    for (std::size_t node = 0; node < num_nodes; ++node) {
      //unsigned node_GID = *** need to fix this ***;

      scatter_value(cell,node) = 3.0*sum;
    }
  }
}

}

#endif
