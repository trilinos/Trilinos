// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_TEST_SCATTER_IMPL_HPP
#define PANZER_TEST_SCATTER_IMPL_HPP

namespace panzer {

template <typename EvalT,typename TRAITS>
int panzer::TestScatter<EvalT, TRAITS>::offset = 0;

template<typename EvalT, typename Traits>
TestScatter<EvalT, Traits>::
TestScatter(
  const Teuchos::ParameterList& p)
{
  std::string test_name     = p.get<std::string>("Test Name");
  std::string test_name_res = p.get<std::string>("Test Name Residual");
  Teuchos::RCP<PHX::DataLayout> dl = p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  value = PHX::MDField<const ScalarT,Cell,NODE>(p.get<std::string>("Test Name"), dl);
  scatter_value = PHX::MDField<ScalarT,Cell,NODE>(test_name_res, dl);

  this->addDependentField(value);
  this->addEvaluatedField(scatter_value);

  localOffset = offset;

  if(offset==0) offset = 10000;
  else offset *= 10;

  std::string n = scatter_value.fieldTag().name();
  this->setName(n);
}

template<typename EvalT, typename Traits>
void
TestScatter<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData /* setupData */,
  PHX::FieldManager<Traits>& /* fm */)
{
  num_nodes = scatter_value.extent(1);
}

template<typename EvalT, typename Traits>
void
TestScatter<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
 // for (int i=0; i < scatter_value.size(); ++i)
 //   scatter_value[i] = 0.0;
  Kokkos::deep_copy(scatter_value.get_static_view(), ScalarT(0.0));

  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
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
