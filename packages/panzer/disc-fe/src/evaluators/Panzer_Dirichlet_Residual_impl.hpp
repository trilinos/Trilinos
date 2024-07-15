// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DIRICHLET_RESIDUAL_IMPL_HPP
#define PANZER_DIRICHLET_RESIDUAL_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
DirichletResidual<EvalT, Traits>::
DirichletResidual(
  const Teuchos::ParameterList& p)
{
  std::string residual_name = p.get<std::string>("Residual Name");
  std::string dof_name = p.get<std::string>("DOF Name"); 
  std::string value_name = p.get<std::string>("Value Name");

  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout"); 

  residual = PHX::MDField<ScalarT>(residual_name, data_layout);
  dof = PHX::MDField<const ScalarT>(dof_name, data_layout);
  value = PHX::MDField<const ScalarT>(value_name, data_layout);
  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
  this->addDependentField(value);
 
  std::string n = "Dirichlet Residual Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DirichletResidual<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  /* fm */)
{
  cell_data_size = residual.fieldTag().dataLayout().extent(1);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DirichletResidual<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  auto residual_v = residual.get_static_view();
  auto dof_v = dof.get_static_view();
  auto value_v = value.get_static_view();
  auto local_cell_data_size = cell_data_size;
  Kokkos::parallel_for (workset.num_cells, KOKKOS_LAMBDA (index_t i) {
    for (std::size_t j = 0; j < local_cell_data_size; ++j)
      residual_v(i,j)=dof_v(i,j)-value_v(i,j);
  });
}

//**********************************************************************

}

#endif
