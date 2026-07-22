// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_WEAKDIRICHLET_RESIDUAL_IMPL_HPP
#define PANZER_WEAKDIRICHLET_RESIDUAL_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
WeakDirichletResidual<EvalT, Traits>::
WeakDirichletResidual(
  const Teuchos::ParameterList& p)
{
  std::string residual_name = p.get<std::string>("Residual Name");
  std::string flux_name = p.get<std::string>("Flux Name");
  std::string normal_name = p.get<std::string>("Normal Name");
  std::string normal_dot_flux_name = normal_name + " dot " + flux_name;
  std::string dof_name = p.get<std::string>("DOF Name"); 
  std::string value_name = p.get<std::string>("Value Name");
  std::string sigma_name = p.get<std::string>("Sigma Name");
  
  const Teuchos::RCP<const panzer::PureBasis> basis =
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis");

  const Teuchos::RCP<const panzer::IntegrationRule> ir = 
    p.get< Teuchos::RCP<const panzer::IntegrationRule> >("IR");


  residual = PHX::MDField<ScalarT>(residual_name, basis->functional);
  normal_dot_flux_plus_pen = PHX::MDField<ScalarT>(normal_dot_flux_name, ir->dl_scalar);
  flux = PHX::MDField<const ScalarT>(flux_name, ir->dl_vector);
  normal = PHX::MDField<const ScalarT>(normal_name, ir->dl_vector);
  dof = PHX::MDField<const ScalarT>(dof_name, ir->dl_scalar);
  value = PHX::MDField<const ScalarT>(value_name, ir->dl_scalar);
  sigma = PHX::MDField<const ScalarT>(sigma_name, ir->dl_scalar);

  this->addEvaluatedField(residual);
  this->addEvaluatedField(normal_dot_flux_plus_pen);
  this->addDependentField(normal);
  this->addDependentField(flux);
  this->addDependentField(dof);
  this->addDependentField(value);
  this->addDependentField(sigma);
 
  basis_name = panzer::basisIRLayout(basis,*ir)->name();

  std::string n = "Weak Dirichlet Residual Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
WeakDirichletResidual<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  num_ip = flux.extent(1);
  num_dim = flux.extent(2);

  TEUCHOS_ASSERT(flux.extent(1) == normal.extent(1));
  TEUCHOS_ASSERT(flux.extent(2) == normal.extent(2));

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
WeakDirichletResidual<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (std::size_t ip = 0; ip < num_ip; ++ip) {
      normal_dot_flux_plus_pen(cell,ip) = ScalarT(0.0);
      for (std::size_t dim = 0; dim < num_dim; ++dim) {
	normal_dot_flux_plus_pen(cell,ip) += normal(cell,ip,dim) * flux(cell,ip,dim);
      }
      normal_dot_flux_plus_pen(cell,ip) += sigma(cell,ip) * (dof(cell,ip) - value(cell,ip)); 
    }
  }

  if(workset.num_cells>0)
    Intrepid2::FunctionSpaceTools<PHX::exec_space>::
      integrate(residual.get_view(),
                normal_dot_flux_plus_pen.get_view(), 
                (this->wda(workset).bases[basis_index])->weighted_basis_scalar.get_view());
  
}

//**********************************************************************

}

#endif
