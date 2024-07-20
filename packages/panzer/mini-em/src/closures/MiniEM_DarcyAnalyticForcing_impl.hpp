// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_DARCYANALYTICFORCING_IMPL_HPP
#define MINIEM_DARCYANALYTICFORCING_IMPL_HPP

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

#include "Panzer_Traits.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_MathematicalConstants.hpp"
#include "Kokkos_MathematicalFunctions.hpp"


namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
DarcyAnalyticForcing<EvalT,Traits>::DarcyAnalyticForcing(const std::string & name,
                                           const panzer::IntegrationRule & ir,
                                           const panzer::FieldLayoutLibrary & fl,
                                           const double kappa,
                                           const std::string& basisName)
{
  using Teuchos::RCP;

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis(basisName);

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  source = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);
  this->addEvaluatedField(source);

  std::string n = "Darcy Analytic Forcing";
  this->setName(n);

  kappa_ = kappa;

}

//**********************************************************************
template <typename EvalT,typename Traits>
void DarcyAnalyticForcing<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,
                                                               PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void DarcyAnalyticForcing<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  const double time = workset.time;
  const auto pi = Kokkos::numbers::pi_v<double>;

  const double kappa = kappa_;
  const int dim = ir_dim;

  const auto coords = workset.int_rules[ir_index]->ip_coordinates.get_static_view();
  auto tmp_source = source.get_static_view();

  if (ir_dim == 3) {
    Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,source.extent_int(1)});
    Kokkos::parallel_for("panzer:DarcyAnalyticForcing 3D",policy,KOKKOS_LAMBDA (const int cell,const int point) {
      auto x = coords(cell,point,0);
      auto y = coords(cell,point,1);
      auto z = coords(cell,point,2);
      tmp_source(cell,point) = (pi*Kokkos::cos(pi*time) + pi*pi*kappa*dim*Kokkos::sin(pi*time)) * Kokkos::sin(pi*x) * Kokkos::sin(pi*y) * Kokkos::sin(pi*z);
    });
  } else {
    Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,source.extent_int(1)});
    Kokkos::parallel_for("panzer:DarcyAnalyticForcing 2D",policy,KOKKOS_LAMBDA (const int cell,const int point) {
      auto x = coords(cell,point,0);
      auto y = coords(cell,point,1);
      tmp_source(cell,point) = (pi*Kokkos::cos(pi*time) + pi*pi*kappa*dim*Kokkos::sin(pi*time)) * Kokkos::sin(pi*x) * Kokkos::sin(pi*y);
    });
  }
}

//**********************************************************************
}

#endif
