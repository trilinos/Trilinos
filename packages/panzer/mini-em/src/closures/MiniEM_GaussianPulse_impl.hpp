// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MINIEM_GAUSSIAN_PULSE_IMPL_HPP
#define MINIEM_GAUSSIAN_PULSE_IMPL_HPP

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

namespace mini_em {

//**********************************************************************
template <typename EvalT,typename Traits>
GaussianPulse<EvalT,Traits>::GaussianPulse(const std::string & name,
                                           const panzer::IntegrationRule & ir,
                                           const panzer::FieldLayoutLibrary & fl,
                                           const double & dt)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_vector;
  ir_degree = ir.cubature_degree;
  ir_dim = ir.spatial_dimension;

  current = PHX::MDField<ScalarT,Cell,Point,Dim>(name, data_layout);
  this->addEvaluatedField(current);

  alpha = 1.0;
  beta  = 5.0*dt;

  Teuchos::RCP<const panzer::PureBasis> basis = fl.lookupBasis("E_edge");

  std::string n = "Gaussian Pulse";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void GaussianPulse<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,
                                                        PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void GaussianPulse<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  double time = workset.time;

  const double factor = std::exp(-(time-2.0*beta)*(time-2.0*beta)/beta/beta);
  const double scale = 1.0/alpha/alpha;
  const auto coords = workset.int_rules[ir_index]->ip_coordinates.get_static_view();
  auto tmp_current = current;
  if (ir_dim == 3) {
    Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,current.extent_int(1)});
    Kokkos::parallel_for("panzer:GaussianPulse 3D",policy,KOKKOS_LAMBDA (const int cell,const int point) {
      auto x = coords(cell,point,0);
      auto y = coords(cell,point,1);
      auto z = coords(cell,point,2);
      auto r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5);
      tmp_current(cell,point,0) = 0.0;
      tmp_current(cell,point,1) = 0.0;
      tmp_current(cell,point,2) = std::exp(-r2*scale)*factor;
    });
  } else {
    Kokkos::MDRangePolicy<PHX::exec_space,Kokkos::Rank<2>> policy({0,0},{workset.num_cells,current.extent_int(1)});
    Kokkos::parallel_for("panzer:GaussianPulse 2D",policy,KOKKOS_LAMBDA (const int cell,const int point) {
      auto x = coords(cell,point,0);
      auto y = coords(cell,point,1);
      auto r2 = (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5);
      tmp_current(cell,point,0) = 0.0;
      tmp_current(cell,point,1) = std::exp(-r2*scale)*factor;
    });
  }
}

//**********************************************************************
}

#endif
