// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Example_SineSource_impl_hpp__
#define __Example_SineSource_impl_hpp__

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_HierarchicParallelism.hpp"

namespace Example {

//**********************************************************************
template <typename EvalT,typename Traits>
SineSource<EvalT,Traits>::SineSource(const std::string & name,
                                         const panzer::IntegrationRule & ir)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
  ir_degree = ir.cubature_degree;

  source = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);

  this->addEvaluatedField(source);

  std::string n = "Sine Source";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SineSource<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,
                                                       PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename Scalar>
class SineSourceFunctor {
  PHX::View<Scalar**> source;           // source at ip
  PHX::View<const double***> ip_coords; // coordinates
public:
  SineSourceFunctor(const PHX::View<Scalar**>& in_source,
		    const PHX::View<const double***>& in_ip_coords)
    : source(in_source), ip_coords(in_ip_coords) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
  {
    const int cell = team.league_rank();
    const int num_points = static_cast<int>(source.extent(1));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int& pt) {
      const double & x = ip_coords(cell,pt,0);
      const double & y = ip_coords(cell,pt,1);
      const double & z = ip_coords(cell,pt,2);
      source(cell,pt) = -12.0*M_PI*M_PI*std::sin(2.0*M_PI*x)*std::sin(2*M_PI*y)*std::sin(2.0*M_PI*z);
    });
  }
};

//**********************************************************************
template <typename EvalT,typename Traits>
void SineSource<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  Example::SineSourceFunctor<ScalarT> ssf(source.get_static_view(),
					  this->wda(workset).int_rules[ir_index]->ip_coordinates.get_static_view());

  auto policy = panzer::HP::inst().teamPolicy<ScalarT>(workset.num_cells);

  Kokkos::parallel_for("MixedPoisson SineSource",policy,ssf);
}

//**********************************************************************
}

#endif
