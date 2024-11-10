// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_T_SQUARED_MODEL_IMPL_HPP
#define USER_APP_T_SQUARED_MODEL_IMPL_HPP

#include "Panzer_HierarchicParallelism.hpp"

//**********************************************************************
template<typename EvalT, typename Traits>
user_app::TSquaredModel<EvalT, Traits>::
TSquaredModel(const Teuchos::ParameterList& p)
  : k_( p.get<double>("Multiplier") ),
    t_( p.get<std::string>("Source Name"),
        p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
    k_t_squared_( p.get<std::string>("Target Name"),
                  p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addDependentField(t_);
  this->addEvaluatedField(k_t_squared_);
  std::string n = "user_app::TSquared: " + k_t_squared_.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
user_app::TSquaredModel<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  const auto num_cells = workset.num_cells;
  const int num_ip = static_cast<int>(t_.extent(1));
  const double k = k_;
  const auto t = t_;
  const auto k_t_squared = k_t_squared_;
  auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::ExecSpace>(num_cells);
  Kokkos::parallel_for("TSquared evaluator",policy,KOKKOS_LAMBDA(const Kokkos::TeamPolicy<PHX::ExecSpace>::member_type& team){
    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_ip), [&] (const int pt) {
      k_t_squared(cell,pt) = k * t(cell,pt)*t(cell,pt);
    });
  });
}

//**********************************************************************

#endif
