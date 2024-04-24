// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
