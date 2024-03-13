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
// Redistribution and use in solution and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of solution code must retain the above copyright
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

#ifndef __Example_SimpleSolution_impl_hpp__
#define __Example_SimpleSolution_impl_hpp__

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace Example {

//**********************************************************************
template <typename EvalT,typename Traits>
SimpleSolution<EvalT,Traits>::SimpleSolution(const std::string & name,
                                             const panzer::IntegrationRule & ir,
                                             const bool curvilinear)
  : curvilinear(curvilinear)
{
  using Teuchos::RCP;

  Teuchos::RCP<PHX::DataLayout> data_layout_scalar = ir.dl_scalar;
  Teuchos::RCP<PHX::DataLayout> data_layout_vector = ir.dl_vector;
  ir_degree = ir.cubature_degree;

  solution = PHX::MDField<ScalarT,Cell,Point>(name, data_layout_scalar);
  solution_grad = PHX::MDField<ScalarT,Cell,Point,Dim>("GRAD_"+name, data_layout_vector);

  this->addEvaluatedField(solution);
  this->addEvaluatedField(solution_grad);
  
  std::string n = "Simple Solution";
  this->setName(n);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SimpleSolution<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
                                                       PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template <typename EvalT,typename Traits>
void SimpleSolution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{ 
  using panzer::index_t;
  auto ip_coordinates = this->wda(workset).int_rules[ir_index]->ip_coordinates.get_static_view();
  auto solution_v = solution.get_static_view();
  auto solution_grad_v = solution_grad.get_static_view();

  const bool curv = this->curvilinear;

  Kokkos::parallel_for (workset.num_cells, KOKKOS_LAMBDA (const index_t cell) {
      for (int point = 0; point < solution_v.extent_int(1); ++point) {

	const double & x = ip_coordinates(cell,point,0);
	const double & y = ip_coordinates(cell,point,1);

  if (!curv) {
	  solution_v(cell,point) = std::sin(2*M_PI*x)*std::sin(2*M_PI*y);
	  solution_grad_v(cell,point,0) = 2.0*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
	  solution_grad_v(cell,point,1) = 2.0*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
  } else {

    const double & theta = std::atan2(y,x);
    const double & r = std::sqrt(x*x + y*y);

    // homogenous solution is T_inner + ( T_inner - T_outer ) / ln (r_inner/r_outer) * ln(r/r_inner)
    
    solution_v(cell,point) = 5.0 + (5.0 - 1.0) / std::log(.5/1.0) * std::log(r/.5);

    const double & fac = (5.0 - 1.0) / std::log(.5/1.0) * 1./r;
    solution_grad_v(cell,point,0) = std::cos(theta) * fac;
    solution_grad_v(cell,point,1) = std::sin(theta) * fac;

  }
      }
    } );

  Kokkos::fence();
}

//**********************************************************************
}

#endif
