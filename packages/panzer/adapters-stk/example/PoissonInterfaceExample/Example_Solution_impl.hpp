// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace Example {

template <typename EvalT,typename Traits>
Solution<EvalT,Traits>::Solution(const std::string& name,
                                 const panzer::IntegrationRule& ir,
                                 const bool in_linear_Robin)
  : linear_Robin(in_linear_Robin)
{
  Teuchos::RCP<PHX::DataLayout> data_layout = ir.dl_scalar;
  ir_degree = ir.cubature_degree;

  solution = PHX::MDField<ScalarT,Cell,Point>(name, data_layout);

  this->addEvaluatedField(solution);
  
  std::string n = "Solution";
  this->setName(n);
}

template <typename EvalT,typename Traits>
void Solution<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd,           
                                                   PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

template <typename EvalT,typename Traits>
void Solution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;
  auto ip_coordinates = workset.int_rules[ir_index]->ip_coordinates.get_static_view();
  auto solution_v = solution.get_static_view();
  auto l_linear_Robin = linear_Robin;
  auto workset_block_0 = (workset.block_id[7] == '0');
  Kokkos::parallel_for ("SimpleSolutin", workset.num_cells, KOKKOS_LAMBDA (const index_t cell) {
    for (int point = 0; point < solution_v.extent_int(1); ++point) {
      const double& x = ip_coordinates(cell,point,0);
      const double& y = ip_coordinates(cell,point,1);

      if (l_linear_Robin) {
        if (ip_coordinates.extent(2) == 2) {
          solution_v(cell,point) = 0.5 - 0.8*x + 0.5*sin(2*M_PI*x)*cos(2*M_PI*y);
        } else {
          const double & z = ip_coordinates(cell,point,2);
          solution_v(cell,point) = 0.5 - 0.8*x + sin(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z)/3.0;
        }
      } else {
        if (workset_block_0)
          solution_v(cell,point) =  0.5 - 0.4*x;
        else
          solution_v(cell,point) = 0.1 - 0.4*x;
      }
    }
  });
  Kokkos::fence();
}

}
