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

#include <cmath>

#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace Example {

template <typename EvalT,typename Traits>
Solution<EvalT,Traits>::Solution(const std::string& name,
                                 const panzer::IntegrationRule& ir,
                                 const bool linear_Robin)
  : linear_Robin(linear_Robin)
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
                                                   PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(solution,fm);

  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0], this->wda);
}

template <typename EvalT,typename Traits>
void Solution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int point = 0; point < solution.extent_int(1); ++point) {
      const double& x = this->wda(workset).int_rules[ir_index]->ip_coordinates(cell,point,0);
      const double& y = this->wda(workset).int_rules[ir_index]->ip_coordinates(cell,point,1);

      if (linear_Robin) {
        if (this->wda(workset).int_rules[ir_index]->ip_coordinates.extent(2) == 2) {
          solution(cell,point) = 0.5 - 0.8*x + 0.5*sin(2*M_PI*x)*cos(2*M_PI*y);
        } else {
          const double & z = this->wda(workset).int_rules[ir_index]->ip_coordinates(cell,point,2);
          solution(cell,point) = 0.5 - 0.8*x + sin(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z)/3.0;
        }
      } else {
        if (workset.block_id[7] == '0')
          solution(cell,point) =  0.5 - 0.4*x;
        else
          solution(cell,point) = 0.1 - 0.4*x;
      }
    }
  }
}

}
