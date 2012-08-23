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

#ifndef PANZER_TEST_SCATTER_IMPL_HPP
#define PANZER_TEST_SCATTER_IMPL_HPP

#include "Panzer_DOFManager.hpp"

template <typename EvalT,typename Traits>
int panzer::TestScatter<EvalT, Traits>::offset = 0;

namespace panzer {

PHX_EVALUATOR_CTOR(TestScatter,p)
{
  std::string test_name     = p.get<std::string>("Test Name");
  std::string test_name_res = p.get<std::string>("Test Name Residual");
  Teuchos::RCP<PHX::DataLayout> dl = p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  value = PHX::MDField<ScalarT,Cell,NODE>(p.get<std::string>("Test Name"), dl);
  scatter_value = PHX::MDField<ScalarT,Cell,NODE>(test_name_res, dl);

  this->addDependentField(value);
  this->addEvaluatedField(scatter_value);

  localOffset = offset;

  if(offset==0) offset = 10000;
  else offset *= 10;

  std::string n = scatter_value.fieldTag().name();
  this->setName(n);
}

PHX_POST_REGISTRATION_SETUP(TestScatter,setupData,fm)
{
  this->utils.setFieldData(scatter_value,fm);
  this->utils.setFieldData(value,fm);

  num_nodes = scatter_value.dimension(1);
}

PHX_EVALUATE_FIELDS(TestScatter,workset)
{ 
  for (int i=0; i < scatter_value.size(); ++i)
    scatter_value[i] = 0.0;

  for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
    ScalarT sum = 0.0;
    for (std::size_t node = 0; node < num_nodes; ++node) 
       sum += value(cell,node);
    sum = sum / double(num_nodes);

    for (std::size_t node = 0; node < num_nodes; ++node) {
      //unsigned node_GID = *** need to fix this ***;

      scatter_value(cell,node) = 3.0*sum;
    }
  }
}

}

#endif
