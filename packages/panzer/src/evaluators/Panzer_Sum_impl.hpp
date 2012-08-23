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

#ifndef PANZER_SUM_IMPL_HPP
#define PANZER_SUM_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Sum,p)
{
  std::string sum_name = p.get<std::string>("Sum Name");
  Teuchos::RCP<std::vector<std::string> > value_names = 
    p.get<Teuchos::RCP<std::vector<std::string> > >("Values Names");
  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");
  
  sum = PHX::MDField<ScalarT>(sum_name, data_layout);
  
  this->addEvaluatedField(sum);
 
  values.resize(value_names->size());
  for (std::size_t i=0; i < value_names->size(); ++i) {
    values[i] = PHX::MDField<ScalarT>( (*value_names)[i], data_layout);
    this->addDependentField(values[i]);
  }
 
  std::string n = "Sum Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Sum,worksets,fm)
{
  this->utils.setFieldData(sum,fm);
  for (std::size_t i=0; i < values.size(); ++i)
    this->utils.setFieldData(values[i],fm);

  cell_data_size = sum.size() / sum.fieldTag().dataLayout().dimension(0);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Sum,workset)
{ 
  std::size_t length = workset.num_cells * cell_data_size;

  for (std::size_t i = 0; i < length; ++i) {
    sum[i] = 0.0;
    for (std::size_t j = 0; j < values.size(); ++j)
      sum[i] += (values[j])[i];
  }

}

//**********************************************************************

}

#endif
