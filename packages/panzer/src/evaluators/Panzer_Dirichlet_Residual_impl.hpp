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

#ifndef PANZER_DIRICHLET_RESIDUAL_IMPL_HPP
#define PANZER_DIRICHLET_RESIDUAL_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(DirichletResidual,p)
{
  std::string residual_name = p.get<std::string>("Residual Name");
  std::string dof_name = p.get<std::string>("DOF Name"); 
  std::string value_name = p.get<std::string>("Value Name");

  Teuchos::RCP<PHX::DataLayout> data_layout = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout"); 

  residual = PHX::MDField<ScalarT>(residual_name, data_layout);
  dof = PHX::MDField<ScalarT>(dof_name, data_layout);
  value = PHX::MDField<ScalarT>(value_name, data_layout);
  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
  this->addDependentField(value);
 
  std::string n = "Dirichlet Residual Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DirichletResidual,worksets,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(dof,fm);
  this->utils.setFieldData(value,fm);

  cell_data_size = 
    residual.size() / residual.fieldTag().dataLayout().dimension(0);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DirichletResidual,workset)
{ 
  std::size_t length = workset.num_cells * cell_data_size;

  for (std::size_t i = 0; i < length; ++i)
    residual[i] = dof[i] - value[i];
}

//**********************************************************************

}

#endif
