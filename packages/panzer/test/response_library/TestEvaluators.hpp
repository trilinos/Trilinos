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

#ifndef UNIT_VALUE_EVALUATOR
#define UNIT_VALUE_EVALUATOR

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace panzer {

PHX_EVALUATOR_CLASS(TestEvaluator)

  PHX::MDField<ScalarT,Cell> dogValues;
  PHX::MDField<ScalarT,Cell> hrsValues;

PHX_EVALUATOR_CLASS_END

//**********************************************************************
PHX_EVALUATOR_CTOR(TestEvaluator,p)
{
  // Read from parameters
  int worksetSize = p.get<int>("Workset Size");

  Teuchos::RCP<PHX::DataLayout> dl = Teuchos::rcp(new PHX::MDALayout<Cell>(worksetSize));
  // grab information from quadrature rule
  dogValues = PHX::MDField<ScalarT,Cell>("Dog", dl);
  hrsValues = PHX::MDField<ScalarT,Cell>("Horse", dl);

  this->addEvaluatedField(dogValues);
  this->addEvaluatedField(hrsValues);
  
  std::string n = "TestEvaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(TestEvaluator,sd,fm)
{
  this->utils.setFieldData(dogValues,fm);
  this->utils.setFieldData(hrsValues,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(TestEvaluator,workset)
{ 
   double extra = 0.0;
   if(workset.block_id=="block_1")
      extra = 44.3;

   for(int i=0;i<dogValues.dimension(0);i++) {
      dogValues(i) = double(i) + 1.0 + extra;
      hrsValues(i) = -double(i) - 5.5 + extra;
   }
}

//**********************************************************************

}

#endif
