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

#ifndef POINT_EVALUATOR
#define POINT_EVALUATOR

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

template <typename ScalarT>
class PointEvaluation {
public:
   virtual void evaluateContainer(const Intrepid::FieldContainer<double> & points,
                                  PHX::MDField<ScalarT> & field) const = 0;
};

PHX_EVALUATOR_CLASS(RandomFieldEvaluator)

  PHX::MDField<ScalarT> field;

public:
  RandomFieldEvaluator(const std::string & name,
                       const Teuchos::RCP<PHX::DataLayout> & dl)
     : field(name,dl) { this->addEvaluatedField(field); }
PHX_EVALUATOR_CLASS_END

//**********************************************************************
PHX_EVALUATOR_CTOR(RandomFieldEvaluator,p)
{
   TEUCHOS_ASSERT(false); // don't do this
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(RandomFieldEvaluator,sd,fm)
{
  this->utils.setFieldData(field,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(RandomFieldEvaluator,workset)
{
   for(int i=0;i<field.size();i++)
      field[i] = double(std::rand())/double(RAND_MAX);
}

//**********************************************************************

#endif
