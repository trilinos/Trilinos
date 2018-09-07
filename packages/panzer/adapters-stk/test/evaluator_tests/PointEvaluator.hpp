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
#include "Phalanx_MDField.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_Evaluator_Macros.hpp"

template <typename ScalarT>
class PointEvaluation {
public:
   // virtual void evaluateContainer(const Kokkos::DynRankView<double,PHX::Device> & points,
   virtual void evaluateContainer(const PHX::MDField<ScalarT,panzer::Cell,panzer::IP,panzer::Dim> & points,
                                  PHX::MDField<ScalarT> & field) const = 0;
};

template<typename EvalT, typename Traits>
class PointEvaluator
  :
  public panzer::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    PointEvaluator(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<ScalarT> scalar;
  PHX::MDField<ScalarT> vectorField;

  bool isVector;
  int quad_order;
  int quad_index;
  Teuchos::RCP<const PointEvaluation<ScalarT> > function;
}; // end of class PointEvaluator


//**********************************************************************
template<typename EvalT, typename Traits>
PointEvaluator<EvalT, Traits>::
PointEvaluator(
  const Teuchos::ParameterList& p)
{
  // Read from parameters
  const std::string name = p.get<std::string>("Name");
  Teuchos::RCP<panzer::IntegrationRule> ir
     = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  isVector = p.get<bool>("Is Vector");
  function = p.get<Teuchos::RCP<const PointEvaluation<ScalarT> > >("Point Evaluator");

  quad_order = ir->cubature_degree;

  // grab information from quadrature rule
  if(isVector) {
     vectorField = PHX::MDField<ScalarT>(name, ir->dl_vector);
     this->addEvaluatedField(vectorField);
  }
  else {
     scalar = PHX::MDField<ScalarT>(name, ir->dl_scalar);
     this->addEvaluatedField(scalar);
  }


  std::string n = "PointEvaluator: " + name;
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
PointEvaluator<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  quad_index =  panzer::getIntegrationRuleIndex(quad_order,(*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
PointEvaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
   if(isVector)
      function->evaluateContainer(this->wda(workset).int_rules[quad_index]->ip_coordinates,vectorField);
   else 
      function->evaluateContainer(this->wda(workset).int_rules[quad_index]->ip_coordinates,vectorField);
}

//**********************************************************************

#endif
