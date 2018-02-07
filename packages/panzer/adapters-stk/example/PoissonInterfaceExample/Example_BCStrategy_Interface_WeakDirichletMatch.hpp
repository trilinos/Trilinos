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

#ifndef PANZER_EXAMPLE_BCSTRATEGY_INTERFACE_WEAKDIRICHLETMATCH_HPP
#define PANZER_EXAMPLE_BCSTRATEGY_INTERFACE_WEAKDIRICHLETMATCH_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Interface_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"

namespace Example {

template <typename EvalT>
class BCStrategy_Interface_WeakDirichletMatch : public panzer::BCStrategy_Interface_DefaultImpl<EvalT> {
public:
  BCStrategy_Interface_WeakDirichletMatch(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

  void setup(const panzer::PhysicsBlock& side_pb,
             const Teuchos::ParameterList& user_data);

  void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                  const panzer::PhysicsBlock& pb,
                                  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                  const Teuchos::ParameterList& models,
                                  const Teuchos::ParameterList& user_data) const;

  virtual void buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                              const panzer::PhysicsBlock& side_pb,
                                                              const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                              const Teuchos::ParameterList& user_data) const;

  virtual void postRegistrationSetup(typename panzer::Traits::SetupData d,
                                     PHX::FieldManager<panzer::Traits>& vm);

  virtual void evaluateFields(typename panzer::Traits::EvalData d);

private:
  std::vector<std::string> paramName;
  double value;
  double temp;
  std::string other_dof_name;

  static void setSumValues(Teuchos::ParameterList& p,
                           const std::string value_name1, const double scalar1,
                           const std::string value_name2, const double scalar2);
};

}

#include "Example_BCStrategy_Interface_WeakDirichletMatch_impl.hpp"

#endif
