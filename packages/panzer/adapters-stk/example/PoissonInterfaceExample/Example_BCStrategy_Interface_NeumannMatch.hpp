// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EXAMPLE_BCSTRATEGY_INTERFACE_NEUMANNMATCH_HPP
#define PANZER_EXAMPLE_BCSTRATEGY_INTERFACE_NEUMANNMATCH_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Interface_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"

namespace Example {

template <typename EvalT>
class BCStrategy_Interface_NeumannMatch : public panzer::BCStrategy_Interface_DefaultImpl<EvalT> {
public:    
    
  BCStrategy_Interface_NeumannMatch(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);
    
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

};

}

#include "Example_BCStrategy_Interface_NeumannMatch_impl.hpp"

#endif
