// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHX_EVALUATION_CONTAINER_BASE_HPP
#define PHX_EVALUATION_CONTAINER_BASE_HPP

#include <cstddef>
#include <string>
#include <map>
#include "Phalanx_DAG_Manager.hpp"

namespace PHX {

  template<typename Traits> class FieldManager;
  class MemoryManager;

  template<typename Traits>
  class EvaluationContainerBase {

  public:

    EvaluationContainerBase();

    virtual ~EvaluationContainerBase();

    virtual void requireField(const PHX::FieldTag& v);

    virtual void aliasField(const PHX::FieldTag& aliasedField,
                            const PHX::FieldTag& targetField) = 0;
    
    virtual void 
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);

    virtual void postRegistrationSetup(typename Traits::SetupData d,
				       PHX::FieldManager<Traits>& vm,
                                       const bool& buildDeviceDAG,
                                       const bool& minimizeDAGMemoryUse,
                                       const PHX::MemoryManager* const memoryManager) = 0;

    virtual void evaluateFields(typename Traits::EvalData d) = 0;

    virtual void preEvaluate(typename Traits::PreEvalData d) = 0;

    virtual void postEvaluate(typename Traits::PostEvalData d) = 0;

    virtual void writeGraphvizFile(const std::string filename,
				   bool writeEvaluatedFields,
				   bool writeDependentFields,
				   bool debugRegisteredEvaluators) const;

    virtual const std::string evaluationType() const = 0;

    virtual void print(std::ostream& os) const = 0;
    
  protected:
    
    PHX::DagManager<Traits> dag_manager_;

  };

  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::EvaluationContainerBase<Traits>& sc);
  
}

#include "Phalanx_EvaluationContainer_Base_Def.hpp"

#endif 
