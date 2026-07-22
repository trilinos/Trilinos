// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATORS_REGISTRAR_HPP
#define PANZER_EVALUATORS_REGISTRAR_HPP

#include "Phalanx_FieldManager.hpp"

#include "Panzer_Evaluator_WithBaseImpl.hpp"

namespace panzer {

/** Classes that call PHX::FieldManager::registerEvaluator on
  * panzer::EvaluatorWithBaseImpl objects inherit from this class to wrap the
  * registerEvaluator call. This class injects the WorksetDetails index into the
  * evaluator.
 */
class EvaluatorsRegistrar {
public:
  //! Set the WorksetDetails index in all evaluators registered through
  //! EquationSetBase::registerEvaluator. The details index can be set multiple
  //! times. The current value applies at registration. Return the previous
  //! value.
  int setDetailsIndex(const int details_index) {
    int old_di = details_index_;
    details_index_ = details_index;
    return old_di;
  }
  //! Get the WorksetDetails index.
  int getDetailsIndex() const { return details_index_; }

protected:
  //! Default ctor initializes WorksetDetails index to 0.
  EvaluatorsRegistrar() : details_index_(0) {}
  virtual ~EvaluatorsRegistrar() {}

  //! Register the evaluator and initialize it with any information in
  //! panzer::EvaluatorWithBaseImpl.
  template <typename EvalT>
  void registerEvaluator(PHX::FieldManager<panzer::Traits>& fm,
                         const Teuchos::RCP< PHX::Evaluator<panzer::Traits> >& op) const;

private:
  int details_index_;
};

template<typename EvalT>
void EvaluatorsRegistrar::
registerEvaluator(PHX::FieldManager<panzer::Traits>& fm,
                  const Teuchos::RCP< PHX::Evaluator<panzer::Traits> >& op) const
{
  Teuchos::RCP< panzer::EvaluatorWithBaseImpl<panzer::Traits> >
    pop = Teuchos::rcp_dynamic_cast< panzer::EvaluatorWithBaseImpl<panzer::Traits> >(op);
  // Temporarily allow casting failure so that Charon continues to work.
#if 0
  TEUCHOS_TEST_FOR_EXCEPTION(pop.is_null(), std::runtime_error,
                             op->getName() + " does not inherit from panzer::EvaluatorWithBaseImpl.");
  pop->setDetailsIndex(details_index_);
#else
  if (Teuchos::nonnull(pop))
    pop->setDetailsIndex(details_index_);
#endif
  fm.template registerEvaluator<EvalT>(op);
}

}

#endif
