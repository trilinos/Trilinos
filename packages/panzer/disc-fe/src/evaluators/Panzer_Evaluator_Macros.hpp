// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_MACROS_HPP
#define PANZER_EVALUATOR_MACROS_HPP

#include "Panzer_Evaluator_WithBaseImpl.hpp"

// **********************************************************************
//! Macro definition of an evaluator class
#define PANZER_EVALUATOR_CLASS(NAME)					\
                                                                        \
  template<typename EvalT, typename Traits>				\
  class NAME : public panzer::EvaluatorWithBaseImpl<Traits>,		\
	       public PHX::EvaluatorDerived<EvalT, Traits>  {		\
    									\
  public:								\
    									\
    NAME(const Teuchos::ParameterList& p);				\
    									\
    void postRegistrationSetup(typename Traits::SetupData d,            \
                               PHX::FieldManager<Traits>& fm);   	\
    									\
    void evaluateFields(typename Traits::EvalData d);			\
    									\
  private:								\
    									\
    typedef typename EvalT::ScalarT ScalarT;

// **********************************************************************
//! Macro definition of an evaluator class with pre/post evaluate methods
#define PANZER_EVALUATOR_CLASS_PP(NAME)					\
  									\
  template<typename EvalT, typename Traits>				\
  class NAME : public panzer::EvaluatorWithBaseImpl<Traits>,		\
	       public PHX::EvaluatorDerived<EvalT, Traits>  {		\
    									\
  public:								\
    									\
    NAME(const Teuchos::ParameterList& p);				\
    									\
    void postRegistrationSetup(typename Traits::SetupData d,            \
                               PHX::FieldManager<Traits>& fm);		\
    									\
    void evaluateFields(typename Traits::EvalData d);			\
									\
    void preEvaluate(typename Traits::PreEvalData d);			\
									\
    void postEvaluate(typename Traits::PostEvalData d);			\
									\
  private:								\
									\
    typedef typename EvalT::ScalarT ScalarT;

// **********************************************************************
//! Macro definition to end an evaluator class declaration
#define PANZER_EVALUATOR_CLASS_END };

#endif
