// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHX_EVALUATOR_MACROS_H
#define PHX_EVALUATOR_MACROS_H

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_ParameterList.hpp"

// **********************************************************************
//! Macro definition of an evaluator class
#define PHX_EVALUATOR_CLASS(NAME)					\
  									\
  template<typename EvalT, typename Traits>				\
  class NAME : public PHX::EvaluatorWithBaseImpl<Traits>,		\
	       public PHX::EvaluatorDerived<EvalT, Traits>  {		\
    									\
  public:								\
    									\
    NAME(const Teuchos::ParameterList& p);				\
    									\
    void postRegistrationSetup(PHX::FieldManager<Traits>& fm);		\
    									\
    void evaluateFields(typename Traits::EvalData d);			\
    									\
  private:								\
    									\
    typedef typename EvalT::ScalarT ScalarT;

// **********************************************************************
//! Macro definition of an evaluator class with pre/post evaluate methods
#define PHX_EVALUATOR_CLASS_PP(NAME)					\
  									\
  template<typename EvalT, typename Traits>				\
  class NAME : public PHX::EvaluatorWithBaseImpl<Traits>,		\
	       public PHX::EvaluatorDerived<EvalT, Traits>  {		\
    									\
  public:								\
    									\
    NAME(const Teuchos::ParameterList& p);				\
    									\
    void postRegistrationSetup(PHX::FieldManager<Traits>& fm);		\
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
#define PHX_EVALUATOR_CLASS_END };


// **********************************************************************
//! Macro definition of an evaluator constructor
#define PHX_EVALUATOR_CTOR(NAME,PLIST)					\
  template<typename EvalT, typename Traits>				\
  NAME <EvalT, Traits>::NAME(const Teuchos::ParameterList& PLIST)


// **********************************************************************
//! Macro definition of an evaluator constructor for a namespaced object
#define PHX_EVALUATOR_CTOR_NAMESPACE(NAMESPACE,NAME,PLIST)		\
  template<typename EvalT, typename Traits>				\
  NAMESPACE::NAME<EvalT, Traits>::NAME(const Teuchos::ParameterList& PLIST)


// **********************************************************************
//! Macro definition for the evaluator postRegistrationSetup method
#define PHX_POST_REGISTRATION_SETUP(NAME,FIELD_MANAGER)			\
  template<typename EvalT, typename Traits>				\
  void NAME<EvalT, Traits>::						\
  postRegistrationSetup(PHX::FieldManager<Traits>& FIELD_MANAGER)


// **********************************************************************
//! Macro definition for the evaluator evaluateFields method
#define PHX_EVALUATE_FIELDS(NAME,EVAL_DATA)                             \
  template<typename EvalT, typename Traits>				\
  void NAME<EvalT, Traits>::						\
  evaluateFields(typename Traits::EvalData EVAL_DATA)


// **********************************************************************
//! Macro definition for the evaluator evaluateFields method
#define PHX_PRE_EVALUATE_FIELDS(NAME,PRE_EVAL_DATA)			\
  template<typename EvalT, typename Traits>				\
  void NAME<EvalT, Traits>::						\
  preEvaluate(typename Traits::PreEvalData PRE_EVAL_DATA)


// **********************************************************************
//! Macro definition for the evaluator evaluateFields method
#define PHX_POST_EVALUATE_FIELDS(NAME,POST_EVAL_DATA)			\
  template<typename EvalT, typename Traits>				\
  void NAME<EvalT, Traits>::						\
  postEvaluate(typename Traits::PostEvalData POSTEVAL_DATA)

// **********************************************************************

#endif 
