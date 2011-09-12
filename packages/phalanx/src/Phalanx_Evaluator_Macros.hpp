// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
#define PHX_POST_REGISTRATION_SETUP(NAME,SETUP_DATA,FIELD_MANAGER)	\
  template<typename EvalT, typename Traits>				\
  void NAME<EvalT, Traits>::						\
  postRegistrationSetup(typename Traits::SetupData SETUP_DATA,          \
			PHX::FieldManager<Traits>& FIELD_MANAGER)


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
