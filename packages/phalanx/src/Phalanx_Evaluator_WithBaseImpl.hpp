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


#ifndef PHX_EVALUATOR_WITHBASEIMPL_H
#define PHX_EVALUATOR_WITHBASEIMPL_H

#include <vector>

#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_MDField.hpp"

namespace PHX {

  /*! @brief Class that implements helper functions for the pure virtual PHX::Evaluator class.
   
      This class implements code that would essentially be repeated in
      each Evaluator class, making it quicker for developers to add
      new evaluators.  All field evaluators should inherit from this
      class if possible instead of the base class so they don't have
      to code the same boilerplate in all evaluators, but this is not
      mandatory.
  */
  template <typename Traits>
  class EvaluatorWithBaseImpl : public PHX::Evaluator<Traits> {

  public:

    EvaluatorWithBaseImpl(const std::string& evaluator_name);

    EvaluatorWithBaseImpl();

    virtual ~EvaluatorWithBaseImpl();

    virtual void addEvaluatedField(const PHX::FieldTag& ft);

    template<typename DataT>
    void addEvaluatedField(const PHX::Field<DataT>& f);

    template<typename DataT,
	     typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	     typename Tag4, typename Tag5, typename Tag6, typename Tag7>
    void addEvaluatedField(const PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,
			   Tag4,Tag5,Tag6,Tag7>& f);

    virtual void addDependentField(const PHX::FieldTag& ft);

    template<typename DataT>
    void addDependentField(const PHX::Field<DataT>& f);

    template<typename DataT,
	     typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	     typename Tag4, typename Tag5, typename Tag6, typename Tag7>
    void addDependentField(const PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,
			   Tag4,Tag5,Tag6,Tag7>& f);

    virtual void setName(const std::string& name);

    virtual void postRegistrationSetup(typename Traits::SetupData d,
				       PHX::FieldManager<Traits>& vm) = 0;

    virtual const std::vector< Teuchos::RCP<FieldTag> >& 
    evaluatedFields() const;

    virtual const std::vector< Teuchos::RCP<FieldTag> >& 
    dependentFields() const;

    virtual void evaluateFields(typename Traits::EvalData d) = 0;

    virtual void preEvaluate(typename Traits::PreEvalData d);

    virtual void postEvaluate(typename Traits::PostEvalData d);

    virtual const std::string& getName() const;

  private:

    std::vector< Teuchos::RCP<FieldTag> > evaluated_;

    std::vector< Teuchos::RCP<FieldTag> > required_;

    std::string name_;
  };

}

#include "Phalanx_Evaluator_WithBaseImpl_Def.hpp"

#endif
