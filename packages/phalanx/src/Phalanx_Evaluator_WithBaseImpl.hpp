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

    virtual void postRegistrationSetup(PHX::FieldManager<Traits>& vm) = 0;

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
