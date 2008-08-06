
#ifndef PHX_EVALUATOR_WITHBASEIMPL_H
#define PHX_EVALUATOR_WITHBASEIMPL_H

#include <vector>

#include "Phalanx_Evaluator.hpp"
#include "Phalanx_Field.hpp"

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

    virtual void addEvaluatedField(const PHX::FieldTag& v);

    template<typename DataT>
    void addEvaluatedField(const PHX::Field<DataT>& h);

    virtual void addDependentField(const PHX::FieldTag& v);

    template<typename DataT>
    void addDependentField(const PHX::Field<DataT>& h);

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
