#ifndef PHX_SCALAR_CONTAINER_BASE_HPP
#define PHX_SCALAR_CONTAINER_BASE_HPP

#include "Phalanx_Evaluator_Manager.hpp"

namespace PHX {

  template<typename Traits> class FieldManager;

  template<typename Traits>
  class ScalarContainerBase {

  public:

    ScalarContainerBase();

    virtual ~ScalarContainerBase();

    virtual void requireField(const PHX::FieldTag& v);

    virtual void 
    registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& p);

    virtual void postRegistrationSetup(std::size_t max_num_cells,
				       PHX::FieldManager<Traits>& vm) = 0;

    virtual void evaluateFields(typename Traits::EvalData d) = 0;

    virtual void preEvaluate(typename Traits::PreEvalData d) = 0;

    virtual void postEvaluate(typename Traits::PostEvalData d) = 0;

    virtual void print(std::ostream& os) const = 0;
    
  protected:
    
    PHX::EvaluatorManager<Traits> vp_manager_;

  };

  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
		  const PHX::ScalarContainerBase<Traits>& sc);
  
}

#include "Phalanx_EvaluationContainer_Base_Def.hpp"

#endif 
