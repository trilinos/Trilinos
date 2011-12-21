#ifndef PANZER_EVALUATOR_PARAMETER_HPP
#define PANZER_EVALUATOR_PARAMETER_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {
    
  template <typename EvalT> class ScalarParameterEntry;

//! Constant parameter from sacado parameter library
  template<typename EvalT, typename Traits>
  class Parameter : 
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits> {
  
  public:
    
    Parameter(const std::string name,
	      const Teuchos::RCP<PHX::DataLayout>& data_layout,
	      const double initial_value);
    
    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& vm);
    
    void evaluateFields(typename Traits::EvalData ud);
    
  private:
    
    typedef typename EvalT::ScalarT ScalarT;
    
    PHX::MDField<ScalarT, Cell, Point> target_field;
    
    std::size_t cell_data_size;
    
    ScalarT initial_value;
    
    Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > param;
    
  };
  
}

#include "Panzer_ParameterT.hpp"

#endif
