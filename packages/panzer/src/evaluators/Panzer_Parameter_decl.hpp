#ifndef PANZER_EVALUATOR_PARAMETER_DECL_HPP
#define PANZER_EVALUATOR_PARAMETER_DECL_HPP

#include "Panzer_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_RCP.hpp"

#include "Panzer_ParameterLibrary.hpp"

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
	      const double initial_value,
	      panzer::ParamLib& param_lib);

    #ifdef HAVE_STOKHOS
    /** Setup a stochastic parameter.
      *
      * \param[in] name Name of parameter and evaluated field
      * \param[in] data_layout Data layout for evaluated field, sized (Cell,Point)
      * \param[in] sg_initial_value Initial value for stochastic parameters
      * \param[in] expansion Expansion to use when constructing the stochastic scalar
      * \param[in] param_lib Parameter library to register the scalar parameter with
      */
    Parameter(const std::string name,
	      const Teuchos::RCP<PHX::DataLayout>& data_layout,
	      const std::vector<double> & sg_initial_value,
              const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
	      panzer::ParamLib& param_lib);
    #endif
    
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

#endif
