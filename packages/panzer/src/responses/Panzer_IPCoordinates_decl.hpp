#ifndef PANZER_EVALUATOR_IP_COORDINATES_DECL_HPP
#define PANZER_EVALUATOR_IP_COORDINATES_DECL_HPP

#include "Panzer_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {
  
  /** Evaluates the coordinates on fist evaluation
   *
   *  NOTE: This assumes a static mesh and does not recompute the
   *  coordinates on each response evaluation.  It is simple to change
   *  in the future if needed.
   */
  template<typename EvalT, typename Traits>
  class IPCoordinates : 
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits> {
  
  public:
    
    struct Coordinate {
      int lid;
      std::vector<double> val;
    };

    IPCoordinates(int ir_order,
		  const Teuchos::RCP<std::vector<Coordinate> >& coords);
    
    void postRegistrationSetup(typename Traits::SetupData d,
			       PHX::FieldManager<Traits>& vm);
    
    void preEvaluate(typename Traits::PreEvalData data);
    void evaluateFields(typename Traits::EvalData ud);
    void postEvaluate(typename Traits::PostEvalData data);

    typedef typename EvalT::ScalarT ScalarT;
    
    const PHX::MDField<ScalarT,Cell,IP,Dim> getEvaluatedField() const;

  private:
    
    PHX::MDField<ScalarT,Cell,IP,Dim> dummy_field;
 
    bool first_evaluation;

    //! integration rule order
    int ir_order;

    //! integration rule index
    int ir_index;

    //! Coordinate vector from ResponseData object
    Teuchos::RCP<std::vector<Coordinate> > coords;
  };
  
}

#endif
