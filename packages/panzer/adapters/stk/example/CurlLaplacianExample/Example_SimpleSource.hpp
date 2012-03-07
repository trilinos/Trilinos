#ifndef EXAMPLE_SIMPLE_SOURCE_DECL_HPP
#define EXAMPLE_SIMPLE_SOURCE_DECL_HPP

#include "Panzer_config.hpp"

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

#include <string>

namespace Example {
    
  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

/** A source for the curl Laplacian that results in the solution
  */
template<typename EvalT, typename Traits>
class SimpleSource : public PHX::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>  {

public:
    SimpleSource(const std::string & name,
                       const panzer::IntegrationRule & ir);
                                                                        
    void postRegistrationSetup(typename Traits::SetupData d,           
                               PHX::FieldManager<Traits>& fm);        
                                                                     
    void evaluateFields(typename Traits::EvalData d);               


private:
  typedef typename EvalT::ScalarT ScalarT;

  // Simulation source
  PHX::MDField<ScalarT,Cell,Point,Dim> source;
  int ir_degree, ir_index;
};

}

#include "Example_SimpleSource_impl.hpp"

#endif
