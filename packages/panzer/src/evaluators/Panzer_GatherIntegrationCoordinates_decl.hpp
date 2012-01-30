#ifndef PANZER_EVALUATOR_GATHER_INTEGRATION_COORDINATES_DECL_HPP
#define PANZER_EVALUATOR_GATHER_INTEGRATION_COORDINATES_DECL_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_PureBasis.hpp"

namespace panzer {

/** \brief Gathers coordinates for the quadrature from the
    workset and stores them in the field manager.
*/
template<typename EvalT, typename Traits>
class GatherIntegrationCoordinates
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  
  GatherIntegrationCoordinates(const panzer::IntegrationRule & quad);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  static std::string fieldName(int degree);

private:
  typedef typename EvalT::ScalarT ScalarT;

  int quadDegree_;
  std::vector<int>::size_type quadIndex_; 
  PHX::MDField<ScalarT,Cell,Point,Dim> quadCoordinates_;

  GatherIntegrationCoordinates();
};

}

// **************************************************************
#endif
