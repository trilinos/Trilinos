#ifndef PANZER_EVALUATOR_GATHER_BASIS_COORDINATES_DECL_HPP
#define PANZER_EVALUATOR_GATHER_BASIS_COORDINATES_DECL_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CloneableEvaluator.hpp"

#include "Panzer_PureBasis.hpp"

namespace panzer {

/** \brief Gathers coordinates for the basis function from the
    workset and stores them in the field manager.
*/
template<typename EvalT, typename Traits>
class GatherBasisCoordinates
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits> {
public:
  
  GatherBasisCoordinates(const panzer::PureBasis & basis);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

  static std::string fieldName(const std::string & basisName);

private:
  typedef typename EvalT::ScalarT ScalarT;

  std::string basisName_;
  std::vector<std::string>::size_type basisIndex_; 
  PHX::MDField<ScalarT,Cell,BASIS,Dim> basisCoordinates_;

  GatherBasisCoordinates();
};

}

// **************************************************************
#endif
