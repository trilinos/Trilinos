#ifndef __PANZER_STK_GatherFields_HPP__
#define __PANZER_STK_GatherFields_HPP__

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"

namespace panzer_stk {

template<typename EvalT, typename Traits> class GatherFields;

// **************************************************************
// **************************************************************
// * Specializations
// **************************************************************
// **************************************************************


// **************************************************************
// Residual 
// **************************************************************
template<typename Traits>
class GatherFields<panzer::Traits::Residual,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Residual, Traits> {
   
  
public:
  
  GatherFields(const Teuchos::RCP<const STK_Interface> & mesh,const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

private:
  typedef typename panzer::Traits::Residual::ScalarT ScalarT;
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;

  std::vector< PHX::MDField<ScalarT,panzer::Cell,panzer::NODE> > gatherFields_;
  std::vector<VariableField*> stkFields_;
 
  Teuchos::RCP<const STK_Interface> mesh_;

  GatherFields();
};

// **************************************************************
// Jacobian
// **************************************************************
template<typename Traits>
class GatherFields<panzer::Traits::Jacobian,Traits>
  : public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<panzer::Traits::Jacobian, Traits> {
  
public:
  GatherFields(const Teuchos::RCP<const STK_Interface> & mesh,const Teuchos::ParameterList& p);
  
  void postRegistrationSetup(typename Traits::SetupData d,
			     PHX::FieldManager<Traits>& vm);
  
  void evaluateFields(typename Traits::EvalData d);

private:

  typedef typename panzer::Traits::Jacobian::ScalarT ScalarT;

  std::vector<std::string> fieldIds_; // field names needing reading

  std::vector< PHX::MDField<ScalarT,panzer::Cell,panzer::NODE> > gatherFields_;

  GatherFields();
};

}

// **************************************************************

#include "Panzer_STK_GatherFieldsT.hpp"

// **************************************************************
#endif
