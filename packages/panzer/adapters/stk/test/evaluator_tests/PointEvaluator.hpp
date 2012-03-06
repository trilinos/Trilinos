#ifndef POINT_EVALUATOR
#define POINT_EVALUATOR

#include <string>
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_MDField.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

template <typename ScalarT>
class PointEvaluation {
public:
   virtual void evaluateContainer(const Intrepid::FieldContainer<double> & points,
                                  PHX::MDField<ScalarT> & field) const = 0;
};

PHX_EVALUATOR_CLASS(PointEvaluator)

  PHX::MDField<ScalarT> scalar;
  PHX::MDField<ScalarT> vectorField;

  bool isVector;
  int quad_order;
  int quad_index;
  Teuchos::RCP<const PointEvaluation<ScalarT> > function;
PHX_EVALUATOR_CLASS_END

//**********************************************************************
PHX_EVALUATOR_CTOR(PointEvaluator,p)
{
  // Read from parameters
  const std::string name = p.get<std::string>("Name");
  Teuchos::RCP<panzer::IntegrationRule> ir
     = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  isVector = p.get<bool>("Is Vector");
  function = p.get<Teuchos::RCP<const PointEvaluation<ScalarT> > >("Point Evaluator");

  quad_order = ir->cubature_degree;

  // grab information from quadrature rule
  if(isVector) {
     vectorField = PHX::MDField<ScalarT>(name, ir->dl_vector);
     this->addEvaluatedField(vectorField);
  }
  else {
     scalar = PHX::MDField<ScalarT>(name, ir->dl_scalar);
     this->addEvaluatedField(scalar);
  }


  std::string n = "PointEvaluator: " + name;
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(PointEvaluator,sd,fm)
{
  if(isVector)
     this->utils.setFieldData(vectorField,fm);
  else
     this->utils.setFieldData(scalar,fm);

  quad_index =  panzer::getIntegrationRuleIndex(quad_order,(*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(PointEvaluator,workset)
{
   if(isVector)
      function->evaluateContainer(workset.int_rules[quad_index]->ip_coordinates,vectorField);
   else 
      function->evaluateContainer(workset.int_rules[quad_index]->ip_coordinates,vectorField);
}

//**********************************************************************

#endif
