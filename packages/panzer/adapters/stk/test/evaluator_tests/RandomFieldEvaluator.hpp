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

PHX_EVALUATOR_CLASS(RandomFieldEvaluator)

  PHX::MDField<ScalarT> field;

public:
  RandomFieldEvaluator(const std::string & name,
                       const Teuchos::RCP<PHX::DataLayout> & dl)
     : field(name,dl) { this->addEvaluatedField(field); }
PHX_EVALUATOR_CLASS_END

//**********************************************************************
PHX_EVALUATOR_CTOR(RandomFieldEvaluator,p)
{
   TEUCHOS_ASSERT(false); // don't do this
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(RandomFieldEvaluator,sd,fm)
{
  this->utils.setFieldData(field,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(RandomFieldEvaluator,workset)
{
   for(int i=0;i<field.size();i++)
      field[i] = double(std::rand())/double(RAND_MAX);
}

//**********************************************************************

#endif
