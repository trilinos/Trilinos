#ifndef PANZER_PointValues_Evaluator_IMPL_HPP
#define PANZER_PointValues_Evaluator_IMPL_HPP

#include <algorithm>
#include "Panzer_PointRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(PointValues_Evaluator,p)
{

  Teuchos::RCP<const panzer::PointRule> pointRule 
     = p.get< Teuchos::RCP<const panzer::PointRule> >("Point Rule");
  Teuchos::RCP<const Intrepid::FieldContainer<double> > userArray
     = p.get<Teuchos::RCP<const Intrepid::FieldContainer<double> > >("Point Array");

  initialize(pointRule,userArray);
}

//**********************************************************************
template <typename EvalT, typename TraitsT>
PointValues_Evaluator<EvalT,TraitsT>::PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                            const Teuchos::RCP<const Intrepid::FieldContainer<double> > & userArray)
{
  initialize(pointRule,userArray);
}

//**********************************************************************
template <typename EvalT, typename TraitsT>
void PointValues_Evaluator<EvalT,TraitsT>::initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                      const Teuchos::RCP<const Intrepid::FieldContainer<double> > & userArray)
{
  TEUCHOS_ASSERT(userArray->rank()==2);

  panzer::MDFieldArrayFactory af(pointRule->getName()+"_");
       
  // copy user array data
  refPointArray = Intrepid::FieldContainer<double>(userArray->dimension(0),userArray->dimension(1));
  TEUCHOS_ASSERT(refPointArray.size()==userArray->size());
  for(int i=0;i<userArray->size();i++)
     refPointArray[i] = (*userArray)[i]; 

  // setup all fields to be evaluated and constructed
  pointValues.setupArrays(pointRule,af);

  // the field manager will allocate all of these field
  this->addEvaluatedField(pointValues.coords_ref);
  this->addEvaluatedField(pointValues.node_coordinates);
  this->addEvaluatedField(pointValues.jac);
  this->addEvaluatedField(pointValues.jac_inv);
  this->addEvaluatedField(pointValues.jac_det);
  this->addEvaluatedField(pointValues.point_coords);

  std::string n = "PointValues_Evaluator: " + pointRule->getName();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(PointValues_Evaluator,sd,fm)
{
  // setup the pointers for evaluation
  this->utils.setFieldData(pointValues.coords_ref,fm);
  this->utils.setFieldData(pointValues.node_coordinates,fm);
  this->utils.setFieldData(pointValues.jac,fm);
  this->utils.setFieldData(pointValues.jac_inv,fm);
  this->utils.setFieldData(pointValues.jac_det,fm);
  this->utils.setFieldData(pointValues.point_coords,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(PointValues_Evaluator,workset)
{ 
  // evaluate the point values (construct jacobians etc...)
  pointValues.evaluateValues(workset.cell_vertex_coordinates,refPointArray);
}

//**********************************************************************

}

#endif
