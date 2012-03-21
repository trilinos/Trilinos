#ifndef PANZER_BasisValues_Evaluator_IMPL_HPP
#define PANZER_BasisValues_Evaluator_IMPL_HPP

#include <algorithm>
#include "Panzer_PointRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(BasisValues_Evaluator,p)
{
  Teuchos::RCP<const panzer::PointRule> pointRule 
     = p.get< Teuchos::RCP<const panzer::PointRule> >("Point Rule");
  Teuchos::RCP<const panzer::PureBasis> inBasis
     = p.get<Teuchos::RCP<const panzer::PureBasis> >("Basis");

  initialize(pointRule,inBasis);
}

//**********************************************************************
template <typename EvalT, typename TraitsT>
BasisValues_Evaluator<EvalT,TraitsT>::BasisValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                            const Teuchos::RCP<const panzer::PureBasis> & inBasis)
{
  initialize(pointRule,inBasis);
}

//**********************************************************************
template <typename EvalT, typename TraitsT>
void BasisValues_Evaluator<EvalT,TraitsT>::initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                      const Teuchos::RCP<const panzer::PureBasis> & inBasis)
{
  basis = inBasis;

  panzer::MDFieldArrayFactory af_pv(pointRule->getName()+"_");
  panzer::MDFieldArrayFactory af_bv(basis->name()+"_"+pointRule->getName()+"_");
       
  // setup all fields to be evaluated and constructed
  pointValues.setupArrays(pointRule,af_pv);

  // the field manager will allocate all of these field
  this->addDependentField(pointValues.coords_ref);
  this->addDependentField(pointValues.node_coordinates);
  this->addDependentField(pointValues.jac);
  this->addDependentField(pointValues.jac_inv);
  this->addDependentField(pointValues.jac_det);
  this->addDependentField(pointValues.point_coords);

  // setup all fields to be evaluated and constructed
  Teuchos::RCP<panzer::BasisIRLayout> layout = Teuchos::rcp(new panzer::BasisIRLayout(basis,*pointRule));
  basisValues.setupArrays(layout,af_bv);

  // the field manager will allocate all of these field

  this->addEvaluatedField(basisValues.basis_ref);      
  this->addEvaluatedField(basisValues.basis);           

  if(basis->getElementSpace()==panzer::PureBasis::HGRAD) {
    this->addEvaluatedField(basisValues.grad_basis_ref);   
    this->addEvaluatedField(basisValues.grad_basis);        
  }

  if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
    this->addEvaluatedField(basisValues.curl_basis_ref);     
    this->addEvaluatedField(basisValues.curl_basis);          
  }

  std::string n = "BasisValues_Evaluator: " +basis->name() + "_" + pointRule->getName();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(BasisValues_Evaluator,sd,fm)
{
  // setup the pointers for the point values data structure
  this->utils.setFieldData(pointValues.coords_ref,fm);
  this->utils.setFieldData(pointValues.node_coordinates,fm);
  this->utils.setFieldData(pointValues.jac,fm);
  this->utils.setFieldData(pointValues.jac_inv,fm);
  this->utils.setFieldData(pointValues.jac_det,fm);
  this->utils.setFieldData(pointValues.point_coords,fm);

  // setup the pointers for the basis values data structure
  this->utils.setFieldData(basisValues.basis_ref,fm);      
  this->utils.setFieldData(basisValues.basis,fm);           

  if(basis->getElementSpace()==panzer::PureBasis::HGRAD) {
    this->utils.setFieldData(basisValues.grad_basis_ref,fm);   
    this->utils.setFieldData(basisValues.grad_basis,fm);        
  }

  if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
    this->utils.setFieldData(basisValues.curl_basis_ref,fm);     
    this->utils.setFieldData(basisValues.curl_basis,fm);          
  }
}

//**********************************************************************
PHX_EVALUATE_FIELDS(BasisValues_Evaluator,workset)
{ 
  // evaluate the point values (construct jacobians etc...)
  basisValues.evaluateValues(pointValues.coords_ref,
                             pointValues.jac,
                             pointValues.jac_det,
                             pointValues.jac_inv);
}

//**********************************************************************

}

#endif
