#ifndef PANZER_CONVECTION_T_HPP
#define PANZER_CONVECTION_T_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace user_app {

//**********************************************************************
PHX_EVALUATOR_CTOR(Convection,p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;

  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");

  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  
  conv = MDField<ScalarT,Cell,Point>(p.get<string>("Operator Name"),scalar);

  a = MDField<ScalarT,Cell,Point,Dim>(p.get<string>("A Name"),vector);

  grad_x = 
    MDField<ScalarT,Cell,Point,Dim>(p.get<string>("Gradient Name"),vector);

  multiplier = p.get<double>("Multiplier");

  this->addEvaluatedField(conv);

  this->addDependentField(a);
  this->addDependentField(grad_x);
  
  std::string n = "Convection: " + conv.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Convection,worksets,fm)
{
  this->utils.setFieldData(conv,fm);
  this->utils.setFieldData(a,fm);
  this->utils.setFieldData(grad_x,fm);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Convection,workset)
{ 
  typedef typename PHX::MDField<ScalarT,Cell,Point>::size_type size_type;
  
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {    
    for (size_type point = 0; point < conv.dimension(1); ++point) {
      
      conv(cell,point) = 0.0;
	
      for (size_type dim = 0; dim < a.dimension(2); ++dim)
	conv(cell,point) += a(cell,point,dim) * grad_x(cell,point,dim);
      
      conv(cell,point) *= multiplier;
      
    }
  }

}

//**********************************************************************

}

#endif
