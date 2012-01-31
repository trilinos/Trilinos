#ifndef USER_APP_CONSTANT_MODEL_HPP
#define USER_APP_CONSTANT_MODEL_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_Field.hpp"

namespace user_app {
    
PHX_EVALUATOR_CLASS(ConstantModel)
  
  ScalarT value;
  
  PHX::Field<ScalarT> constant;
  
PHX_EVALUATOR_CLASS_END


#ifdef HAVE_STOKHOS

template<typename Traits>                        
class ConstantModel<typename Traits::SGResidual,Traits> : public PHX::EvaluatorWithBaseImpl<Traits>,         
             public PHX::EvaluatorDerived<typename Traits::SGResidual, Traits>  {    
public:                                                      
  ConstantModel(const Teuchos::ParameterList& p);                    
                                                           
  void postRegistrationSetup(typename Traits::SetupData d,
                             PHX::FieldManager<Traits>& fm);      
                                                                   
  void evaluateFields(typename Traits::EvalData d) {}
                                                                     
private:                                                              
                                                                      
  typedef typename Traits::SGResidual::ScalarT ScalarT;

  ScalarT value;
  
  PHX::Field<ScalarT> constant;
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion;
};

template<typename Traits>                        
class ConstantModel<typename Traits::SGJacobian,Traits> : public PHX::EvaluatorWithBaseImpl<Traits>,         
             public PHX::EvaluatorDerived<typename Traits::SGJacobian, Traits>  {    
public:                                                      
  ConstantModel(const Teuchos::ParameterList& p);                    
                                                           
  void postRegistrationSetup(typename Traits::SetupData d,
                             PHX::FieldManager<Traits>& fm);      
                                                                   
  void evaluateFields(typename Traits::EvalData d) {}
                                                                     
private:                                                              
                                                                      
  typedef typename Traits::SGJacobian::ScalarT ScalarT;

  ScalarT value;
  
  PHX::Field<ScalarT> constant;
  Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion;
};

#endif

}

#include "user_app_ConstantModel_impl.hpp"

#endif
