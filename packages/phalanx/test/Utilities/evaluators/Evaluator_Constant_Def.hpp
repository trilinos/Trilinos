// @HEADER
// @HEADER

#include "Teuchos_ParameterList.hpp"

//**********************************************************************
template<typename EvalT, typename Traits>
Constant<EvalT, Traits>::Constant(Teuchos::ParameterList& p) :
  value( p.get<double>("Value") ),
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);
  
  std::string n = "Constant Provider: " + constant.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
Constant<EvalT, Traits>::~Constant()
{ }

//**********************************************************************
template<typename EvalT, typename Traits>
void Constant<EvalT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,vm);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Constant<EvalT, Traits>::evaluateFields(typename Traits::EvalData d)
{ }

//**********************************************************************
