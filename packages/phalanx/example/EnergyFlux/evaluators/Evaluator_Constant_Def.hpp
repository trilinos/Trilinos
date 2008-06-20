
//**********************************************************************
template<typename ScalarT, typename Traits>
Constant<ScalarT, Traits>::Constant(Teuchos::ParameterList& p) :
  value( p.get<double>("Value") ),
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);
  
  std::string n = "Constant Provider: " + constant.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename ScalarT, typename Traits>
Constant<ScalarT, Traits>::~Constant()
{ }

//**********************************************************************
template<typename ScalarT, typename Traits>
void Constant<ScalarT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  using namespace PHX;
  vm.template setFieldData(constant);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
}

//**********************************************************************
template<typename ScalarT, typename Traits>
void Constant<ScalarT, Traits>::evaluateFields(typename Traits::EvalData d)
{ }

//**********************************************************************
