//**********************************************************************
template<typename ScalarT, typename Traits> Density<ScalarT, Traits>::
Density(const Teuchos::ParameterList& p) :
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{ 
  this->addEvaluatedField(density);
  this->addDependentField(temp);
  this->setName("Density");
}

//**********************************************************************
template<typename ScalarT, typename Traits>
Density<ScalarT, Traits>::~Density()
{ }

//**********************************************************************
template<typename ScalarT, typename Traits>
void Density<ScalarT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  vm.setFieldData(density);
  vm.setFieldData(temp);

  data_layout_size = density.fieldTag().dataLayout()->size();
}

//**********************************************************************
template<typename ScalarT, typename Traits>
void Density<ScalarT, Traits>::evaluateFields(typename Traits::EvalData d)
{ 
  std::size_t size = d.size() * data_layout_size;
  
  for (std::size_t i = 0; i < size; ++i)
    density[i] =  temp[i] * temp[i];
}

//**********************************************************************
