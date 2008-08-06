//**********************************************************************
template<typename EvalT, typename Traits> Density<EvalT, Traits>::
Density(const Teuchos::ParameterList& p) :
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") ),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{ 
  this->addEvaluatedField(density);
  this->addDependentField(temp);
  this->setName("Density");
}

//**********************************************************************
template<typename EvalT, typename Traits>
Density<EvalT, Traits>::~Density()
{ }

//**********************************************************************
template<typename EvalT, typename Traits>
void Density<EvalT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  this->utils.setFieldData(density,vm);
  this->utils.setFieldData(temp,vm);

  data_layout_size = density.fieldTag().dataLayout().size();
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Density<EvalT, Traits>::evaluateFields(typename Traits::EvalData d)
{ 
  std::size_t size = d.size() * data_layout_size;
  
  for (std::size_t i = 0; i < size; ++i)
    density[i] =  temp[i] * temp[i];
}

//**********************************************************************
