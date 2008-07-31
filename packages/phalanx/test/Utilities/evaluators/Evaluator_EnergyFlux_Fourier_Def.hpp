//**********************************************************************
template<typename EvalT, typename Traits> 
Fourier<EvalT, Traits>::
Fourier(const Teuchos::ParameterList& p) :
  flux("Energy_Flux", 
       p.get< Teuchos::RCP<PHX::DataLayout> >("Vector Data Layout") ),
  density("Density", 
	  p.get< Teuchos::RCP<PHX::DataLayout> >("Scalar Data Layout") ),
  dc("Diffusion Coefficient", 
     p.get< Teuchos::RCP<PHX::DataLayout> >("Scalar Data Layout") ),
  grad_temp("Temperature Gradient", 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Vector Data Layout") )
{ 
  this->addEvaluatedField(flux);
  this->addDependentField(density);
  this->addDependentField(dc);
  this->addDependentField(grad_temp);

  this->setName("Fourier");
}

//**********************************************************************
template<typename EvalT, typename Traits> 
Fourier<EvalT, Traits>::~Fourier()
{ }

//**********************************************************************
template<typename EvalT, typename Traits> 
void Fourier<EvalT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  this->utils.setFieldData(flux,vm);
  this->utils.setFieldData(density,vm);
  this->utils.setFieldData(dc,vm);
  this->utils.setFieldData(grad_temp,vm);

  data_layout_size = flux.fieldTag().dataLayout()->size();
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Fourier<EvalT, Traits>::evaluateFields(typename Traits::EvalData d)
{ 
  std::size_t size = d.size() * data_layout_size;

  for (std::size_t i = 0; i < size; ++i)
    flux[i] = - density[i] * dc[i] * grad_temp[i];
}

//**********************************************************************
