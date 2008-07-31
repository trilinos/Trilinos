#ifndef PHX_EXAMPLE_VP_NONLINEAR_SOURCE_DEF_HPP
#define PHX_EXAMPLE_VP_NONLINEAR_SOURCE_DEF_HPP

//**********************************************************************
template<typename EvalT, typename Traits> NonlinearSource<EvalT, Traits>::
NonlinearSource(const Teuchos::ParameterList& p) :
  source("Nonlinear Source", 
	 p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout")),
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout")),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout"))
{ 
  this->addEvaluatedField(source);
  this->addDependentField(density);
  this->addDependentField(temp);

  this->setName("NonlinearSource");
}

//**********************************************************************
template<typename EvalT, typename Traits> 
NonlinearSource<EvalT, Traits>::~NonlinearSource()
{ }

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
postRegistrationSetup(PHX::FieldManager<Traits>& vm)
{
  this->utils.setFieldData(source,vm);
  this->utils.setFieldData(density,vm);
  this->utils.setFieldData(temp,vm);

  data_layout_size = source.fieldTag().dataLayout()->size();
}

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
evaluateFields(typename Traits::EvalData d)
{ 
  std::size_t size = d.size() * data_layout_size;
  
  for (std::size_t i = 0; i < size; ++i)
    source[i] = density[i] * temp[i] * temp[i];
}

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
preEvaluate(typename Traits::PreEvalData d)
{ 
  using namespace std;
  cout << "In Density Pre Op" << endl;
}

//**********************************************************************
template< typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
postEvaluate(typename Traits::PostEvalData d)
{ 
  using namespace std;
  cout << "In Density Post Op" << endl;
}

//**********************************************************************

#endif
