#ifndef USER_APP_CONSTANT_MODEL_T_HPP
#define USER_APP_CONSTANT_MODEL_T_HPP

//**********************************************************************
PHX_EVALUATOR_CTOR_NAMESPACE(user_app,ConstantModel,p) :
  value( p.get<double>("Value") ),
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);
  
  std::string n = "Constant: " + constant.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(user_app::ConstantModel,worksets,fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
}

//**********************************************************************
PHX_EVALUATE_FIELDS(user_app::ConstantModel,d)
{ }

//**********************************************************************


#ifdef HAVE_STOKHOS

namespace user_app {

template <typename Traits>
ConstantModel<typename Traits::SGResidual,Traits>::ConstantModel(const Teuchos::ParameterList & p) :
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);

  if(p.isParameter("Expansion")) {
     expansion = p.get<Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > >("Expansion");
     value.reset(expansion);
     value.copyForWrite();
     value.fastAccessCoeff(0) = p.get<double>("Value");
     value.fastAccessCoeff(1) = p.get<double>("UQ");
  }
  else {
     value = p.get<double>("Value");
  }

  std::string n = "Constant: " + constant.fieldTag().name();
  this->setName(n);
}

template <typename Traits>
void ConstantModel<typename Traits::SGResidual,Traits>::postRegistrationSetup(typename Traits::SetupData d,
                                                                              PHX::FieldManager<Traits>& fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
}

template <typename Traits>
ConstantModel<typename Traits::SGJacobian,Traits>::ConstantModel(const Teuchos::ParameterList & p) :
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);

  if(p.isParameter("Expansion")) {
     expansion = p.get<Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > >("Expansion");
     value.val().reset(expansion);
     value.val().copyForWrite();
     value.val().fastAccessCoeff(0) = p.get<double>("Value");
     value.val().fastAccessCoeff(1) = p.get<double>("UQ");
  }
  else {
     value.val() = p.get<double>("Value");
  }

  std::string n = "Constant: " + constant.fieldTag().name();
  this->setName(n);
}

template <typename Traits>
void ConstantModel<typename Traits::SGJacobian,Traits>::postRegistrationSetup(typename Traits::SetupData d,
                                                                              PHX::FieldManager<Traits>& fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
}

}

#endif

#endif
