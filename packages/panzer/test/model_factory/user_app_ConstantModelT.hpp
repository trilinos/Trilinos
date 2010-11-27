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

#endif
