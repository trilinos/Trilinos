// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef USER_APP_CONSTANT_MODEL_T_HPP
#define USER_APP_CONSTANT_MODEL_T_HPP

//**********************************************************************
template<typename EvalT, typename Traits>
user_app::ConstantModel<EvalT, Traits>::
ConstantModel(
  const Teuchos::ParameterList& p) :
  value( p.get<double>("Value") ),
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);

  // Make this unshared so that it is not overwritten
  this->addUnsharedField(constant.fieldTag().clone());
  
  std::string n = "user_app::Constant: " + constant.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
user_app::ConstantModel<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData worksets,
  PHX::FieldManager<Traits>& fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);

  for (std::size_t i = 0; i < static_cast<std::size_t>(constant.size()); ++i)
    constant[i] = value;
  TEUCHOS_ASSERT(false);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
user_app::ConstantModel<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData d)
{}

//**********************************************************************

#endif
