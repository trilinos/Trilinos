// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CONSTANT_IMPL_HPP
#define PANZER_CONSTANT_IMPL_HPP

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Constant<EvalT, Traits>::
Constant(
  const Teuchos::ParameterList& p) :
  value( p.get<double>("Value") ),
  constant( p.get<std::string>("Name"), 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout") )
{
  this->addEvaluatedField(constant);

  // Make this unshared so that it is not overwritten
  this->addUnsharedField(constant.fieldTag().clone());

  std::string n = "Constant: " + constant.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Constant<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData  /* worksets */,
  PHX::FieldManager<Traits>&  fm)
{
  using namespace PHX;
  this->utils.setFieldData(constant,fm);
  
  constant.deep_copy(value);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Constant<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  /* d */)
{}

//**********************************************************************

}

#endif
