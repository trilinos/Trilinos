// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
void Constant<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData /* d */,
		      PHX::FieldManager<Traits>& /* vm */)
{
  num_points = static_cast<int>(constant.extent(0));
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Constant<EvalT, Traits>::evaluateFields(typename Traits::EvalData /* d */)
{
  constant.deep_copy(value);
}

//**********************************************************************
// Needed for task parallel interface
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void Constant<EvalT, Traits>::operator () (const int cell) const
{
  for (int pt=0; pt < num_points; ++pt)
    constant(cell,pt) = value; 
}

//**********************************************************************
