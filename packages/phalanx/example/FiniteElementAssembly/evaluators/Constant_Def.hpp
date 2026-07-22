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
Constant<EvalT,Traits>::
Constant(const std::string& field_name,
         const Teuchos::RCP<PHX::DataLayout>& layout,
         const double& val) :
  value(val),
  constant(field_name,layout)
{
  this->addEvaluatedField(constant);
  std::string n = "Constant: " + constant.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Constant<EvalT,Traits>::
postRegistrationSetup(typename Traits::SetupData ,
                      PHX::FieldManager<Traits>& )
{ constant.deep_copy(value); }

//**********************************************************************
template<typename EvalT, typename Traits>
void Constant<EvalT,Traits>::evaluateFields(typename Traits::EvalData )
{ }

//**********************************************************************
