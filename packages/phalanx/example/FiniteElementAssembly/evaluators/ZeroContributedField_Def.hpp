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
ZeroContributedField<EvalT,Traits>::
ZeroContributedField(const std::string& field_name,
                     const Teuchos::RCP<PHX::DataLayout>& layout) :
  field(field_name,layout)
{
  // "Evalauted" is always called before "Contributed" for the same field
  this->addEvaluatedField(field);
  this->setName("ZeroContributedField: " + field.fieldTag().identifier());
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
ZeroContributedField<EvalT,Traits>::evaluateFields(typename Traits::EvalData )
{
  field.deep_copy(ScalarT(0.0));
}

//**********************************************************************
