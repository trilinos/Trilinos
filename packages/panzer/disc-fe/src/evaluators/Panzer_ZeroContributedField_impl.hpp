// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   PANZER_ZEROCONTRIBUTEDFIELD_IMPL_HPP
#define   PANZER_ZEROCONTRIBUTEDFIELD_IMPL_HPP

namespace panzer
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  ZeroContributedField<EvalT, Traits>::
  ZeroContributedField(
    const std::string& fieldName,
    PHX::DataLayout&   layout)
  {
    using PHX::MDField;
    using Teuchos::rcpFromRef;
    field_ = MDField<ScalarT>(fieldName, rcpFromRef(layout));
    this->addEvaluatedField(field_);
    this->setName("ZeroContributedField:  " + field_.fieldTag().identifier());
  } // end of Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  evaluateFields()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  ZeroContributedField<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData /* d */)
  {
    field_.deep_copy(ScalarT(0.0));
  } // end of evaluateFields()

} // end of namespace panzer

#endif // PANZER_ZEROCONTRIBUTEDFIELD_IMPL_HPP
