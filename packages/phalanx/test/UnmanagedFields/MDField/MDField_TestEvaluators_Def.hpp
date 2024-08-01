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
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

class CELL;
class BASIS;

namespace PHX {

  PHX_EVALUATOR_CTOR(EvalUnmanaged,/* plist */)

  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::MDALayout<CELL,BASIS>> dl = 
      rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

    a = PHX::MDField<double,CELL,BASIS>("a",dl);
    b = PHX::MDField<const double,CELL,BASIS>("b",dl);
    c = PHX::MDField<double>("c",dl);
    d = PHX::MDField<const double>("d",dl);

    // static
    this->addEvaluatedField(a);
    this->addDependentField(b);

    // dynamic
    this->addEvaluatedField(c);
    this->addDependentField(d);
  }

  PHX_POST_REGISTRATION_SETUP(EvalUnmanaged,/* data */,fm)
  {
    this->utils.setFieldData(a,fm);
    this->utils.setFieldData(b,fm);
    this->utils.setFieldData(c,fm);
    this->utils.setFieldData(d,fm);
  }

  PHX_EVALUATE_FIELDS(EvalUnmanaged,/* data */)
  {
    a.deep_copy(b);
    c.deep_copy(d);
  }

  PHX_EVALUATOR_CTOR(EvalDummy,/* plist */)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::MDALayout<CELL,BASIS>> dl = 
      rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

    b = PHX::MDField<double,CELL,BASIS>("b",dl);
    d = PHX::MDField<double>("d",dl);

    // static
    this->addEvaluatedField(b);

    // dynamic
    this->addEvaluatedField(d);
  }

  PHX_POST_REGISTRATION_SETUP(EvalDummy,/* data */,fm)
  {
    this->utils.setFieldData(b,fm);
    this->utils.setFieldData(d,fm);
  }

  PHX_EVALUATE_FIELDS(EvalDummy,/* data */)
  {
    // do nothing - they are unmanaged - data values are set by user
  }

} 
