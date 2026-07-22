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
#include "Phalanx_DataLayout_DynamicLayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"

struct CELL;
struct BASIS;

namespace PHX {

  PHX_EVALUATOR_CTOR(EvalUnmanaged,/* plist */)

  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::DataLayout> dl = 
      rcp(new PHX::Layout("H-Grad<CELL,BASIS>",100,4));

    a = PHX::Field<double,2>("a",dl);
    b = PHX::Field<const double,2>("b",dl);
    c = PHX::Field<double,2>("c",dl);
    d = PHX::Field<const double,2>("d",dl);

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

    RCP<PHX::DataLayout> dl = 
      rcp(new PHX::Layout("H-Grad<CELL,BASIS>",100,4));

    b = PHX::Field<double,2>("b",dl);
    d = PHX::Field<double,2>("d",dl);

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
