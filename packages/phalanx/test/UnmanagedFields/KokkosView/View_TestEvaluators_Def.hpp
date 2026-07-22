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

  PHX_EVALUATOR_CTOR(EvalUnmanaged,plist) :
    tag_a("a",plist.get<Teuchos::RCP<PHX::DataLayout>>("dl")),
    tag_b("b",plist.get<Teuchos::RCP<PHX::DataLayout>>("dl")),
    tag_c("c",plist.get<Teuchos::RCP<PHX::DataLayout>>("dl")),
    tag_d("d",plist.get<Teuchos::RCP<PHX::DataLayout>>("dl"))  
  {
    // static
    this->addEvaluatedField(tag_a,a);
    this->addDependentField(tag_b,b);

    // dynamic
    this->addEvaluatedField(tag_c,c);
    this->addDependentField(tag_d,d);
  }

  PHX_POST_REGISTRATION_SETUP(EvalUnmanaged,/* data */,fm)
  {
    this->utils.setFieldData(tag_a,a,fm);
    this->utils.setFieldData(tag_b,b,fm);
    this->utils.setFieldData(tag_c,c,fm);
    this->utils.setFieldData(tag_d,d,fm);
  }

  PHX_EVALUATE_FIELDS(EvalUnmanaged,/* data */)
  {
    Kokkos::deep_copy(a,b);
    Kokkos::deep_copy(c,d);
  }

  PHX_EVALUATOR_CTOR(EvalDummy,plist) :
      tag_b("b",plist.get<Teuchos::RCP<PHX::DataLayout>>("dl")),
      tag_d("d",plist.get<Teuchos::RCP<PHX::DataLayout>>("dl"))  
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    // static
    this->addEvaluatedField(tag_b,b);

    // dynamic
    this->addEvaluatedField(tag_d,d);
  }

  PHX_POST_REGISTRATION_SETUP(EvalDummy,/* data */,fm)
  {
    this->utils.setFieldData(tag_b,b,fm);
    this->utils.setFieldData(tag_d,d,fm);
  }

  PHX_EVALUATE_FIELDS(EvalDummy,/* data */)
  {
    // do nothing - they are unmanaged - data values are set by user
  }

} 
