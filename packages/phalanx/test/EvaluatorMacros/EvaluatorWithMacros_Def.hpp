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

namespace PHX {

  // Implementation #1
  PHX_EVALUATOR_CTOR(EvaluatorWithMacros1,plist)
  {
    std::cout << "CTOR1: " << plist.name() << std::endl;
  }

  PHX_POST_REGISTRATION_SETUP(EvaluatorWithMacros1,data,fm)
  {
    std::vector<PHX::index_size_type> ddims(1,0);
    fm.template setKokkosExtendedDataTypeDimensions<EvalT>(ddims);
    std::cout << "postRegistrationSetup1: " << data << std::endl;
  }

  PHX_EVALUATE_FIELDS(EvaluatorWithMacros1,data)
  {
    std::cout << "evaluateFields1: " << data << std::endl;
  }

  template<typename EvalT,typename Traits>
  void EvaluatorWithMacros1<EvalT,Traits>::evaluates(const std::string& n)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::MDALayout<CELL,BASIS>> dl = 
      rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
    PHX::Tag<typename EvalT::ScalarT> tag(n,dl);
    this->addEvaluatedField(tag);
  }

  template<typename EvalT,typename Traits>
  void EvaluatorWithMacros1<EvalT,Traits>::depends(const std::string& n)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::MDALayout<CELL,BASIS>> dl = 
      rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
    PHX::Tag<typename EvalT::ScalarT> tag(n,dl);
    this->addDependentField(tag);
  }

  template<typename EvalT,typename Traits>
  void EvaluatorWithMacros1<EvalT,Traits>::
  bindField(const PHX::FieldTag& , const std::any& )
  {
    // DO NOTHING! This overrides the requirement for a pointer to
    // fields for unit testing.
  }
  
  // Implementation #2
  PHX_EVALUATOR_CTOR(EvaluatorWithMacros2,plist)
  {
    std::cout << "CTOR2: " << plist.name() << std::endl;
  }

  PHX_POST_REGISTRATION_SETUP(EvaluatorWithMacros2,data,fm)
  {
    std::vector<PHX::index_size_type> ddims(1,0);
    fm.template setKokkosExtendedDataTypeDimensions<EvalT>(ddims);
    std::cout << "postRegistrationSetup2: " << data << std::endl;
  }

  PHX_EVALUATE_FIELDS(EvaluatorWithMacros2,data)
  {
    std::cout << "evaluateFields2: " << data << std::endl;
  }

  PHX_PRE_EVALUATE_FIELDS(EvaluatorWithMacros2,data)
  {
    std::cout << "preEvaluateFields2: " << data << std::endl;
  }

  PHX_POST_EVALUATE_FIELDS(EvaluatorWithMacros2,data)
  {
    std::cout << "postEvaluateFields2: " << data << std::endl;
  }

  template<typename EvalT,typename Traits>
  void EvaluatorWithMacros2<EvalT,Traits>::evaluates(const std::string& n)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::MDALayout<CELL,BASIS>> dl = 
      rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
    PHX::Tag<typename EvalT::ScalarT> tag(n,dl);
    this->addEvaluatedField(tag);
  }

  template<typename EvalT,typename Traits>
  void EvaluatorWithMacros2<EvalT,Traits>::depends(const std::string& n)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::MDALayout<CELL,BASIS>> dl = 
      rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
    PHX::Tag<typename EvalT::ScalarT> tag(n,dl);
    this->addDependentField(tag);
  }

  template<typename EvalT,typename Traits>
  void EvaluatorWithMacros2<EvalT,Traits>::
  bindField(const PHX::FieldTag& , const std::any& )
  {
    // DO NOTHING! This overrides the requirement for a pointer to
    // fields for unit testing.
  }

} 
