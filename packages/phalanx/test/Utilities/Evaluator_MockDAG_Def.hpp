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
#include "Phalanx_DataLayout_DynamicLayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_MDField.hpp"

class CELL;
class BASIS;

namespace PHX {

  template<typename EvalT,typename Traits>
  void MockDAG<EvalT,Traits>::
  postRegistrationSetup(typename Traits::SetupData ,
			PHX::FieldManager<Traits>& )
  {
    //this->utils.setFieldData(flux,fm);
  }

  template<typename EvalT,typename Traits>
  void MockDAG<EvalT,Traits>::evaluateFields(typename Traits::EvalData /* d */){}

  template<typename EvalT,typename Traits>
  void MockDAG<EvalT,Traits>::evaluates(const std::string& n,
                                        const bool use_dynamic_layout)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::DataLayout> dl;
    if (use_dynamic_layout)
      dl = rcp(new PHX::Layout("H-Grad"));
    else
      dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

    PHX::Tag<typename EvalT::ScalarT> tag(n,dl);
    this->addEvaluatedField(tag);
  }

  template<typename EvalT,typename Traits>
  void MockDAG<EvalT,Traits>::depends(const std::string& n,
                                      const bool use_dynamic_layout)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::DataLayout> dl;
    if (use_dynamic_layout)
      dl = rcp(new PHX::Layout("H-Grad"));
    else
      dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

    PHX::Tag<typename EvalT::ScalarT> tag(n,dl);
    this->addDependentField(tag);
  }

  template<typename EvalT,typename Traits>
  void MockDAG<EvalT,Traits>::contributes(const std::string& n,
                                          const bool use_dynamic_layout)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::DataLayout> dl;
    if (use_dynamic_layout)
      dl = rcp(new PHX::Layout("H-Grad"));
    else
      dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

    PHX::Tag<typename EvalT::ScalarT> tag(n,dl);
    this->addContributedField(tag);
  }

  template<typename EvalT,typename Traits>
  void MockDAG<EvalT,Traits>::unshared(const std::string& n)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<PHX::MDALayout<CELL,BASIS>> dl =
      rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
    RCP<PHX::Tag<typename EvalT::ScalarT>> tag =
      Teuchos::rcp(new PHX::Tag<typename EvalT::ScalarT>(n,dl));
    this->addUnsharedField(tag);
  }

}
