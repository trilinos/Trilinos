// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
  void EvaluatorWithMacros1<EvalT,Traits>::requires(const std::string& n)
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
  bindField(const PHX::FieldTag& , const PHX::any& )
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
  void EvaluatorWithMacros2<EvalT,Traits>::requires(const std::string& n)
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
  bindField(const PHX::FieldTag& , const PHX::any& )
  {
    // DO NOTHING! This overrides the requirement for a pointer to
    // fields for unit testing.
  }

} 
