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
  void MockDAG<EvalT,Traits>::requires(const std::string& n,
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
