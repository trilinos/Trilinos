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
