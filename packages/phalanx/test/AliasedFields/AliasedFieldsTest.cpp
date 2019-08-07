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

#include <limits>

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDField_UnmanagedAllocator.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_TypeStrings.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "MyTraits.hpp"

PHX_EXTENT(CELL)
PHX_EXTENT(BASIS)

#include "AliasedFields_TestEvaluators.hpp"

// ****************************************
// Check that alias field works correctly
// ****************************************
TEUCHOS_UNIT_TEST(aliased_fields, basic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Put the field manager outside all evalaution code so that we can
  // test that memory management is handled correctly.
  FieldManager<MyTraits> fm;

  RCP<MDALayout<CELL,BASIS>> dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

  // Register evaluators
  {
    RCP<EvalA<MyTraits::Residual,MyTraits>> a = rcp(new EvalA<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(a);
  }
  {
    RCP<EvalC<MyTraits::Residual,MyTraits>> c = rcp(new EvalC<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(c);
  }

  // Require field a
  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag_a("a",dl);
    fm.requireField<MyTraits::Residual>(tag_a);
  }
  
  // Alias "b" to "c"
  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag_b("b",dl);
    PHX::Tag<MyTraits::Residual::ScalarT> tag_c("c",dl);
    fm.aliasFieldForAllEvaluationTypes(tag_b,tag_c);
  }
  
  fm.postRegistrationSetup(0);
  fm.evaluateFields<MyTraits::Residual>(0);
  
  MDField<double,CELL,BASIS> b("b",dl);
  MDField<double,CELL,BASIS> c("c",dl);
  fm.getFieldData<MyTraits::Residual,double,CELL,BASIS>(b);
  fm.getFieldData<MyTraits::Residual,double,CELL,BASIS>(c);

  const auto tol = std::numeric_limits<double>::epsilon() * 100.0;

  for (int cell = 0; cell < b.extent_int(0); ++cell) {
    for (int basis = 0; basis < c.extent_int(1); ++basis) {
      TEST_FLOATING_EQUALITY(b(cell,basis),4.0,tol);
    }
  }
}

// ****************************************
// Check that we catch bad data layouts in the aliased field
// ****************************************
TEUCHOS_UNIT_TEST(aliased_fields, throw_on_bad_layout)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Put the field manager outside all evalaution code so that we can
  // test that memory management is handled correctly.
  FieldManager<MyTraits> fm;

  RCP<MDALayout<CELL,BASIS>> dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));
  RCP<MDALayout<CELL,BASIS>> bad_dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",50,4));

  // Register evaluators
  {
    RCP<EvalA<MyTraits::Residual,MyTraits>> a = rcp(new EvalA<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(a);
  }
  {
    RCP<EvalC<MyTraits::Residual,MyTraits>> c = rcp(new EvalC<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(c);
  }

  // Require field a
  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag_a("a",dl);
    fm.requireField<MyTraits::Residual>(tag_a);
  }
  
  // Alias "b" to "c"
  {
    // use an incorrectly sized data layout
    PHX::Tag<MyTraits::Residual::ScalarT> tag_b("b",bad_dl);
    PHX::Tag<MyTraits::Residual::ScalarT> tag_c("c",dl);
    TEST_THROW(fm.aliasFieldForAllEvaluationTypes(tag_b,tag_c),std::runtime_error)
  }
  
}

// ****************************************
// Check that we catch bad scalar types in the aliased field
// ****************************************
TEUCHOS_UNIT_TEST(aliased_fields, throw_on_bad_scalar_type)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Put the field manager outside all evalaution code so that we can
  // test that memory management is handled correctly.
  FieldManager<MyTraits> fm;

  RCP<MDALayout<CELL,BASIS>> dl = rcp(new PHX::MDALayout<CELL,BASIS>("H-Grad",100,4));

  // Register evaluators
  {
    RCP<EvalA<MyTraits::Residual,MyTraits>> a = rcp(new EvalA<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(a);
  }
  {
    RCP<EvalC<MyTraits::Residual,MyTraits>> c = rcp(new EvalC<MyTraits::Residual,MyTraits>(dl));
    fm.registerEvaluator<MyTraits::Residual>(c);
  }

  // Require field a
  {
    PHX::Tag<MyTraits::Residual::ScalarT> tag_a("a",dl);
    fm.requireField<MyTraits::Residual>(tag_a);
  }
  
  // Alias "b" to "c"
  {
    // inject a float instead of double
    PHX::Tag<float> tag_b("b",dl);
    PHX::Tag<MyTraits::Residual::ScalarT> tag_c("c",dl);
    TEST_THROW(fm.aliasFieldForAllEvaluationTypes(tag_b,tag_c),std::runtime_error)
  }
  
}
