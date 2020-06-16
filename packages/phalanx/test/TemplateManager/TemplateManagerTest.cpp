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

#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_TemplateManager.hpp"
#include "MyTraits.hpp"
#include "Sacado_mpl_placeholders.hpp"

class DummyBase {
public:
  virtual ~DummyBase() = default;
  int val_;
};

template<typename T>
class DummyDerived : public DummyBase {
  T tmp_;
};

template<typename Traits>
class DummyTemplateManager :
  public PHX::TemplateManager<typename Traits::EvalTypes,
                              DummyBase,
                              DummyDerived<Sacado::mpl::placeholders::_> > {
public:
  DummyTemplateManager() = default;
  ~DummyTemplateManager() = default;
};

TEUCHOS_UNIT_TEST(template_manager_test, basic)
{
  using namespace PHX;
  using R = MyTraits::Residual;
  using J = MyTraits::Jacobian;

  DummyTemplateManager<MyTraits> tm;

  tm.disableType<J>();

  tm.buildObjects();

  auto r = tm.getAsObject<R>();
  auto j = tm.getAsObject<J>();

  // Check ability to disable type
  TEST_ASSERT(nonnull(r));
  TEST_ASSERT(j.is_null());

  // Check ability to delete type
  tm.deleteType<R>();
  auto r2 = tm.getAsObject<R>();
  TEST_ASSERT(r2.is_null());
}
