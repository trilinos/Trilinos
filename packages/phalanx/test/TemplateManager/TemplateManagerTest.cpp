// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
