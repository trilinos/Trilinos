// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>

// Teko-Package includes
#include "Teko_Utilities.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_SpmdVectorBase.hpp"

// Test-rig

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;

TEUCHOS_UNIT_TEST(tUtilitiesTests, clipping) {
  RCP<Thyra::VectorSpaceBase<double> > vs = Thyra::defaultSpmdVectorSpace<double>(10);
  RCP<Thyra::MultiVectorBase<double> > x  = Thyra::createMembers<double>(vs, 2);

  Thyra::assign(x.ptr(), 1.0);

  // try to clip lower values
  Teko::clipLower(x, 2.0);

  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 0), 10.0 * 2.0, 1e-16);
  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 1), 10.0 * 2.0, 1e-16);

  Teko::clipUpper(x, 1.0);

  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 0), 10.0 * 1.0, 1e-16);
  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 1), 10.0 * 1.0, 1e-16);

  Teuchos::ArrayRCP<double> col0, col1;
  rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(0))->getNonconstLocalData(
      Teuchos::ptrFromRef(col0));
  rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(1))->getNonconstLocalData(
      Teuchos::ptrFromRef(col1));

  TEST_EQUALITY(col0.size(), 10);
  TEST_EQUALITY(col1.size(), 10);

  for (int i = 0; i < 10; i++) {
    col0[i] = i - 5;
    col1[i] = 5 - i;
  }

  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 0), 25.0, 1e-16);
  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 1), 25.0, 1e-16);

  Teko::clipLower(x, 0.0);

  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 0), 10.0, 1e-16);
  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 1), 15.0, 1e-16);

  for (int i = 0; i < 10; i++) {
    col0[i] = i - 5;
    col1[i] = 5 - i;
  }

  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 0), 25.0, 1e-16);
  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 1), 25.0, 1e-16);

  Teko::clipUpper(x, 0.0);

  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 0), 15.0, 1e-16);
  TEST_FLOATING_EQUALITY(Teko::norm_1(x, 1), 10.0, 1e-16);
}

TEUCHOS_UNIT_TEST(tUtilitiesTests, replaceValues) {
  RCP<Thyra::VectorSpaceBase<double> > vs = Thyra::defaultSpmdVectorSpace<double>(10);
  RCP<Thyra::MultiVectorBase<double> > x  = Thyra::createMembers<double>(vs, 2);

  Teuchos::ArrayRCP<double> col0, col1;
  rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(0))->getNonconstLocalData(
      Teuchos::ptrFromRef(col0));
  rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(1))->getNonconstLocalData(
      Teuchos::ptrFromRef(col1));

  for (int i = 0; i < 10; i++) {
    col0[i] = i;
    col1[i] = -i;
  }

  Teko::replaceValue(x, 0.0, 99.0);

  TEST_EQUALITY(col0[0], 99.0);
  TEST_EQUALITY(col1[0], 99.0);
  for (int i = 1; i < 10; i++) {
    TEST_EQUALITY(col0[i], double(i));
    TEST_EQUALITY(col1[i], double(-i));
  }
}

TEUCHOS_UNIT_TEST(tUtilitiesTests, averages) {
  RCP<Thyra::VectorSpaceBase<double> > vs = Thyra::defaultSpmdVectorSpace<double>(9);
  RCP<Thyra::MultiVectorBase<double> > x  = Thyra::createMembers<double>(vs, 2);

  Teuchos::ArrayRCP<double> col0, col1;
  rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(0))->getNonconstLocalData(
      Teuchos::ptrFromRef(col0));
  rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(x->col(1))->getNonconstLocalData(
      Teuchos::ptrFromRef(col1));

  for (int i = 0; i < 9; i++) {
    col0[i] = 2.3 + (i - 4.0);
    col1[i] = -3.7 + (i - 4.0);
  }

  std::vector<double> averages;
  Teko::columnAverages(x, averages);

  TEST_EQUALITY(averages.size(), 2);
  TEST_FLOATING_EQUALITY(averages[0], 2.3, 1e-14);
  TEST_FLOATING_EQUALITY(averages[1], -3.7, 1e-14);

  double avg = Teko::average(x);
  TEST_FLOATING_EQUALITY(avg, (9.0 * 2.3 - 9.0 * 3.7) / 18.0, 1e-14);
}
