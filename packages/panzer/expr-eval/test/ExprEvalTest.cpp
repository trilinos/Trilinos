// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include "Panzer_ExprEval_impl.hpp"

namespace panzer {

TEUCHOS_UNIT_TEST(ExprEval, test_1d_system)
{
  Expr::Eval<double*> eval;
  Expr::set_cmath_functions(eval);
  Teuchos::any result;
  eval.read_string(result, "1+1", "one plus one");
  Kokkos::View<double const> x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  auto h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 2.0);
  eval.set("a", 5.0);
  eval.read_string(result, "a+a", "a plus a");
  x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 10.0);
  eval.read_string(result, "b=4;\na+b", "a plus b");
  x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 9.0);
  eval.read_string(result, "a-b", "a minus b");
  x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 1.0);
  eval.read_string(result, "abs(b-a)", "abs(b-a)");
  x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 1.0);
  eval.read_string(result, "a^2 - b^2", "a^2 - b^2");
  x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 9.0);
  auto y = Kokkos::View<double*>("y", 3);
  auto h_y = Kokkos::create_mirror_view(y);
  h_y(0) = 1.0;
  h_y(1) = 2.0;
  h_y(2) = 3.0;
  Kokkos::deep_copy(y, h_y);
  eval.set("y", y);
  eval.read_string(result, "2 + y", "2 + y");
  auto z = Teuchos::any_cast<Kokkos::View<double const*>>(result);
  auto h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0) == 3.0);
  TEUCHOS_ASSERT(h_z(1) == 4.0);
  TEUCHOS_ASSERT(h_z(2) == 5.0);
  eval.read_string(result, "y * a", "y * a");
  z = Teuchos::any_cast<Kokkos::View<double const*>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0) == 5.0);
  TEUCHOS_ASSERT(h_z(1) == 10.0);
  TEUCHOS_ASSERT(h_z(2) == 15.0);
  eval.read_string(result, "y + y^2", "y + y^2");
  z = Teuchos::any_cast<Kokkos::View<double const*>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0) == 2.0);
  TEUCHOS_ASSERT(h_z(1) == 6.0);
  TEUCHOS_ASSERT(h_z(2) == 12.0);
  eval.read_string(result, "0.5 < 1.0 ? 0.5 : 1.0", "0.5 < 1.0 ? 0.5 : 1.0");
  x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 0.5);
  eval.read_string(result, "y < 1.5 ? 0.0 : 1.0", "y < 1.5 ? 0.0 : 1.0");
  z = Teuchos::any_cast<Kokkos::View<double const*>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0) == 0.0);
  TEUCHOS_ASSERT(h_z(1) == 1.0);
  TEUCHOS_ASSERT(h_z(2) == 1.0);
  eval.read_string(result, "(-y) + (-a)", "(-y) + (-a)");
  z = Teuchos::any_cast<Kokkos::View<double const*>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0) == -6.0);
  TEUCHOS_ASSERT(h_z(1) == -7.0);
  TEUCHOS_ASSERT(h_z(2) == -8.0);
  eval.read_string(result, "abs(y - a)", "abs(y - a)");
  z = Teuchos::any_cast<Kokkos::View<double const*>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0) == 4.0);
  TEUCHOS_ASSERT(h_z(1) == 3.0);
  TEUCHOS_ASSERT(h_z(2) == 2.0);
}

TEUCHOS_UNIT_TEST(ExprEval, test_2d_system)
{
  Expr::Eval<double**> eval;
  Teuchos::any result;
  eval.set("a", 5.0);
  eval.read_string(result, "a+a", "a plus a");
  Kokkos::View<double const> x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  auto h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 10.0);
  eval.read_string(result, "b=4;\na+b", "a plus b");
  x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 9.0);
  auto y = Kokkos::View<double**>("y", 3, 2);
  auto h_y = Kokkos::create_mirror_view(y);
  h_y(0, 0) = h_y(0, 1) = 1.0;
  h_y(1, 0) = h_y(1, 1) = 2.0;
  h_y(2, 0) = h_y(2, 1) = 3.0;
  Kokkos::deep_copy(y, h_y);
  eval.set("y", y);
  eval.read_string(result, "2 + y", "2 + y");
  auto z = Teuchos::any_cast<Kokkos::View<double const**>>(result);
  auto h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0, 0) == 3.0);
  TEUCHOS_ASSERT(h_z(1, 0) == 4.0);
  TEUCHOS_ASSERT(h_z(2, 0) == 5.0);
  TEUCHOS_ASSERT(h_z(0, 1) == 3.0);
  TEUCHOS_ASSERT(h_z(1, 1) == 4.0);
  TEUCHOS_ASSERT(h_z(2, 1) == 5.0);
  eval.read_string(result, "y * a", "y * a");
  z = Teuchos::any_cast<Kokkos::View<double const**>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0, 0) == 5.0);
  TEUCHOS_ASSERT(h_z(1, 0) == 10.0);
  TEUCHOS_ASSERT(h_z(2, 0) == 15.0);
  TEUCHOS_ASSERT(h_z(0, 1) == 5.0);
  TEUCHOS_ASSERT(h_z(1, 1) == 10.0);
  TEUCHOS_ASSERT(h_z(2, 1) == 15.0);
  eval.read_string(result, "y + y^2", "y + y^2");
  z = Teuchos::any_cast<Kokkos::View<double const**>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0, 0) == 2.0);
  TEUCHOS_ASSERT(h_z(1, 0) == 6.0);
  TEUCHOS_ASSERT(h_z(2, 0) == 12.0);
  TEUCHOS_ASSERT(h_z(0, 1) == 2.0);
  TEUCHOS_ASSERT(h_z(1, 1) == 6.0);
  TEUCHOS_ASSERT(h_z(2, 1) == 12.0);
  eval.read_string(result, "0.5 < 1.0 ? 0.5 : 1.0", "0.5 < 1.0 ? 0.5 : 1.0");
  x = Teuchos::any_cast<Kokkos::View<double const>>(result);
  h_x = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(h_x, x);
  TEUCHOS_ASSERT(h_x() == 0.5);
  eval.read_string(result, "y < 1.5 ? 0.0 : 1.0", "y < 1.5 ? 0.0 : 1.0");
  z = Teuchos::any_cast<Kokkos::View<double const**>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0, 0) == 0.0);
  TEUCHOS_ASSERT(h_z(1, 0) == 1.0);
  TEUCHOS_ASSERT(h_z(2, 0) == 1.0);
  TEUCHOS_ASSERT(h_z(0, 1) == 0.0);
  TEUCHOS_ASSERT(h_z(1, 1) == 1.0);
  TEUCHOS_ASSERT(h_z(2, 1) == 1.0);
  eval.read_string(result, "(-y) + (-a)", "(-y) + (-a)");
  z = Teuchos::any_cast<Kokkos::View<double const**>>(result);
  h_z = Kokkos::create_mirror_view(z);
  Kokkos::deep_copy(h_z, z);
  TEUCHOS_ASSERT(h_z(0, 0) == -6.0);
  TEUCHOS_ASSERT(h_z(1, 0) == -7.0);
  TEUCHOS_ASSERT(h_z(2, 0) == -8.0);
  TEUCHOS_ASSERT(h_z(0, 1) == -6.0);
  TEUCHOS_ASSERT(h_z(1, 1) == -7.0);
  TEUCHOS_ASSERT(h_z(2, 1) == -8.0);
}

}
