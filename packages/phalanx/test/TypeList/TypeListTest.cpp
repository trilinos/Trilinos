// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// The sacado_vector_compatibility test below intentionally uses the
// deprecated Sacado::mpl::vector path. The warning is reported inside
// Phalanx_TypeList.hpp, so it must be disabled before the include.
#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_TypeList.hpp"
#include "Sacado_mpl_find.hpp"
#include "Sacado_mpl_for_each.hpp"
#include <string>
#include <vector>

namespace {

  struct Residual {};
  struct Jacobian {};
  struct Tangent {};
  struct Stochastic {}; // never in EvalTypes

  using namespace PHX::mpl;
  using EvalTypes = type_list<Residual, Jacobian, Tangent>;

  // size, at, contains, find
  static_assert(size_v<EvalTypes> == 3);
  static_assert(size_v<type_list<>> == 0);
  static_assert(std::is_same_v<at_t<EvalTypes,1>, Jacobian>);
  static_assert(contains_v<Jacobian, EvalTypes>);
  static_assert(!contains_v<Stochastic, EvalTypes>);
  static_assert(find_v<Residual, EvalTypes> == 0);
  static_assert(find_v<Tangent, EvalTypes> == 2);
  static_assert(found_v<Jacobian, EvalTypes>);
  static_assert(!found_v<Stochastic, EvalTypes>);
  static_assert(find_v<Stochastic, EvalTypes> == npos<EvalTypes>);
  static_assert(find_v<Residual, type_list<>> == npos<type_list<>>);

  // index_v rejects types not in the list at compile time
  template <typename T, typename L>
  concept HasIndex = requires { index_v<T, L>; };
  static_assert(index_v<Tangent, EvalTypes> == 2);
  static_assert(HasIndex<Jacobian, EvalTypes>);
  static_assert(!HasIndex<Stochastic, EvalTypes>);

  // find_if
  using NumericTypes = type_list<int, float, double>;
  static_assert(find_if_v<std::is_floating_point, NumericTypes> == 1);
  static_assert(find_if_v<std::is_pointer, NumericTypes> == npos<NumericTypes>);

  // push/pop
  static_assert(std::is_same_v<push_back_t<EvalTypes,Stochastic>, type_list<Residual,Jacobian,Tangent,Stochastic>>);
  static_assert(std::is_same_v<push_front_t<EvalTypes,Stochastic>, type_list<Stochastic,Residual,Jacobian,Tangent>>);
  static_assert(std::is_same_v<push_back_t<type_list<>,int>, type_list<int>>);
  static_assert(std::is_same_v<pop_front_t<EvalTypes>, type_list<Jacobian,Tangent>>);
  static_assert(std::is_same_v<pop_back_t<EvalTypes>, type_list<Residual,Jacobian>>);
  static_assert(std::is_same_v<pop_front_t<type_list<int>>, type_list<>>);
  static_assert(std::is_same_v<pop_back_t<type_list<int>>, type_list<>>);

  // concat, remove, unique
  static_assert(std::is_same_v<concat_t<type_list<int>,type_list<>,type_list<float,double>>, NumericTypes>);
  using WithDuplicates = type_list<int, float, int, double, int>;
  static_assert(std::is_same_v<remove_t<Jacobian,EvalTypes>, type_list<Residual,Tangent>>);
  static_assert(std::is_same_v<remove_t<Stochastic,EvalTypes>, EvalTypes>);
  static_assert(std::is_same_v<remove_t<int,WithDuplicates>, type_list<float,double>>);
  static_assert(std::is_same_v<unique_t<WithDuplicates>, NumericTypes>);
  static_assert(std::is_same_v<unique_t<EvalTypes>, EvalTypes>);
  static_assert(std::is_same_v<unique_t<type_list<>>, type_list<>>);

  // concepts
  static_assert(TypeList<EvalTypes>);
  static_assert(!TypeList<int>);
  static_assert(TypeSequence<EvalTypes>);
  static_assert(TypeSequence<Sacado::mpl::vector<int>>);
  static_assert(!TypeSequence<int>);

  template <typename T> std::string name();
  template <> std::string name<Residual>() { return "R"; }
  template <> std::string name<Jacobian>() { return "J"; }
  template <> std::string name<Tangent>() { return "T"; }

}

TEUCHOS_UNIT_TEST(type_list, for_each_type)
{
  std::vector<std::string> names;
  for_each_type<EvalTypes>([&]<typename EvalT>() { names.push_back(name<EvalT>()); });
  TEST_EQUALITY(names.size(), 3);
  TEST_EQUALITY(names[0], "R");
  TEST_EQUALITY(names[1], "J");
  TEST_EQUALITY(names[2], "T");

  int count = 0;
  for_each_type<type_list<>>([&]<typename T>() { ++count; });
  TEST_EQUALITY(count, 0);
}

TEUCHOS_UNIT_TEST(type_list, for_each_index)
{
  std::vector<std::size_t> indices;
  for_each_index<size_v<EvalTypes>>([&]<std::size_t I>() {
    static_assert(find_v<at_t<EvalTypes,I>, EvalTypes> == I);
    indices.push_back(I);
  });
  TEST_EQUALITY(indices.size(), 3);
  for (std::size_t i=0; i < indices.size(); ++i)
    TEST_EQUALITY(indices[i], i);
}

// Sacado::mpl::vector is still accepted by the query algorithms, but deprecated.
TEUCHOS_UNIT_TEST(type_list, sacado_vector_compatibility)
{
  using SacadoTypes = Sacado::mpl::vector<Residual, Jacobian, Tangent>;
  static_assert(std::is_same_v<to_type_list_t<SacadoTypes>, EvalTypes>);
  static_assert(size_v<SacadoTypes> == 3);
  static_assert(std::is_same_v<at_t<SacadoTypes,2>, Tangent>);
  static_assert(find_v<Jacobian, SacadoTypes> == 1);
  static_assert(!found_v<Stochastic, SacadoTypes>);

  std::vector<std::string> names;
  for_each_type<SacadoTypes>([&]<typename EvalT>() { names.push_back(name<EvalT>()); });
  TEST_EQUALITY(names.size(), 3);
  TEST_EQUALITY(names[2], "T");
}

// Sacado::mpl algorithms also work on a type_list, for downstream code
// that applies them to a Traits::EvalTypes.
namespace {
  struct AppendName {
    std::vector<std::string>& names;
    template <typename T> void operator()(T) const { names.push_back(name<T>()); }
  };
}

TEUCHOS_UNIT_TEST(type_list, sacado_mpl_interop)
{
  static_assert(Sacado::mpl::size<EvalTypes>::value == 3);
  static_assert(Sacado::mpl::find<EvalTypes,Jacobian>::value == 1);
  static_assert(std::is_same_v<Sacado::mpl::at<EvalTypes,2>::type, Tangent>);

  std::vector<std::string> names;
  Sacado::mpl::for_each_no_kokkos<EvalTypes>(AppendName{names});
  TEST_EQUALITY(names.size(), 3);
  TEST_EQUALITY(names[0], "R");
  TEST_EQUALITY(names[1], "J");
  TEST_EQUALITY(names[2], "T");
}
