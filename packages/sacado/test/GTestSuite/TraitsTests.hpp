// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef TRAITSTESTS_HPP
#define TRAITSTESTS_HPP

#include <type_traits>

// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Random.hpp"
#include "Sacado_mpl_apply.hpp"

// gtest includes
#include <gtest/gtest.h>

// A class for testing Sacado::Traits definitions for Sacado AD types
template <class ADType>
class TraitsTests : public ::testing::Test {
protected:

  typedef typename Sacado::mpl::apply<ADType,double>::type ad1_t;
  typedef typename Sacado::mpl::apply<ADType,ad1_t>::type ad2_t;

  // Random number generator
  Sacado::Random<double> urand;

  // GTest creates the test fixture as a child class, so ad1_t, ad2_t are
  // not visible.  Create some data members of these types to get the types
  // in the test cases
  ad1_t ad1;
  ad2_t ad2;

  TraitsTests() : urand(), ad1(), ad2() {}
  ~TraitsTests() {}

}; // class TraitsTests

TYPED_TEST_SUITE_P(TraitsTests);

TYPED_TEST_P(TraitsTests, testScalarType) {
  typedef decltype(this->ad1) ad1_t;
  typedef decltype(this->ad2) ad2_t;

  bool same = std::is_same< typename Sacado::ScalarType<ad1_t>::type, double >::value;
  ASSERT_TRUE(same == true);

same = std::is_same< typename Sacado::ScalarType<ad2_t>::type,double >::value;
  ASSERT_TRUE(same == true);
}

TYPED_TEST_P(TraitsTests, testValueType) {
  typedef decltype(this->ad1) ad1_t;
  typedef decltype(this->ad2) ad2_t;

  bool same = std::is_same< typename Sacado::ValueType<ad1_t>::type,double >::value;
  ASSERT_TRUE(same == true);

  same = std::is_same< typename Sacado::ValueType<ad2_t>::type,ad1_t >::value;
  ASSERT_TRUE(same == true);
}

TYPED_TEST_P(TraitsTests, testIsADType) {
  typedef decltype(this->ad1) ad1_t;
  typedef decltype(this->ad2) ad2_t;

  ASSERT_TRUE(Sacado::IsADType<ad1_t>::value == true);
  ASSERT_TRUE(Sacado::IsADType<ad2_t>::value == true);
}

TYPED_TEST_P(TraitsTests, testIsScalarType) {
  typedef decltype(this->ad1) ad1_t;
  typedef decltype(this->ad2) ad2_t;

  ASSERT_TRUE(Sacado::IsScalarType<ad1_t>::value == false);
  ASSERT_TRUE(Sacado::IsScalarType<ad2_t>::value == false);
}

TYPED_TEST_P(TraitsTests, testValue) {
  typedef decltype(this->ad1) ad1_t;
  typedef decltype(this->ad2) ad2_t;

  double val = this->urand.number();
  ad1_t a(val);
  ASSERT_TRUE(Sacado::Value<ad1_t>::eval(a) == val);

  ad2_t b(a);
  ASSERT_TRUE(Sacado::Value<ad2_t>::eval(b) == a);
}

TYPED_TEST_P(TraitsTests, testScalarValue) {
  typedef decltype(this->ad1) ad1_t;
  typedef decltype(this->ad2) ad2_t;

  double val = this->urand.number();
  ad1_t a(val);
  ASSERT_TRUE(Sacado::ScalarValue<ad1_t>::eval(a) == val);

  ad2_t b(a);
  ASSERT_TRUE(Sacado::ScalarValue<ad2_t>::eval(b) == val);
}

TYPED_TEST_P(TraitsTests, testStringName) {
  typedef decltype(this->ad1) ad1_t;
  typedef decltype(this->ad2) ad2_t;

  // Currently we can't check the string name, here we are just making sure
  // it compiles
  Sacado::StringName<ad1_t>::eval();
  Sacado::StringName<ad2_t>::eval();
  // ASSERT_TRUE(Sacado::StringName<ad1_t>::eval() == name + "< double, double >");
  // ASSERT_TRUE(Sacado::StringName<ad2_t>::eval() == name + "< " + name + "< double, double >, double >");
}

REGISTER_TYPED_TEST_SUITE_P(
  TraitsTests,
  testScalarType,
  testValueType,
  testIsADType,
  testIsScalarType,
  testValue,
  testScalarValue,
  testStringName
  );

#endif // TRAITSTESTS_HPP
