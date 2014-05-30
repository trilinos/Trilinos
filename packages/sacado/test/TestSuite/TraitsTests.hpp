// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef TRAITSTESTS_HPP
#define TRAITSTESTS_HPP

// Sacado includes
#include "Sacado.hpp"
#include "Sacado_Random.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_mpl_is_same.hpp"

// Cppunit includes
#include <cppunit/extensions/HelperMacros.h>

// A class for testing Sacado::Traits definitions for Sacado AD types
template <class ADType>
class TraitsTests : public CppUnit::TestFixture {

  CPPUNIT_TEST_SUITE( TraitsTests );
  
  CPPUNIT_TEST(testScalarType);
  CPPUNIT_TEST(testValueType);
  CPPUNIT_TEST(testIsADType);
  CPPUNIT_TEST(testIsScalarType);
  CPPUNIT_TEST(testValue);
  CPPUNIT_TEST(testScalarValue);
  CPPUNIT_TEST(testStringName);

  CPPUNIT_TEST_SUITE_END();

public:

  TraitsTests();
  ~TraitsTests() {}

  void setUp() {}

  void tearDown() {}

  void testScalarType();
  void testValueType();
  void testIsADType();
  void testIsScalarType();
  void testValue();
  void testScalarValue();
  void testStringName();

protected:

  typedef typename Sacado::mpl::apply<ADType,double>::type ad1_t;
  typedef typename Sacado::mpl::apply<ADType,ad1_t>::type ad2_t;

  // Random number generator
  Sacado::Random<double> urand;

  // Memory pools for DMFad
  Sacado::Fad::MemPoolManager<double> poolManager;
  Sacado::Fad::MemPoolManager< Sacado::Fad::DMFad<double> > poolManager2;

}; // class TraitsTests

template <class ADType>
TraitsTests<ADType>::
TraitsTests() :
  urand(),
  poolManager(1),
  poolManager2(1)
{
  Sacado::Fad::MemPool *pool = poolManager.getMemoryPool(1);
  Sacado::Fad::DMFad<double>::setDefaultPool(pool);

  Sacado::Fad::MemPool *pool2 = poolManager2.getMemoryPool(1);
  Sacado::Fad::DMFad< Sacado::Fad::DMFad<double> >::setDefaultPool(pool2);
}

template <class ADType>
void
TraitsTests<ADType>::
testScalarType() {
  bool same = Sacado::mpl::is_same< typename Sacado::ScalarType<ad1_t>::type, double >::value;
  CPPUNIT_ASSERT(same == true);
  
  same = Sacado::mpl::is_same< typename Sacado::ScalarType<ad2_t>::type,double >::value;
  CPPUNIT_ASSERT(same == true);
}

template <class ADType>
void
TraitsTests<ADType>::
testValueType() {
  bool same = Sacado::mpl::is_same< typename Sacado::ValueType<ad1_t>::type,double >::value;
  CPPUNIT_ASSERT(same == true);

  same = Sacado::mpl::is_same< typename Sacado::ValueType<ad2_t>::type,ad1_t >::value;
  CPPUNIT_ASSERT(same == true);
}

template <class ADType>
void
TraitsTests<ADType>::
testIsADType() {
  CPPUNIT_ASSERT(Sacado::IsADType<ad1_t>::value == true);
  CPPUNIT_ASSERT(Sacado::IsADType<ad2_t>::value == true);
}

template <class ADType>
void
TraitsTests<ADType>::
testIsScalarType() {
  CPPUNIT_ASSERT(Sacado::IsScalarType<ad1_t>::value == false);
  CPPUNIT_ASSERT(Sacado::IsScalarType<ad2_t>::value == false);
}

template <class ADType>
void
TraitsTests<ADType>::
testValue() {
  double val = urand.number();
  ad1_t a(val);
  CPPUNIT_ASSERT(Sacado::Value<ad1_t>::eval(a) == val);

  ad2_t b(a);
  CPPUNIT_ASSERT(Sacado::Value<ad2_t>::eval(b) == a);
}

template <class ADType>
void
TraitsTests<ADType>::
testScalarValue() {
  double val = urand.number();
  ad1_t a(val);
  CPPUNIT_ASSERT(Sacado::ScalarValue<ad1_t>::eval(a) == val);

  ad2_t b(a);
  CPPUNIT_ASSERT(Sacado::ScalarValue<ad2_t>::eval(b) == val);
}

template <class ADType>
void
TraitsTests<ADType>::
testStringName() {
  // Currently we can't check the string name, here we are just making sure
  // it compiles
  Sacado::StringName<ad1_t>::eval();
  Sacado::StringName<ad2_t>::eval();
  // CPPUNIT_ASSERT(Sacado::StringName<ad1_t>::eval() == name + "< double, double >");
  // CPPUNIT_ASSERT(Sacado::StringName<ad2_t>::eval() == name + "< " + name + "< double, double >, double >");
}

#endif // TRAITSTESTS_HPP
