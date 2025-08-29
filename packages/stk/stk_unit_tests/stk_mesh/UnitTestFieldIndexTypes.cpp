// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <gtest/gtest.h>
#include <stk_mesh/base/FieldIndexTypes.hpp>

namespace {

void strong_function_component(stk::mesh::ComponentIdx /*comp*/) {}
void strong_function_copy(stk::mesh::CopyIdx /*copy*/) {}
void strong_function_entity(stk::mesh::EntityIdx /*entity*/) {}

//------------------------------------------------------------------------------
TEST(FieldIndexTypes, constructionAndArgumentPassing_component)
{
  stk::mesh::ComponentIdx componentIdx0(0);     // Same underlying type regardless of construction type
  stk::mesh::ComponentIdx componentIdx0u(0u);
  stk::mesh::ComponentIdx componentIdx0l(0l);
  stk::mesh::ComponentIdx componentIdx0ul(0ul);

  // strong_function_component(0);           // Does not compile
  strong_function_component(0_comp);
  // strong_function_component(0_copy);      // Does not compile
  // strong_function_component(0_entity);    // Does not compile
  strong_function_component(componentIdx0);
  strong_function_component(componentIdx0u);
  strong_function_component(componentIdx0l);
  strong_function_component(componentIdx0ul);

  strong_function_component(x_comp);
  strong_function_component(y_comp);
  strong_function_component(z_comp);

  EXPECT_EQ(x_comp, stk::mesh::ComponentIdx(0));
  EXPECT_EQ(y_comp, stk::mesh::ComponentIdx(1));
  EXPECT_EQ(z_comp, stk::mesh::ComponentIdx(2));
}

//------------------------------------------------------------------------------
TEST(FieldIndexTypes, constructionAndArgumentPassing_copy)
{
  stk::mesh::CopyIdx copyIdx0(0);     // Same underlying type regardless of construction type
  stk::mesh::CopyIdx copyIdx0u(0u);
  stk::mesh::CopyIdx copyIdx0l(0l);
  stk::mesh::CopyIdx copyIdx0ul(0ul);

  // strong_function_copy(0);              // Does not compile
  // strong_function_copy(0_comp);         // Does not compile
  strong_function_copy(0_copy);
  // strong_function_copy(0_entity);       // Does not compile
  strong_function_copy(copyIdx0);
  strong_function_copy(copyIdx0u);
  strong_function_copy(copyIdx0l);
  strong_function_copy(copyIdx0ul);
}

//------------------------------------------------------------------------------
TEST(FieldIndexTypes, constructionAndArgumentPassing_entity)
{
  stk::mesh::EntityIdx entityIdx0(0);     // Same underlying type regardless of construction type
  stk::mesh::EntityIdx entityIdx0u(0u);
  stk::mesh::EntityIdx entityIdx0l(0l);
  stk::mesh::EntityIdx entityIdx0ul(0ul);

  strong_function_entity(0);  // Don't love this, but must be implicitly convertible from int for Kokkos support
  // strong_function_entity(0_comp);         // Does not compile
  // strong_function_entity(0_copy);         // Does not compile
  strong_function_entity(0_entity);
  strong_function_entity(entityIdx0);
  strong_function_entity(entityIdx0u);
  strong_function_entity(entityIdx0l);
  strong_function_entity(entityIdx0ul);
}


using TestTypes = ::testing::Types<stk::mesh::ComponentIdx, stk::mesh::CopyIdx, stk::mesh::EntityIdx>;
using TestTypesInterop = ::testing::Types<std::pair<stk::mesh::ComponentIdx, int>,
                                        std::pair<stk::mesh::ComponentIdx, unsigned int>,
                                        std::pair<stk::mesh::ComponentIdx, long>,
                                        std::pair<stk::mesh::ComponentIdx, unsigned int>,

                                        std::pair<stk::mesh::CopyIdx, int>,
                                        std::pair<stk::mesh::CopyIdx, unsigned int>,
                                        std::pair<stk::mesh::CopyIdx, long>,
                                        std::pair<stk::mesh::CopyIdx, unsigned int>,

                                        std::pair<stk::mesh::EntityIdx, int>,
                                        std::pair<stk::mesh::EntityIdx, unsigned int>,
                                        std::pair<stk::mesh::EntityIdx, long>,
                                        std::pair<stk::mesh::EntityIdx, unsigned int>>;


template <typename T>
class TestFieldIndexTypes : public testing::Test {};
TYPED_TEST_SUITE(TestFieldIndexTypes, TestTypes,);

template <typename T>
class TestFieldIndexTypesInterop : public testing::Test {};
TYPED_TEST_SUITE(TestFieldIndexTypesInterop, TestTypesInterop,);

TYPED_TEST(TestFieldIndexTypes, implicitConversionOperators)
{
  using IndexType = TypeParam;

  int varInt = IndexType(123);
  unsigned varUnsigned = IndexType(123);
  int varLong = IndexType(123);
  unsigned long varUnsignedLong = IndexType(123);

  EXPECT_EQ(varInt, 123);
  EXPECT_EQ(varUnsigned, 123u);
  EXPECT_EQ(varLong, 123l);
  EXPECT_EQ(varUnsignedLong, 123lu);
}

TYPED_TEST(TestFieldIndexTypes, incrementAndDecrementOperators)
{
  using IndexType = TypeParam;

  IndexType index(0);
  EXPECT_EQ(static_cast<int>(index), 0);

  EXPECT_EQ(static_cast<int>(++index), 1);
  EXPECT_EQ(static_cast<int>(index++), 1);
  EXPECT_EQ(static_cast<int>(index), 2);

  EXPECT_EQ(static_cast<int>(--index), 1);
  EXPECT_EQ(static_cast<int>(index--), 1);
  EXPECT_EQ(static_cast<int>(index), 0);

  index += 5;
  EXPECT_EQ(static_cast<int>(index), 5);
  index -= 5;
  EXPECT_EQ(static_cast<int>(index), 0);

  index += 50u;
  EXPECT_EQ(static_cast<int>(index), 50);
  index -= 50u;
  EXPECT_EQ(static_cast<int>(index), 0);

  index += 500l;
  EXPECT_EQ(static_cast<int>(index), 500);
  index -= 500l;
  EXPECT_EQ(static_cast<int>(index), 0);

  index += 5000lu;
  EXPECT_EQ(static_cast<int>(index), 5000);
  index -= 5000lu;
  EXPECT_EQ(static_cast<int>(index), 0);

  IndexType indexDelta(12);
  index += indexDelta;
  EXPECT_EQ(static_cast<int>(index), 12);
  index -= indexDelta;
  EXPECT_EQ(static_cast<int>(index), 0);

  enum { DUMMY_VALUE = 3 };  // Needed for Kokkos support
  index += DUMMY_VALUE;
  EXPECT_EQ(static_cast<int>(index), 3);
  index -= DUMMY_VALUE;
  EXPECT_EQ(static_cast<int>(index), 0);
}

TYPED_TEST(TestFieldIndexTypes, copyAndMoveConstructorsAndAssignment)
{
  using IndexType = TypeParam;

  IndexType base(5);

  IndexType copyConstructed = base;
  EXPECT_EQ(static_cast<int>(copyConstructed), 5);

  IndexType moveConstructed = IndexType(10);
  EXPECT_EQ(static_cast<int>(moveConstructed), 10);

  IndexType copyAssigned(0);
  copyAssigned = base;
  EXPECT_EQ(static_cast<int>(copyAssigned), 5);

  IndexType moveAssigned(0);
  moveAssigned = IndexType(10);
  EXPECT_EQ(static_cast<int>(moveAssigned), 10);

  IndexType assigned(0);
  assigned = 10;
  EXPECT_EQ(static_cast<int>(assigned), 10);

  assigned = 20u;
  EXPECT_EQ(static_cast<int>(assigned), 20);

  assigned = 30l;
  EXPECT_EQ(static_cast<int>(assigned), 30);

  assigned = 40lu;
  EXPECT_EQ(static_cast<int>(assigned), 40);
}

TYPED_TEST(TestFieldIndexTypes, arithmeticOperators)
{
  using IndexType = TypeParam;

  IndexType value(1);
  IndexType rhs(12);

  EXPECT_EQ(static_cast<int>(1), 1);

  value = value + rhs;
  EXPECT_EQ(static_cast<int>(value), 13);

  value = value - rhs;
  EXPECT_EQ(static_cast<int>(value), 1);

  value = value * rhs;
  EXPECT_EQ(static_cast<int>(value), 12);

  value = value / rhs;
  EXPECT_EQ(static_cast<int>(value), 1);


  value = value + 20;
  EXPECT_EQ(static_cast<int>(value), 21);

  value = value - 20;
  EXPECT_EQ(static_cast<int>(value), 1);

  value = value * 20;
  EXPECT_EQ(static_cast<int>(value), 20);

  value = value / 20;
  EXPECT_EQ(static_cast<int>(value), 1);


  value = 20 + value;
  EXPECT_EQ(static_cast<int>(value), 21);

  value = 22 - value;
  EXPECT_EQ(static_cast<int>(value), 1);

  value = 20 * value;
  EXPECT_EQ(static_cast<int>(value), 20);

  value = 20 / value;
  EXPECT_EQ(static_cast<int>(value), 1);
}

TYPED_TEST(TestFieldIndexTypesInterop, arithmeticOperators)
{
  using IndexType   = typename TypeParam::first_type;
  using BuiltinType = typename TypeParam::second_type;

  IndexType value(3);
  BuiltinType value2 = 4;
  BuiltinType value3 = 2;


  EXPECT_EQ(value  + value2, 7);
  EXPECT_EQ(value2 + value,  7);

  EXPECT_EQ(value  - value3, 1);
  EXPECT_EQ(value2 - value,  1);

  EXPECT_EQ(value  * value2, 12);
  EXPECT_EQ(value2 * value,  12);

  EXPECT_EQ(value  / value3, 1);
  EXPECT_EQ(value2 / value,  1);

  IndexType value_copy = value;
  BuiltinType value2_copy = value2;
  EXPECT_EQ(value_copy  += value2, IndexType(7));
  EXPECT_EQ(value2_copy += value,  BuiltinType(7));

  value_copy = value;
  value2_copy = value2;
  EXPECT_EQ(value_copy  -= value3, IndexType(1));
  EXPECT_EQ(value2_copy -= value,  BuiltinType(1));

  value_copy = value;
  value2_copy = value2;
  EXPECT_EQ(value_copy  *= value2, IndexType(12));
  EXPECT_EQ(value2_copy *= value,  BuiltinType(12));

  value_copy = value;
  value2_copy = value2;
  EXPECT_EQ(value_copy  /= value3, IndexType(1));
  EXPECT_EQ(value2_copy /= value,  BuiltinType(1));
}

TYPED_TEST(TestFieldIndexTypes, comparisonOperators)
{
  using IndexType = TypeParam;

  IndexType value1(1);
  IndexType value2(2);

  EXPECT_EQ(value1 == value1, true);
  EXPECT_EQ(value1 != value1, false);

  EXPECT_EQ(value1 == value2, false);
  EXPECT_EQ(value1 != value2, true);

  EXPECT_EQ(value1 < value1, false);
  EXPECT_EQ(value1 < value2, true);

  EXPECT_EQ(value2 <= value1, false);
  EXPECT_EQ(value1 <= value1, true);
  EXPECT_EQ(value1 <= value2, true);

  EXPECT_EQ(value1 > value1, false);
  EXPECT_EQ(value2 > value1, true);

  EXPECT_EQ(value1 >= value2, false);
  EXPECT_EQ(value1 >= value1, true);
  EXPECT_EQ(value2 >= value1, true);


  EXPECT_EQ(value1 == 1, true);
  EXPECT_EQ(value1 != 1, false);

  EXPECT_EQ(value1 == 2, false);
  EXPECT_EQ(value1 != 2, true);

  EXPECT_EQ(value1 < 1, false);
  EXPECT_EQ(value1 < 2, true);

  EXPECT_EQ(value2 <= 1, false);
  EXPECT_EQ(value1 <= 1, true);
  EXPECT_EQ(value1 <= 2, true);

  EXPECT_EQ(value1 > 1, false);
  EXPECT_EQ(value2 > 1, true);

  EXPECT_EQ(value1 >= 2, false);
  EXPECT_EQ(value1 >= 1, true);
  EXPECT_EQ(value2 >= 1, true);


  EXPECT_EQ(1 == value1, true);
  EXPECT_EQ(1 != value1, false);

  EXPECT_EQ(1 == value2, false);
  EXPECT_EQ(1 != value2, true);

  EXPECT_EQ(1 < value1, false);
  EXPECT_EQ(1 < value2, true);

  EXPECT_EQ(2 <= value1, false);
  EXPECT_EQ(1 <= value1, true);
  EXPECT_EQ(1 <= value2, true);

  EXPECT_EQ(1 > value1, false);
  EXPECT_EQ(2 > value1, true);

  EXPECT_EQ(1 >= value2, false);
  EXPECT_EQ(1 >= value1, true);
  EXPECT_EQ(2 >= value1, true);
}

TYPED_TEST(TestFieldIndexTypes, traditionalForLooping)
{
  using IndexType = TypeParam;

  IndexType lowerBound(0);
  IndexType upperBound(10);

  int counter = 0;
  for (IndexType index = lowerBound; index < upperBound; ++index) {
    counter += static_cast<int>(index);
  }
  EXPECT_EQ(counter, 45);

  counter = 0;
  for (IndexType index = lowerBound; index < 10; ++index) {
    counter += static_cast<int>(index);
  }
  EXPECT_EQ(counter, 45);


  counter = 0;
  for (IndexType index = upperBound-1; index > lowerBound-1; --index) {
    counter += static_cast<int>(index);
  }
  EXPECT_EQ(counter, 45);

  counter = 0;
  for (IndexType index = upperBound-1; index > -1; --index) {
    counter += static_cast<int>(index);
  }
  EXPECT_EQ(counter, 45);
}


struct ComponentThingy {
  ComponentThingy(int size) : m_size(size) {}
  inline stk::mesh::ComponentIdxProxy components() const { return stk::mesh::ComponentIdxProxy(m_size); }
  int m_size;
};

struct CopyThingy {
  CopyThingy(int size) : m_size(size) {}
  inline stk::mesh::CopyIdxProxy copies() const { return stk::mesh::CopyIdxProxy(m_size); }
  int m_size;
};

struct EntityThingy {
  EntityThingy(int size) : m_size(size) {}
  inline stk::mesh::EntityIdxProxy entities() const { return stk::mesh::EntityIdxProxy(m_size); }
  int m_size;
};


TEST(FieldIndexTypes, rangeBasedForLooping_component) {
  ComponentThingy thingy(10);

  int counter = 0;
  for (stk::mesh::ComponentIdx index : thingy.components()) {
    counter += static_cast<int>(index);
  }
  EXPECT_EQ(counter, 45);
}

TEST(FieldIndexTypes, rangeBasedForLooping_copy) {
  CopyThingy thingy(10);

  int counter = 0;
  for (stk::mesh::CopyIdx index : thingy.copies()) {
    counter += static_cast<int>(index);
  }
  EXPECT_EQ(counter, 45);
}

TEST(FieldIndexTypes, rangeBasedForLooping_entity) {
  EntityThingy thingy(10);

  int counter = 0;
  for (stk::mesh::EntityIdx index : thingy.entities()) {
    counter += static_cast<int>(index);
  }
  EXPECT_EQ(counter, 45);
}

TEST(FieldIndexTypes, iteratorLooping_component) {
  ComponentThingy thingy(10);

  int counter = 0;
  for (stk::mesh::ComponentIdxIterator iter = thingy.components().begin(); iter != thingy.components().end(); ++iter) {
    counter += static_cast<int>(*iter);
  }
  EXPECT_EQ(counter, 45);

  counter = 0;
  for (stk::mesh::ComponentIdxIterator iter = thingy.components().rbegin(); iter != thingy.components().rend(); --iter) {
    counter += static_cast<int>(*iter);
  }
  EXPECT_EQ(counter, 45);
}

TEST(FieldIndexTypes, iteratorLooping_copy) {
  CopyThingy thingy(10);

  int counter = 0;
  for (stk::mesh::CopyIdxIterator iter = thingy.copies().begin(); iter != thingy.copies().end(); ++iter) {
    counter += static_cast<int>(*iter);
  }
  EXPECT_EQ(counter, 45);

  counter = 0;
  for (stk::mesh::CopyIdxIterator iter = thingy.copies().rbegin(); iter != thingy.copies().rend(); --iter) {
    counter += static_cast<int>(*iter);
  }
  EXPECT_EQ(counter, 45);
}

TEST(FieldIndexTypes, iteratorLooping_entity) {
  EntityThingy thingy(10);

  int counter = 0;
  for (stk::mesh::EntityIdxIterator iter = thingy.entities().begin(); iter != thingy.entities().end(); ++iter) {
    counter += static_cast<int>(*iter);
  }
  EXPECT_EQ(counter, 45);

  counter = 0;
  for (stk::mesh::EntityIdxIterator iter = thingy.entities().rbegin(); iter != thingy.entities().rend(); --iter) {
    counter += static_cast<int>(*iter);
  }
  EXPECT_EQ(counter, 45);
}

}
