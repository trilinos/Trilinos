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

#include "gtest/gtest.h"
#include "stk_util/util/CSet.hpp"  // for CSet, stk
#include "stk_unit_test_utils/ObjectLifetimeSpy.hpp"


namespace {

using namespace stk::unit_test_util;

class ObjectLifetimeSpyA : public ObjectLifetimeSpy
{
public:
  ObjectLifetimeSpyA(int value)
    : ObjectLifetimeSpy(value)
  {}
};

class ObjectLifetimeSpyB : public ObjectLifetimeSpy
{
public:
  ObjectLifetimeSpyB(int value)
    : ObjectLifetimeSpy(value)
  {}
};


TEST(TestCSetComparison, NotEqual) {
  const std::type_info * intType = &typeid(int);
  const std::type_info * unsignedType = &typeid(unsigned);

  stk::cset::less_cset compare;
  EXPECT_TRUE(compare(intType, unsignedType));
  EXPECT_FALSE(compare(unsignedType, intType));
}

TEST(TestCSetComparison, Equal) {
  const std::type_info * intType1 = &typeid(int);
  const std::type_info * intType2 = &typeid(int);

  stk::cset::less_cset compare;
  EXPECT_FALSE(compare(intType1, intType2));
  EXPECT_FALSE(compare(intType2, intType1));
}


class TestCSet : public ::testing::Test
{
 protected:
  virtual void SetUp() override {
    objectLifetimeSpy_clearCounts();
  }

  virtual void TearDown() override {
    EXPECT_TRUE(objectLifetimeSpy_checkBalancedConstructionsDestructions());
  }
};


TEST_F(TestCSet, NoDelete_InsertSingleObject)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    const ObjectLifetimeSpy* insertedObj = cs.insert_no_delete<ObjectLifetimeSpy>(&obj);
    EXPECT_EQ(insertedObj->value(), 1);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, NoDelete_InsertTwoObjects)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpyA objA(1);
    ObjectLifetimeSpyB objB(2);
    const ObjectLifetimeSpyA* insertedObjA = cs.insert_no_delete<ObjectLifetimeSpyA>(&objA);
    const ObjectLifetimeSpyB* insertedObjB = cs.insert_no_delete<ObjectLifetimeSpyB>(&objB);
    EXPECT_EQ(insertedObjA->value(), 1);
    EXPECT_EQ(insertedObjB->value(), 2);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  2);
}

TEST_F(TestCSet, NoDelete_DoubleInsert_NoReplace)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj1(1);
    ObjectLifetimeSpy obj2(2);
    const ObjectLifetimeSpy* insertedObjFirst  = cs.insert_no_delete<ObjectLifetimeSpy>(&obj1);
    const ObjectLifetimeSpy* insertedObjSecond = cs.insert_no_delete<ObjectLifetimeSpy>(&obj2);
    EXPECT_EQ(insertedObjFirst->value(), 1);
    EXPECT_EQ(insertedObjSecond->value(), 1);
    EXPECT_EQ(insertedObjFirst, insertedObjSecond);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  2);
}

TEST_F(TestCSet, GetNotFound)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    const ObjectLifetimeSpy* getObj = cs.get<ObjectLifetimeSpy>();
    EXPECT_EQ(getObj, nullptr);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, NoDelete_GetSingleObject)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    const ObjectLifetimeSpy* insertedObj = cs.insert_no_delete<ObjectLifetimeSpy>(&obj);
    EXPECT_EQ(insertedObj->value(), 1);

    const ObjectLifetimeSpy* getObj = cs.get<ObjectLifetimeSpy>();
    EXPECT_EQ(getObj->value(), 1);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, NoDelete_GetTwoObjects)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpyA objA(1);
    ObjectLifetimeSpyB objB(2);
    const ObjectLifetimeSpyA* insertedObjA = cs.insert_no_delete<ObjectLifetimeSpyA>(&objA);
    const ObjectLifetimeSpyB* insertedObjB = cs.insert_no_delete<ObjectLifetimeSpyB>(&objB);
    EXPECT_EQ(insertedObjA->value(), 1);
    EXPECT_EQ(insertedObjB->value(), 2);

    const ObjectLifetimeSpyA* getObjA = cs.get<ObjectLifetimeSpyA>();
    const ObjectLifetimeSpyB* getObjB = cs.get<ObjectLifetimeSpyB>();
    EXPECT_EQ(getObjA->value(), 1);
    EXPECT_EQ(getObjB->value(), 2);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  2);
}

TEST_F(TestCSet, RemoveFromEmpty)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    EXPECT_FALSE(cs.remove<ObjectLifetimeSpy>(&obj));
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, NoDelete_Remove)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    const ObjectLifetimeSpy* insertedObj = cs.insert_no_delete<ObjectLifetimeSpy>(&obj);
    EXPECT_EQ(insertedObj->value(), 1);

    EXPECT_TRUE(cs.remove<ObjectLifetimeSpy>(&obj));
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(), 0);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, NoDelete_RemoveSameTwice)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    const ObjectLifetimeSpy* insertedObj = cs.insert_no_delete<ObjectLifetimeSpy>(&obj);
    EXPECT_EQ(insertedObj->value(), 1);

    EXPECT_TRUE (cs.remove<ObjectLifetimeSpy>(&obj));
    EXPECT_FALSE(cs.remove<ObjectLifetimeSpy>(&obj));
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, NoDelete_RemoveTwoObjects)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpyA objA(1);
    ObjectLifetimeSpyB objB(2);
    const ObjectLifetimeSpyA* insertedObjA = cs.insert_no_delete<ObjectLifetimeSpyA>(&objA);
    const ObjectLifetimeSpyB* insertedObjB = cs.insert_no_delete<ObjectLifetimeSpyB>(&objB);
    EXPECT_EQ(insertedObjA->value(), 1);
    EXPECT_EQ(insertedObjB->value(), 2);

    EXPECT_TRUE(cs.remove<ObjectLifetimeSpyA>(&objA));
    EXPECT_TRUE(cs.remove<ObjectLifetimeSpyB>(&objB));
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(), 0);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  2);
}


TEST_F(TestCSet, WithDelete_InsertSingleObject)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy* obj = new ObjectLifetimeSpy(1);
    const ObjectLifetimeSpy* insertedObj = cs.insert_with_delete<ObjectLifetimeSpy>(obj);
    EXPECT_EQ(insertedObj->value(), 1);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, WithDelete_InsertTwoObjects)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpyA* objA = new ObjectLifetimeSpyA(1);
    ObjectLifetimeSpyB* objB = new ObjectLifetimeSpyB(2);
    const ObjectLifetimeSpyA* insertedObjA = cs.insert_with_delete<ObjectLifetimeSpyA>(objA);
    const ObjectLifetimeSpyB* insertedObjB = cs.insert_with_delete<ObjectLifetimeSpyB>(objB);
    EXPECT_EQ(insertedObjA->value(), 1);
    EXPECT_EQ(insertedObjB->value(), 2);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  2);
}

TEST_F(TestCSet, WithDelete_DoubleInsert_NoReplace)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy* obj1 = new ObjectLifetimeSpy(1);
    ObjectLifetimeSpy* obj2 = new ObjectLifetimeSpy(2);
    const ObjectLifetimeSpy* insertedObjFirst  = cs.insert_with_delete<ObjectLifetimeSpy>(obj1);
    const ObjectLifetimeSpy* insertedObjSecond = cs.insert_with_delete<ObjectLifetimeSpy>(obj2);
    EXPECT_EQ(insertedObjFirst->value(), 1);
    EXPECT_EQ(insertedObjSecond->value(), 1);
    EXPECT_EQ(insertedObjFirst, insertedObjSecond);
    delete obj2;
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  2);
}

TEST_F(TestCSet, WithDelete_GetSingleObject)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy* obj = new ObjectLifetimeSpy(1);
    const ObjectLifetimeSpy* insertedObj = cs.insert_with_delete<ObjectLifetimeSpy>(obj);
    EXPECT_EQ(insertedObj->value(), 1);

    const ObjectLifetimeSpy* getObj = cs.get<ObjectLifetimeSpy>();
    EXPECT_EQ(getObj->value(), 1);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, WithDelete_GetTwoObjects)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpyA* objA = new ObjectLifetimeSpyA(1);
    ObjectLifetimeSpyB* objB = new ObjectLifetimeSpyB(2);
    const ObjectLifetimeSpyA* insertedObjA = cs.insert_with_delete<ObjectLifetimeSpyA>(objA);
    const ObjectLifetimeSpyB* insertedObjB = cs.insert_with_delete<ObjectLifetimeSpyB>(objB);
    EXPECT_EQ(insertedObjA->value(), 1);
    EXPECT_EQ(insertedObjB->value(), 2);

    const ObjectLifetimeSpyA* getObjA = cs.get<ObjectLifetimeSpyA>();
    const ObjectLifetimeSpyB* getObjB = cs.get<ObjectLifetimeSpyB>();
    EXPECT_EQ(getObjA->value(), 1);
    EXPECT_EQ(getObjB->value(), 2);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  2);
}

TEST_F(TestCSet, WithDelete_Remove)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy* obj = new ObjectLifetimeSpy(1);
    const ObjectLifetimeSpy* insertedObj = cs.insert_with_delete<ObjectLifetimeSpy>(obj);
    EXPECT_EQ(insertedObj->value(), 1);

    EXPECT_TRUE(cs.remove<ObjectLifetimeSpy>(obj));
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(), 0);
    delete obj;
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, WithDelete_RemoveSameTwice)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy* obj = new ObjectLifetimeSpy(1);
    const ObjectLifetimeSpy* insertedObj = cs.insert_with_delete<ObjectLifetimeSpy>(obj);
    EXPECT_EQ(insertedObj->value(), 1);

    EXPECT_TRUE (cs.remove<ObjectLifetimeSpy>(obj));
    EXPECT_FALSE(cs.remove<ObjectLifetimeSpy>(obj));
    delete obj;
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  1);
}

TEST_F(TestCSet, WithDelete_RemoveTwoObjects)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpyA* objA = new ObjectLifetimeSpyA(1);
    ObjectLifetimeSpyB* objB = new ObjectLifetimeSpyB(2);
    const ObjectLifetimeSpyA* insertedObjA = cs.insert_with_delete<ObjectLifetimeSpyA>(objA);
    const ObjectLifetimeSpyB* insertedObjB = cs.insert_with_delete<ObjectLifetimeSpyB>(objB);
    EXPECT_EQ(insertedObjA->value(), 1);
    EXPECT_EQ(insertedObjB->value(), 2);

    EXPECT_TRUE(cs.remove<ObjectLifetimeSpyA>(objA));
    EXPECT_TRUE(cs.remove<ObjectLifetimeSpyB>(objB));
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(), 0);
    delete objA;
    delete objB;
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),  2);
}


TEST_F(TestCSet, CopyConstructEmpty)
{
  {
    stk::CSet cs;
    stk::CSet csCopy = cs;
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
}

TEST_F(TestCSet, CopyConstructNoDelete)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    {
      const ObjectLifetimeSpy* insertedObj = cs.insert_no_delete<ObjectLifetimeSpy>(&obj);
      EXPECT_EQ(insertedObj->value(), 1);

      stk::CSet csCopy = cs;

      const ObjectLifetimeSpy* copiedGetObj = csCopy.get<ObjectLifetimeSpy>();
      EXPECT_EQ(copiedGetObj->value(), 1);
      EXPECT_EQ(insertedObj, copiedGetObj);
    }

    EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(TestCSet, CopyConstructWithDelete)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy * obj = new ObjectLifetimeSpy(1);
    {
      const ObjectLifetimeSpy* insertedObj = cs.insert_with_delete<ObjectLifetimeSpy>(obj);
      EXPECT_EQ(insertedObj->value(), 1);

      stk::CSet csCopy = cs;

      const ObjectLifetimeSpy* copiedGetObj = csCopy.get<ObjectLifetimeSpy>();
      EXPECT_EQ(copiedGetObj->value(), 1);
      EXPECT_EQ(insertedObj, copiedGetObj);
    }

    EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}


TEST_F(TestCSet, CopyAssignEmpty)
{
  {
    stk::CSet cs;
    {
      stk::CSet secondCs;
      secondCs = cs;
    }
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
}

TEST_F(TestCSet, CopyAssignNoDelete)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    {
      const ObjectLifetimeSpy* insertedObj = cs.insert_no_delete<ObjectLifetimeSpy>(&obj);
      EXPECT_EQ(insertedObj->value(), 1);

      stk::CSet secondCs;
      secondCs = cs;

      const ObjectLifetimeSpy* assignedGetObj = secondCs.get<ObjectLifetimeSpy>();
      EXPECT_EQ(assignedGetObj->value(), 1);
      EXPECT_EQ(insertedObj, assignedGetObj);
    }

    EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(TestCSet, CopyAssignWithDelete)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy * obj = new ObjectLifetimeSpy(1);
    {
      const ObjectLifetimeSpy* insertedObj = cs.insert_with_delete<ObjectLifetimeSpy>(obj);
      EXPECT_EQ(insertedObj->value(), 1);

      stk::CSet secondCs;
      secondCs = cs;

      const ObjectLifetimeSpy* assignedGetObj = secondCs.get<ObjectLifetimeSpy>();
      EXPECT_EQ(assignedGetObj->value(), 1);
      EXPECT_EQ(insertedObj, assignedGetObj);
    }

    EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}


TEST_F(TestCSet, MoveConstructEmpty)
{
  {
    stk::CSet cs;
    stk::CSet csCopy = std::move(cs);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
}

TEST_F(TestCSet, MoveConstructNoDelete)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    {
      const ObjectLifetimeSpy* insertedObj = cs.insert_no_delete<ObjectLifetimeSpy>(&obj);
      EXPECT_EQ(insertedObj->value(), 1);

      stk::CSet csCopy = std::move(cs);

      const ObjectLifetimeSpy* copiedGetObj = csCopy.get<ObjectLifetimeSpy>();
      EXPECT_EQ(copiedGetObj->value(), 1);
      EXPECT_EQ(insertedObj, copiedGetObj);
    }

    EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(TestCSet, MoveConstructWithDelete)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy * obj = new ObjectLifetimeSpy(1);
    {
      const ObjectLifetimeSpy* insertedObj = cs.insert_with_delete<ObjectLifetimeSpy>(obj);
      EXPECT_EQ(insertedObj->value(), 1);

      stk::CSet csCopy = std::move(cs);

      const ObjectLifetimeSpy* copiedGetObj = csCopy.get<ObjectLifetimeSpy>();
      EXPECT_EQ(copiedGetObj->value(), 1);
      EXPECT_EQ(insertedObj, copiedGetObj);
    }

    EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}


TEST_F(TestCSet, MoveAssignEmpty)
{
  {
    stk::CSet cs;
    {
      stk::CSet secondCs;
      secondCs = std::move(cs);
    }
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
}

TEST_F(TestCSet, MoveAssignNoDelete)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy obj(1);
    {
      const ObjectLifetimeSpy* insertedObj = cs.insert_no_delete<ObjectLifetimeSpy>(&obj);
      EXPECT_EQ(insertedObj->value(), 1);

      stk::CSet secondCs;
      secondCs = std::move(cs);

      const ObjectLifetimeSpy* assignedGetObj = secondCs.get<ObjectLifetimeSpy>();
      EXPECT_EQ(assignedGetObj->value(), 1);
      EXPECT_EQ(insertedObj, assignedGetObj);
    }

    EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(TestCSet, MoveAssignWithDelete)
{
  {
    stk::CSet cs;
    ObjectLifetimeSpy * obj = new ObjectLifetimeSpy(1);
    {
      const ObjectLifetimeSpy* insertedObj = cs.insert_with_delete<ObjectLifetimeSpy>(obj);
      EXPECT_EQ(insertedObj->value(), 1);

      stk::CSet secondCs;
      secondCs = std::move(cs);

      const ObjectLifetimeSpy* assignedGetObj = secondCs.get<ObjectLifetimeSpy>();
      EXPECT_EQ(assignedGetObj->value(), 1);
      EXPECT_EQ(insertedObj, assignedGetObj);
    }

    EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
    EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
    EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
  }

  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveAssignments(),   0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

}
