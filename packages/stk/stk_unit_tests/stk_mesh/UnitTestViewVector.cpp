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
#include "stk_mesh/baseImpl/ViewVector.hpp"
#include "stk_unit_test_utils/ObjectLifetimeSpy.hpp"
#include "stk_ngp_test/ngp_test.hpp"

namespace {

//------------------------------------------------------------------------------
TEST(HostViewVector, defaultConstruction)
{
  stk::mesh::impl::ViewVector<int> viewVector;

  EXPECT_EQ(viewVector.name(), std::string(""));
  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);
  EXPECT_EQ(viewVector.empty(), true);
  EXPECT_EQ(viewVector.data(), nullptr);
}

TEST(HostViewVector, nameConstruction)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector");

  EXPECT_EQ(viewVector.name(), std::string("TestViewVector"));
  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);
  EXPECT_EQ(viewVector.empty(), true);
  EXPECT_NE(viewVector.data(), nullptr);  // Actual allocation of zero size, so not null
}

TEST(HostViewVector, nameAndZeroSizeConstruction)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 0);

  EXPECT_EQ(viewVector.name(), std::string("TestViewVector"));
  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);
  EXPECT_EQ(viewVector.empty(), true);
  EXPECT_NE(viewVector.data(), nullptr);  // Actual allocation of zero size, so not null
}

TEST(HostViewVector, nameAndSizeConstruction)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 5);

  EXPECT_EQ(viewVector.name(), std::string("TestViewVector"));
  EXPECT_EQ(viewVector.size(), 5u);
  EXPECT_EQ(viewVector.capacity(), 5u);
  EXPECT_EQ(viewVector.empty(), false);
  EXPECT_NE(viewVector.data(), nullptr);

  for (unsigned i = 0; i < viewVector.size(); ++i) {
    EXPECT_EQ(viewVector[i], int{});
  }
}

TEST(HostViewVector, indexing)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  viewVector[0] = 1;
  viewVector[1] = 2;
  viewVector[2] = 3;

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
}

TEST(HostViewVector, data)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  viewVector[0] = 1;
  viewVector[1] = 2;
  viewVector[2] = 3;

  EXPECT_EQ(viewVector.data()[0], 1);
  EXPECT_EQ(viewVector.data()[1], 2);
  EXPECT_EQ(viewVector.data()[2], 3);
}

TEST(HostViewVector, reserveUp)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  viewVector[0] = 1;
  viewVector[1] = 2;
  viewVector[2] = 3;

  viewVector.reserve(4);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(4);

  EXPECT_EQ(viewVector.size(), 4u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
  EXPECT_EQ(viewVector[3], 4);
}

TEST(HostViewVector, reserveEqual)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  viewVector[0] = 1;
  viewVector[1] = 2;
  viewVector[2] = 3;

  viewVector.reserve(3);  // Same, so no change
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
}

TEST(HostViewVector, reserveDown)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  viewVector[0] = 1;
  viewVector[1] = 2;
  viewVector[2] = 3;

  viewVector.reserve(2);  // Less, so no change
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
}

TEST(HostViewVector, resizeUp)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  viewVector[0] = 1;
  viewVector[1] = 2;
  viewVector[2] = 3;

  viewVector.resize(4);
  EXPECT_EQ(viewVector.size(), 4u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
  EXPECT_EQ(viewVector[3], int{});
}

TEST(HostViewVector, resizeSame)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  viewVector[0] = 1;
  viewVector[1] = 2;
  viewVector[2] = 3;

  viewVector.resize(3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
}

TEST(HostViewVector, resizeDown)
{
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector[0] = 1;
    viewVector[1] = 2;
    viewVector[2] = 3;

    viewVector.resize(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    EXPECT_EQ(viewVector[0], 1);
    EXPECT_EQ(viewVector[1], 2);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize(0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 3u);
  }
}

TEST(HostViewVector, resizeScaleUp_fromEmptyInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<int>::growth_scale, 2u);
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(4);
    EXPECT_EQ(viewVector.size(), 4u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 8u);
  }
}

TEST(HostViewVector, resizeScaleUp_fromPositiveInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<int>::growth_scale, 2u);
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 8u);
  }

  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(4);
    EXPECT_EQ(viewVector.size(), 4u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 8u);
  }

  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize_scale(4);
    EXPECT_EQ(viewVector.size(), 4u);
    EXPECT_EQ(viewVector.capacity(), 6u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize_scale(6);
    EXPECT_EQ(viewVector.size(), 6u);
    EXPECT_EQ(viewVector.capacity(), 6u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize_scale(7);
    EXPECT_EQ(viewVector.size(), 7u);
    EXPECT_EQ(viewVector.capacity(), 12u);
  }
}

TEST(HostViewVector, resizeScaleSame)
{
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 5u);

    viewVector.resize_scale(5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 5u);
  }
}

TEST(HostViewVector, resizeScaleDown)
{
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 3u);
  }
  {
    stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 4);
    EXPECT_EQ(viewVector.size(), 4u);
    EXPECT_EQ(viewVector.capacity(), 4u);

    viewVector.resize_scale(3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
}

TEST(HostViewVector, clearAlreadyEmpty)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 0);
  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);

  viewVector.clear();

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);
}

TEST(HostViewVector, clear)
{
  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector", 1);
  EXPECT_EQ(viewVector.size(), 1u);
  EXPECT_EQ(viewVector.capacity(), 1u);

  viewVector.clear();

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 1u);
}

TEST(HostViewVector, swap)
{
  stk::mesh::impl::ViewVector<int> viewVector1("ViewVector1", 1);
  stk::mesh::impl::ViewVector<int> viewVector2("ViewVector2", 2);

  EXPECT_EQ(viewVector1.name(), std::string("ViewVector1"));
  EXPECT_EQ(viewVector1.size(), 1u);
  EXPECT_EQ(viewVector1.capacity(), 1u);

  EXPECT_EQ(viewVector2.name(), std::string("ViewVector2"));
  EXPECT_EQ(viewVector2.size(), 2u);
  EXPECT_EQ(viewVector2.capacity(), 2u);

  viewVector1[0] = 10;

  viewVector2[0] = 20;
  viewVector2[1] = 21;

  viewVector1.swap(viewVector2);

  EXPECT_EQ(viewVector1.name(), std::string("ViewVector2"));
  EXPECT_EQ(viewVector1.size(), 2u);
  EXPECT_EQ(viewVector1.capacity(), 2u);

  EXPECT_EQ(viewVector2.name(), std::string("ViewVector1"));
  EXPECT_EQ(viewVector2.size(), 1u);
  EXPECT_EQ(viewVector2.capacity(), 1u);

  viewVector1[0] = 20;
  viewVector1[1] = 21;

  viewVector2[0] = 10;
}

TEST(HostViewVector, stdSwap)
{
  stk::mesh::impl::ViewVector<int> viewVector1("ViewVector1", 1);
  stk::mesh::impl::ViewVector<int> viewVector2("ViewVector2", 2);

  EXPECT_EQ(viewVector1.name(), std::string("ViewVector1"));
  EXPECT_EQ(viewVector1.size(), 1u);
  EXPECT_EQ(viewVector1.capacity(), 1u);

  EXPECT_EQ(viewVector2.name(), std::string("ViewVector2"));
  EXPECT_EQ(viewVector2.size(), 2u);
  EXPECT_EQ(viewVector2.capacity(), 2u);

  viewVector1[0] = 10;

  viewVector2[0] = 20;
  viewVector2[1] = 21;

  std::swap(viewVector1, viewVector2);

  EXPECT_EQ(viewVector1.name(), std::string("ViewVector2"));
  EXPECT_EQ(viewVector1.size(), 2u);
  EXPECT_EQ(viewVector1.capacity(), 2u);

  EXPECT_EQ(viewVector2.name(), std::string("ViewVector1"));
  EXPECT_EQ(viewVector2.size(), 1u);
  EXPECT_EQ(viewVector2.capacity(), 1u);

  viewVector1[0] = 20;
  viewVector1[1] = 21;

  viewVector2[0] = 10;
}

TEST(HostViewVector, pushBack_lValue_fromDefaultInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<int>::growth_scale, 2u);

  stk::mesh::impl::ViewVector<int> viewVector;

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);  // Starts at capacity 0

  int value = 1;
  viewVector.push_back(value);
  EXPECT_EQ(viewVector.size(), 1u);
  EXPECT_EQ(viewVector.capacity(), 1u);

  viewVector.push_back(++value);
  EXPECT_EQ(viewVector.size(), 2u);
  EXPECT_EQ(viewVector.capacity(), 2u);

  viewVector.push_back(++value);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(++value);
  EXPECT_EQ(viewVector.size(), 4u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(++value);
  EXPECT_EQ(viewVector.size(), 5u);
  EXPECT_EQ(viewVector.capacity(), 8u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
  EXPECT_EQ(viewVector[3], 4);
  EXPECT_EQ(viewVector[4], 5);
}

TEST(HostViewVector, pushBack_lValue_fromNamedInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<int>::growth_scale, 2u);

  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector");

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);

  int value = 1;
  viewVector.push_back(value);
  EXPECT_EQ(viewVector.size(), 1u);
  EXPECT_EQ(viewVector.capacity(), 1u);

  viewVector.push_back(++value);
  EXPECT_EQ(viewVector.size(), 2u);
  EXPECT_EQ(viewVector.capacity(), 2u);

  viewVector.push_back(++value);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(++value);
  EXPECT_EQ(viewVector.size(), 4u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(++value);
  EXPECT_EQ(viewVector.size(), 5u);
  EXPECT_EQ(viewVector.capacity(), 8u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
  EXPECT_EQ(viewVector[3], 4);
  EXPECT_EQ(viewVector[4], 5);
}

TEST(HostViewVector, pushBack_rValue_fromDefaultInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<int>::growth_scale, 2u);

  stk::mesh::impl::ViewVector<int> viewVector;

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);

  viewVector.push_back(1);
  EXPECT_EQ(viewVector.size(), 1u);
  EXPECT_EQ(viewVector.capacity(), 1u);

  viewVector.push_back(2);
  EXPECT_EQ(viewVector.size(), 2u);
  EXPECT_EQ(viewVector.capacity(), 2u);

  viewVector.push_back(3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(4);
  EXPECT_EQ(viewVector.size(), 4u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(5);
  EXPECT_EQ(viewVector.size(), 5u);
  EXPECT_EQ(viewVector.capacity(), 8u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
  EXPECT_EQ(viewVector[3], 4);
  EXPECT_EQ(viewVector[4], 5);
}

TEST(HostViewVector, pushBack_rValue_fromNamedInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<int>::growth_scale, 2u);

  stk::mesh::impl::ViewVector<int> viewVector("TestViewVector");

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);

  viewVector.push_back(1);
  EXPECT_EQ(viewVector.size(), 1u);
  EXPECT_EQ(viewVector.capacity(), 1u);

  viewVector.push_back(2);
  EXPECT_EQ(viewVector.size(), 2u);
  EXPECT_EQ(viewVector.capacity(), 2u);

  viewVector.push_back(3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(4);
  EXPECT_EQ(viewVector.size(), 4u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.push_back(5);
  EXPECT_EQ(viewVector.size(), 5u);
  EXPECT_EQ(viewVector.capacity(), 8u);

  EXPECT_EQ(viewVector[0], 1);
  EXPECT_EQ(viewVector[1], 2);
  EXPECT_EQ(viewVector[2], 3);
  EXPECT_EQ(viewVector[3], 4);
  EXPECT_EQ(viewVector[4], 5);
}

struct DummyClass {
  DummyClass(int a, int b)
    : m_a(a),
      m_b(b)
  {}

  DummyClass()
    : m_a(0),
      m_b(0)
  {}

  DummyClass(const DummyClass& other) = default;
  DummyClass(DummyClass&& other) = default;
  DummyClass& operator=(const DummyClass& rhs) = default;
  DummyClass& operator=(DummyClass&& rhs) = default;

  friend bool operator==(const DummyClass& lhs, const DummyClass& rhs) {
    return (lhs.m_a == rhs.m_a) && (lhs.m_b == rhs.m_b);
  }

  int m_a;
  int m_b;
};

TEST(HostViewVector, emplaceBack_fromDefaultInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<DummyClass>::growth_scale, 2u);

  stk::mesh::impl::ViewVector<DummyClass> viewVector;

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);  // Starts at capacity 0

  viewVector.emplace_back(1, 1);
  EXPECT_EQ(viewVector.size(), 1u);
  EXPECT_EQ(viewVector.capacity(), 1u);

  viewVector.emplace_back(2, 2);
  EXPECT_EQ(viewVector.size(), 2u);
  EXPECT_EQ(viewVector.capacity(), 2u);

  viewVector.emplace_back(3, 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.emplace_back(4, 4);
  EXPECT_EQ(viewVector.size(), 4u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  viewVector.emplace_back(5, 5);
  EXPECT_EQ(viewVector.size(), 5u);
  EXPECT_EQ(viewVector.capacity(), 8u);

  EXPECT_EQ(viewVector[0], DummyClass(1, 1));
  EXPECT_EQ(viewVector[1], DummyClass(2, 2));
  EXPECT_EQ(viewVector[2], DummyClass(3, 3));
  EXPECT_EQ(viewVector[3], DummyClass(4, 4));
  EXPECT_EQ(viewVector[4], DummyClass(5, 5));
}


//------------------------------------------------------------------------------
TEST(DeviceViewVector, defaultConstruction)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector;

  EXPECT_EQ(viewVector.name(), std::string(""));
  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);
  EXPECT_EQ(viewVector.empty(), true);
  EXPECT_EQ(viewVector.data(), nullptr);
}

TEST(DeviceViewVector, nameConstruction)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector");

  EXPECT_EQ(viewVector.name(), std::string("TestViewVector"));
  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);
  EXPECT_EQ(viewVector.empty(), true);
  EXPECT_NE(viewVector.data(), nullptr);
}

template <typename VECTOR_T>
void check_default_initialization(const VECTOR_T& viewVector) {
  Kokkos::parallel_for(1,
    KOKKOS_LAMBDA(const int&) {
      for (unsigned i = 0; i < viewVector.size(); ++i) {
        NGP_EXPECT_EQ(viewVector[i], int{});
      }
    });
}

NGP_TEST(DeviceViewVector, nameAndZeroSizeConstruction)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 0);

  EXPECT_EQ(viewVector.name(), std::string("TestViewVector"));
  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);
  EXPECT_EQ(viewVector.empty(), true);
  EXPECT_NE(viewVector.data(), nullptr);
}

NGP_TEST(DeviceViewVector, nameAndSizeConstruction)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 5);

  EXPECT_EQ(viewVector.name(), std::string("TestViewVector"));
  EXPECT_EQ(viewVector.size(), 5u);
  EXPECT_EQ(viewVector.capacity(), 5u);
  EXPECT_EQ(viewVector.empty(), false);
  EXPECT_NE(viewVector.data(), nullptr);

  check_default_initialization(viewVector);
}

template <typename VECTOR_T>
void fill_device_values(const VECTOR_T& viewVector)
{
  Kokkos::parallel_for(1,
    KOKKOS_LAMBDA(const int&) {
      for (unsigned i = 0; i < viewVector.size(); ++i) {
        viewVector[i] = i+1;
      }
    });
}

template <typename VECTOR_T>
void check_device_values(const VECTOR_T& viewVector, unsigned numSetValues = 0)
{
  Kokkos::parallel_for(1,
    KOKKOS_LAMBDA(const int&) {
      for (unsigned i = 0; i < viewVector.size(); ++i) {
        if ((numSetValues == 0u) || ((numSetValues > 0u) && (i < numSetValues))) {
          NGP_EXPECT_EQ(viewVector[i], static_cast<int>(i+1));
        }
        else {
          NGP_EXPECT_EQ(viewVector[i], int{});
        }
      }
    });
}

TEST(DeviceViewVector, indexing)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);

  fill_device_values(viewVector);
  check_device_values(viewVector);
}

template <typename VECTOR_T>
void check_device_values_from_data(const VECTOR_T& viewVector, unsigned numSetValues = 0)
{
  Kokkos::parallel_for(1,
    KOKKOS_LAMBDA(const int&) {
      for (unsigned i = 0; i < viewVector.size(); ++i) {
        if ((numSetValues == 0u) || ((numSetValues > 0u) && (i < numSetValues))) {
          NGP_EXPECT_EQ(viewVector.data()[i], static_cast<int>(i+1));
        }
        else {
          NGP_EXPECT_EQ(viewVector.data()[i], int{});
        }
      }
    });
}

TEST(DeviceViewVector, data)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);

  fill_device_values(viewVector);
  check_device_values_from_data(viewVector);
}

TEST(DeviceViewVector, reserveUp)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);

  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  fill_device_values(viewVector);

  viewVector.reserve(4);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  check_device_values(viewVector);
}

TEST(DeviceViewVector, reserveEqual)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);

  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  fill_device_values(viewVector);

  viewVector.reserve(3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  check_device_values(viewVector);
}

TEST(DeviceViewVector, reserveDown)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);

  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  fill_device_values(viewVector);

  viewVector.reserve(2);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);  // No change

  check_device_values(viewVector);
}

TEST(DeviceViewVector, resizeUp)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  fill_device_values(viewVector);

  viewVector.resize(4);
  EXPECT_EQ(viewVector.size(), 4u);
  EXPECT_EQ(viewVector.capacity(), 4u);

  check_device_values(viewVector, 3);
}

TEST(DeviceViewVector, resizeSame)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  fill_device_values(viewVector);

  viewVector.resize(3);
  EXPECT_EQ(viewVector.size(), 3u);
  EXPECT_EQ(viewVector.capacity(), 3u);

  check_device_values(viewVector);
}

TEST(DeviceViewVector, resizeDown)
{
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    fill_device_values(viewVector);

    viewVector.resize(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    check_device_values(viewVector);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize(0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 3u);
  }
}

TEST(DeviceViewVector, resizeScaleUp_fromEmptyInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<int>::growth_scale, 2u);
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(4);
    EXPECT_EQ(viewVector.size(), 4u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 8u);
  }
}

TEST(DeviceViewVector, resizeScaleUp_fromPositiveInit)
{
  ASSERT_EQ(stk::mesh::impl::ViewVector<int>::growth_scale, 2u);
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 8u);
  }

  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(4);
    EXPECT_EQ(viewVector.size(), 4u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 8u);
  }

  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize_scale(4);
    EXPECT_EQ(viewVector.size(), 4u);
    EXPECT_EQ(viewVector.capacity(), 6u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize_scale(6);
    EXPECT_EQ(viewVector.size(), 6u);
    EXPECT_EQ(viewVector.capacity(), 6u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize_scale(7);
    EXPECT_EQ(viewVector.size(), 7u);
    EXPECT_EQ(viewVector.capacity(), 12u);
  }
}

TEST(DeviceViewVector, resizeScaleSame)
{
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);

    viewVector.resize_scale(0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 5u);

    viewVector.resize_scale(5);
    EXPECT_EQ(viewVector.size(), 5u);
    EXPECT_EQ(viewVector.capacity(), 5u);
  }
}

TEST(DeviceViewVector, resizeScaleDown)
{
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    viewVector.resize_scale(0);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);

    viewVector.resize_scale(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 3u);

    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 3u);
  }
  {
    stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 4);
    EXPECT_EQ(viewVector.size(), 4u);
    EXPECT_EQ(viewVector.capacity(), 4u);

    viewVector.resize_scale(3);
    EXPECT_EQ(viewVector.size(), 3u);
    EXPECT_EQ(viewVector.capacity(), 4u);
  }
}

TEST(DeviceViewVector, clearAlreadyEmpty)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 0);
  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);

  viewVector.clear();

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 0u);
}

TEST(DeviceViewVector, clear)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector("TestViewVector", 1);
  EXPECT_EQ(viewVector.size(), 1u);
  EXPECT_EQ(viewVector.capacity(), 1u);

  viewVector.clear();

  EXPECT_EQ(viewVector.size(), 0u);
  EXPECT_EQ(viewVector.capacity(), 1u);
}

TEST(DeviceViewVector, swap)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector1("ViewVector1", 1);
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector2("ViewVector2", 2);

  EXPECT_EQ(viewVector1.name(), std::string("ViewVector1"));
  EXPECT_EQ(viewVector1.size(), 1u);
  EXPECT_EQ(viewVector1.capacity(), 1u);

  EXPECT_EQ(viewVector2.name(), std::string("ViewVector2"));
  EXPECT_EQ(viewVector2.size(), 2u);
  EXPECT_EQ(viewVector2.capacity(), 2u);

  viewVector1.swap(viewVector2);

  EXPECT_EQ(viewVector1.name(), std::string("ViewVector2"));
  EXPECT_EQ(viewVector1.size(), 2u);
  EXPECT_EQ(viewVector1.capacity(), 2u);

  EXPECT_EQ(viewVector2.name(), std::string("ViewVector1"));
  EXPECT_EQ(viewVector2.size(), 1u);
  EXPECT_EQ(viewVector2.capacity(), 1u);
}

TEST(DeviceViewVector, stdSwap)
{
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector1("ViewVector1", 1);
  stk::mesh::impl::ViewVector<int, stk::ngp::MemSpace> viewVector2("ViewVector2", 2);

  EXPECT_EQ(viewVector1.name(), std::string("ViewVector1"));
  EXPECT_EQ(viewVector1.size(), 1u);
  EXPECT_EQ(viewVector1.capacity(), 1u);

  EXPECT_EQ(viewVector2.name(), std::string("ViewVector2"));
  EXPECT_EQ(viewVector2.size(), 2u);
  EXPECT_EQ(viewVector2.capacity(), 2u);

  std::swap(viewVector1, viewVector2);

  EXPECT_EQ(viewVector1.name(), std::string("ViewVector2"));
  EXPECT_EQ(viewVector1.size(), 2u);
  EXPECT_EQ(viewVector1.capacity(), 2u);

  EXPECT_EQ(viewVector2.name(), std::string("ViewVector1"));
  EXPECT_EQ(viewVector2.size(), 1u);
  EXPECT_EQ(viewVector2.capacity(), 1u);
}


using namespace stk::unit_test_util;

class ViewVectorObjectLifetimes : public ::testing::Test {
protected:
  ViewVectorObjectLifetimes() {}
  void SetUp() override {
    objectLifetimeSpy_clearCounts();
  }
  void TearDown() override {
    EXPECT_TRUE(objectLifetimeSpy_checkBalancedConstructionsDestructions());
  }
};

TEST_F(ViewVectorObjectLifetimes, DefaultConstruction)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
}

TEST_F(ViewVectorObjectLifetimes, NamedConstruction)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector("LifetimeVector");
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 0u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
}

TEST_F(ViewVectorObjectLifetimes, ConstructionSize1)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector("LifetimeVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, CopyConstruction)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector("LifetimeVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    auto copiedViewVector(viewVector);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
    EXPECT_EQ(copiedViewVector.size(), 1u);
    EXPECT_EQ(copiedViewVector.capacity(), 1u);
    EXPECT_EQ(viewVector.data(), copiedViewVector.data());  // Shallow copy
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, MoveConstruction)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector("LifetimeVector", 1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);

    auto movedViewVector(std::move(viewVector));
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
    EXPECT_EQ(movedViewVector.size(), 1u);      // A moved-from View behaves identically to a copied-from View
    EXPECT_EQ(movedViewVector.capacity(), 1u);
    EXPECT_EQ(viewVector.data(), movedViewVector.data());  // Shallow copy
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, CopyAssignment)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVectorA("LifetimeVectorA", 0);
    EXPECT_EQ(viewVectorA.size(), 0u);
    EXPECT_EQ(viewVectorA.capacity(), 0u);
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVectorB("LifetimeVectorB", 1);
    EXPECT_EQ(viewVectorB.size(), 1u);
    EXPECT_EQ(viewVectorB.capacity(), 1u);

    viewVectorA = viewVectorB;

    EXPECT_EQ(viewVectorA.size(), 1u);
    EXPECT_EQ(viewVectorA.capacity(), 1u);
    EXPECT_EQ(viewVectorB.size(), 1u);
    EXPECT_EQ(viewVectorB.capacity(), 1u);
    EXPECT_EQ(viewVectorA.data(), viewVectorB.data());  // Shallow copy
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, MoveAssignment)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVectorA("LifetimeVectorA", 0);
    EXPECT_EQ(viewVectorA.size(), 0u);
    EXPECT_EQ(viewVectorA.capacity(), 0u);
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVectorB("LifetimeVectorB", 1);
    EXPECT_EQ(viewVectorB.size(), 1u);
    EXPECT_EQ(viewVectorB.capacity(), 1u);

    viewVectorA = std::move(viewVectorB);

    EXPECT_EQ(viewVectorA.size(), 1u);
    EXPECT_EQ(viewVectorA.capacity(), 1u);
    EXPECT_EQ(viewVectorB.size(), 0u);
    EXPECT_EQ(viewVectorB.capacity(), 1u);
    EXPECT_EQ(viewVectorA.data(), viewVectorB.data());
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, reserveUpFromEmpty)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.reserve(1);
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     0);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      0);
}

TEST_F(ViewVectorObjectLifetimes, reserveUpFrom1)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.emplace_back(1);
    viewVector.reserve(2);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      2);
}

TEST_F(ViewVectorObjectLifetimes, resizeUpFromEmpty)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.resize(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, resizeUpFrom1)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.emplace_back(1);
    viewVector.resize(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     2);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      3);
}

TEST_F(ViewVectorObjectLifetimes, resizeScaleUpFromEmpty)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.resize_scale(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, resizeScaleUpFrom1)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.emplace_back(1);
    viewVector.resize_scale(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     2);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      3);
}

TEST_F(ViewVectorObjectLifetimes, clear)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector("TestViewVector", 1);
    viewVector.clear();
    EXPECT_EQ(viewVector.size(), 0u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, swap)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVectorA("TestViewVector", 1);
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVectorB("TestViewVector", 2);
    EXPECT_EQ(viewVectorA.size(), 1u);
    EXPECT_EQ(viewVectorA.capacity(), 1u);
    EXPECT_EQ(viewVectorB.size(), 2u);
    EXPECT_EQ(viewVectorB.capacity(), 2u);

    viewVectorA.swap(viewVectorB);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     3);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      3);
}

TEST_F(ViewVectorObjectLifetimes, pushBackRvalue)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.push_back(ObjectLifetimeSpy(1));
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      2);
}

TEST_F(ViewVectorObjectLifetimes, pushBackLvalue)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    ObjectLifetimeSpy spy(1);
    viewVector.push_back(spy);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      2);
}

TEST_F(ViewVectorObjectLifetimes, emplaceBack)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.emplace_back(1);
    EXPECT_EQ(viewVector.size(), 1u);
    EXPECT_EQ(viewVector.capacity(), 1u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     1);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      1);
}

TEST_F(ViewVectorObjectLifetimes, pushBackRvalueTwice_reallocation)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.push_back(ObjectLifetimeSpy(1));
    viewVector.push_back(ObjectLifetimeSpy(2));
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     2);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 3);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      5);
}

TEST_F(ViewVectorObjectLifetimes, pushBackLvalueTwice_reallocation)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    ObjectLifetimeSpy spy1(1);
    viewVector.push_back(spy1);
    ObjectLifetimeSpy spy2(2);
    viewVector.push_back(spy2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     2);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 2);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      5);
}

TEST_F(ViewVectorObjectLifetimes, emplaceBackTwice_reallocation)
{
  {
    stk::mesh::impl::ViewVector<ObjectLifetimeSpy> viewVector;
    viewVector.emplace_back(1);
    viewVector.emplace_back(2);
    EXPECT_EQ(viewVector.size(), 2u);
    EXPECT_EQ(viewVector.capacity(), 2u);
  }
  EXPECT_EQ(objectLifetimeSpy_getNumConstructions(),     2);
  EXPECT_EQ(objectLifetimeSpy_getNumCopyConstructions(), 0);
  EXPECT_EQ(objectLifetimeSpy_getNumMoveConstructions(), 1);
  EXPECT_EQ(objectLifetimeSpy_getNumDestructions(),      3);
}


}
