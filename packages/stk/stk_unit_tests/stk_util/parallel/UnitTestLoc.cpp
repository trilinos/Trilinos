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
#include "stk_util/parallel/MPI.hpp"  // for Loc, create_Loc, operator==
#include <cstdint>                    // for int64_t
#include <utility>                    // for move

namespace sierra
{
namespace MPI
{

TEST(Loc, Constructors)
{
  Loc<double, int64_t> l;
  l.m_value = 1.;
  l.m_loc = 3;

  Loc<double, int64_t> copy(l);
  EXPECT_DOUBLE_EQ(l.m_value, copy.m_value);
  EXPECT_EQ(l.m_loc, copy.m_loc);

  Loc<double, int64_t> move(std::move(l));
  EXPECT_DOUBLE_EQ(l.m_value, move.m_value);
  EXPECT_EQ(l.m_loc, move.m_loc);
}

TEST(Loc, Assignment)
{
  Loc<double, int64_t> l;
  l.m_value = 1.;
  l.m_loc = 3;

  Loc<double, int64_t> assign;
  assign = l;
  EXPECT_DOUBLE_EQ(l.m_value, assign.m_value);
  EXPECT_EQ(l.m_loc, assign.m_loc);

  // Assignment from a dereferenced volatile * is required for use of Loc with
  // Kokkos::atomic_compare_exchange();
  volatile Loc<double, int64_t> * const volatile_ptr = &l;
  Loc<double, int64_t> assign_from_volatile;
  assign_from_volatile = *volatile_ptr;
  EXPECT_DOUBLE_EQ(l.m_value, assign_from_volatile.m_value);
  EXPECT_EQ(l.m_loc, assign_from_volatile.m_loc);

  Loc<double, int64_t> assign_to_volatile;
  volatile Loc<double, int64_t> * const to_volatile_ptr = &assign_to_volatile;
  *to_volatile_ptr = l;
  EXPECT_DOUBLE_EQ(l.m_value, assign_to_volatile.m_value);
  EXPECT_EQ(l.m_loc, assign_to_volatile.m_loc);
}

TEST(Loc, CompareEqual)
{
  double lhs_val = 1.;
  int64_t lhs_id = 1;
  auto lhs = create_Loc(lhs_val, lhs_id);

  EXPECT_TRUE(lhs == lhs);

  {
    double rhs_val = 5.;
    int64_t rhs_id = 1;
    auto rhs = create_Loc(rhs_val, rhs_id);

    EXPECT_FALSE(lhs == rhs);
  }

  {
    double rhs_val = 1.;
    int64_t rhs_id = 2;
    auto rhs = create_Loc(rhs_val, rhs_id);

    EXPECT_FALSE(lhs == rhs);
  }
}

}
}
