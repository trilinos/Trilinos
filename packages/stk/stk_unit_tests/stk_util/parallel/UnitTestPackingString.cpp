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
#include "stk_util/parallel/Parallel.hpp"      // for parallel_machine_size, MPI_COMM_WORLD, Par...
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBroadcast, CommBuffer
#include <cstdlib>                             // for rand, RAND_MAX
#include <functional>                          // for function
#include <map>                                 // for map
#include <string>                              // for string, basic_string
#include <vector>                              // for vector

namespace
{
std::string random_string(const int len)
{
  static const char alphanum[] = "0123456789"
                                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                 "abcdefghijklmnopqrstuvwxyz";

  std::vector<char> s(len);
  for (int i = 0; i < len; ++i)
  {
    s[i] = alphanum[rand() % (sizeof(alphanum) - 1)];
  }

  std::string str;
  str.assign(s.data(), len);
  return str;
}
double random_double(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}
int random_int(int imin, int imax) { return rand() % (imax - imin) + imin; }

template<class T>
void test_map(std::function<T()> generator) {
  stk::ParallelMachine global = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(global) > 1) return;

  for (int i = 0; i < 100; ++i)
  {
    int num_entries = random_int(1,4);
    std::map<std::string, T> smap;
    for(int j = 0; j < num_entries; ++j) {
      std::string key = random_string(random_int(1, 15));
      T val = generator();
      smap[key] = val;
    }

    stk::CommBroadcast bcast(MPI_COMM_WORLD, 0);
    auto& b = bcast.send_buffer();
    b.pack(smap);
    bcast.allocate_buffer();
    b.pack(smap);

    std::map<std::string, T> rmap;
    b.reset();
    EXPECT_ANY_THROW(b.peek(rmap));
    b.unpack(rmap);

    for(auto&& vp : smap) {
      EXPECT_EQ(vp.second, rmap[vp.first]);
    }
  }
}

} // namespace

TEST(StringPacking, packSendAndPeek)
{
  stk::ParallelMachine global = MPI_COMM_WORLD;
  if (stk::parallel_machine_size(global) > 1) return;

  std::string key = "this is a test";

  stk::CommBroadcast bcast(MPI_COMM_WORLD, 0);
  auto& b = bcast.send_buffer();
  b.pack(key);
  bcast.allocate_buffer();
  b.pack(key);

  b.reset();

  std::string rpeek;
  b.peek(rpeek);
  std::string recv;
  b.unpack(recv);

  EXPECT_EQ(key, rpeek);
  EXPECT_EQ(key, recv);
}

TEST(MapPacking, PackingStringMap)
{
  test_map<std::string>([](){ return random_string(random_int(1, 20));});
}

TEST(MapPacking, PackingDoubleMap)
{
  test_map<double>([](){ return random_double(-1e12, 1e12);});
}

TEST(MapPacking, PackingIntMap)
{
  test_map<int>([](){ return rand();});
}
