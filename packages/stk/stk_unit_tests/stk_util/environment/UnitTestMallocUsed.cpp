// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stddef.h>                     // for size_t, ptrdiff_t
#include <stdlib.h>                     // for free, malloc, rand, realloc
#include <string.h>                     // for memset
#include <iomanip>                      // for operator<<, setw
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_util/environment/CPUTime.hpp>  // for cpu_time
#include <gtest/gtest.h>

namespace
{

TEST(UnitTestMallocUsed, Malloc_1_8)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 8;

  size_t start = malloc_used();

  size_t used = 0;

  char *x = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

TEST(UnitTestMallocUsed, Malloc_1_16)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 16;

  size_t start = malloc_used();

  size_t used = 0;

  char *x = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

TEST(UnitTestMallocUsed, Malloc_1_32)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 32;

  size_t start = malloc_used();

  size_t used = 0;

  char *x = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

TEST(UnitTestMallocUsed, Malloc_1_1024)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 1024;

  size_t start = malloc_used();

  size_t used = 0;

  char *x = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

TEST(UnitTestMallocUsed, Malloc_100x1024)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 1024;

  for (size_t i = 0; i != 100; ++i)  {

    size_t start = malloc_used();

    size_t used = 0;

    char *x = static_cast<char *>(malloc(bytes_to_allocate));

    used += malloc_used();

    free(x);

    size_t end = malloc_used();

    EXPECT_LE(bytes_to_allocate, used - start);
    EXPECT_LE(start, end);
    std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
  }
#endif
}

TEST(UnitTestMallocUsed, Malloc_1_1M)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 1024*1024;

  size_t start = malloc_used();

  size_t used = 0;

  char *x = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

TEST(UnitTestMallocUsed, Malloc_1_100M)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 100*1024*1024;

  size_t start = malloc_used();

  size_t used = 0;

  char *x = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

TEST(UnitTestMallocUsed, Malloc_100_32)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 32;

  size_t start = malloc_used();

  char *x[100];
  size_t used = 0;

  for (size_t i = 0; i < 100; ++i)
    x[i] = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  for (size_t i = 0; i < 100; ++i)
    free(x[i]);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate*100, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

TEST(UnitTestMallocUsed, Malloc_100_1024)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 1024;

  size_t start = malloc_used();

  char *x[100];
  size_t used = 0;

  for (size_t i = 0; i < 100; ++i)
    x[i] = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  for (size_t i = 0; i < 100; ++i)
    free(x[i]);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate*100, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

TEST(UnitTestMallocUsed, Malloc_100_1M)
{
#if defined SIERRA_PTMALLOC3_ALLOCATOR || defined SIERRA_PTMALLOC2_ALLOCATOR
  static const size_t bytes_to_allocate = 1024*1024;

  size_t start = malloc_used();

  char *x[100];
  size_t used = 0;

  for (size_t i = 0; i < 100; ++i)
    x[i] = static_cast<char *>(malloc(bytes_to_allocate));

  used += malloc_used();

  for (size_t i = 0; i < 100; ++i)
    free(x[i]);

  size_t end = malloc_used();

  EXPECT_LE(bytes_to_allocate*100, used - start);
  EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

}
