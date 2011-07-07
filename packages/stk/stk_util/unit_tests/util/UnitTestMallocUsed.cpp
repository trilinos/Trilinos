/*------------------------------------------------------------------------*/
/*                 Copyright 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdlib.h>

#include <stk_util/util/MallocUsed.h>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_1_8)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 8;
    
  size_t start = malloc_used();
    
  size_t used = 0;

  char *x = (char *) malloc(bytes_to_allocate);

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_1_16)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 16;
    
  size_t start = malloc_used();
    
  size_t used = 0;

  char *x = (char *) malloc(bytes_to_allocate);

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_1_32)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 32;
    
  size_t start = malloc_used();
    
  size_t used = 0;

  char *x = (char *) malloc(bytes_to_allocate);

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_1_1024)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 1024;
    
  size_t start = malloc_used();
    
  size_t used = 0;

  char *x = (char *) malloc(bytes_to_allocate);

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_100x1024)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 1024;

  for (size_t i = 0; i != 100; ++i)  {
    
    size_t start = malloc_used();
    
    size_t used = 0;
    
    char *x = (char *) malloc(bytes_to_allocate);
    
    used += malloc_used();
    
    free(x);
    
    size_t end = malloc_used();
  
    STKUNIT_EXPECT_LE(bytes_to_allocate, used - start);
    STKUNIT_EXPECT_LE(start, end);
    std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
  }
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_1_1M)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 1024*1024;
    
  size_t start = malloc_used();
    
  size_t used = 0;

  char *x = (char *) malloc(bytes_to_allocate);

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_1_100M)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 100*1024*1024;
    
  size_t start = malloc_used();
    
  size_t used = 0;

  char *x = (char *) malloc(bytes_to_allocate);

  used += malloc_used();

  free(x);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_100_32)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 32;

  size_t start = malloc_used();

  char *x[100];
  size_t used = 0;
    
  for (size_t i = 0; i < 100; ++i)
    x[i] = (char *) malloc(bytes_to_allocate);
    
  used += malloc_used();

  for (size_t i = 0; i < 100; ++i)
    free(x[i]);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate*100, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_100_1024)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 1024;

  size_t start = malloc_used();

  char *x[100];
  size_t used = 0;
    
  for (size_t i = 0; i < 100; ++i)
    x[i] = (char *) malloc(bytes_to_allocate);
    
  used += malloc_used();

  for (size_t i = 0; i < 100; ++i)
    free(x[i]);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate*100, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}

STKUNIT_UNIT_TEST(UnitTestMallocUsed, Malloc_100_1M)
{
#ifdef SIERRA_PTMALLOC3_ALLOCATOR
  static const size_t bytes_to_allocate = 1024*1024;

  size_t start = malloc_used();

  char *x[100];
  size_t used = 0;
    
  for (size_t i = 0; i < 100; ++i)
    x[i] = (char *) malloc(bytes_to_allocate);
    
  used += malloc_used();

  for (size_t i = 0; i < 100; ++i)
    free(x[i]);

  size_t end = malloc_used();

  STKUNIT_EXPECT_LE(bytes_to_allocate*100, used - start);
  STKUNIT_EXPECT_LE(start, end);
  std::cout << "start " << start << ", end " << end << ", used " << used - start << std::endl;
#endif
}
