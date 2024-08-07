// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Test for Zoltan2::InputTraitsBad

// This should only be used for this test - verifies bad forms do not compile
#define Z2_INVERT_STATIC_ASSERT_FOR_UNIT_TESTING

#include <Zoltan2_InputTraits.hpp>

using Zoltan2::InputTraits;

#define BEGIN_CHECK   if( InputTraits<
#define END_CHECK     >::name() == "unused string" ) \
  { std::cout << "Unused - Just make sure it compiles." << std::endl; }

template< class A, class B, class C >
class SomeBadType {
public: SomeBadType() {}
};

int main(int argc, char *argv[])
{
  // makes a general list of 'bad' types - generally just check these all
  // would fail - the #define above inverts the status of the static_asserts
  // that means these should all compile but only in this test and nowhere else

  // scalar ordinal (first slot) must be float, double, int, or long long
  // this validates we would fail for any
  BEGIN_CHECK    Zoltan2::BasicUserTypes<unsigned int, int, long>                       END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<unsigned long, int, long>                      END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<uint32_t, int, long>                           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<unsigned long long, int, long>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<uint64_t, int, long>                           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<size_t, int, long>                             END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<std::complex<int>, int, long>                  END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<std::complex<double>, int, long>               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<std::complex<float>, int, long>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<SomeBadType<int,long, int>, int, long>         END_CHECK

  // local ordinal (second slot) must always be signed
  // this validates that we would fail for any unsigned setting
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, unsigned int, long>                     END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, unsigned long, long>                    END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, uint32_t, long>                         END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, unsigned long long, long>               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, uint64_t, long>                         END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, size_t, long>                           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, float, long>                            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, double, long>                           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, std::complex<int>, long>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, std::complex<float>, long>              END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, std::complex<double>, long>             END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, SomeBadType<int,long, int>, long>       END_CHECK

  // global ordinal (third slot) must be int, long, long long, ssize_t,
  // unsigned int, unsigned long, unsigned long long, size_t
  BEGIN_CHECK    Zoltan2::BasicUserTypes<int, int, std::complex<int>>                   END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<int, int, std::complex<float>>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<int, int, std::complex<double>>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<int, int, SomeBadType<int,long, int>>          END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<int, int, float>                               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<int, int, double>                              END_CHECK

  // we also want to make sure that user typedefs will work just fine
  // lots of redundancy here  to make sure std::is_same is  working as we think
  typedef signed int user_int_t;
  typedef unsigned int user_unsigned_int_t;
  typedef float user_float_t;
  typedef double user_double_t;
  typedef uint32_t user_uint32_t;
  typedef uint64_t user_uint64_t;
  typedef signed long user_long_t;
  typedef unsigned long user_unsigned_long_t;
  typedef unsigned long long user_unsigned_long_long_t;
  typedef size_t user_size_t;
  typedef SomeBadType<int,int, int> user_some_bad_t;

  // scalar ordinal (first slot) must be float, double, or int
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_unsigned_int_t, user_int_t, user_long_t>                                      END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_unsigned_long_t, user_int_t, user_long_t>                                     END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_uint32_t, user_int_t, user_long_t>                                            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_unsigned_long_long_t, user_int_t, user_long_t>                                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_uint64_t, user_int_t, user_long_t>                                            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_size_t, user_int_t, user_long_t>                                              END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<std::complex<user_int_t>, user_int_t, user_long_t>                                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<SomeBadType<user_int_t,user_long_t, user_int_t>, user_int_t, user_long_t>          END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_uint32_t, user_int_t, user_long_t>                                            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_some_bad_t, user_int_t, user_long_t>                                          END_CHECK

  // local ordinal (second slot) must always be signed
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_unsigned_int_t, user_long_t>                                    END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_unsigned_long_t, user_long_t>                                   END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_uint32_t, user_long_t>                                          END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_unsigned_long_long_t, user_long_t>                              END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_uint64_t, user_long_t>                                          END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_size_t, user_long_t>                                            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_float_t, user_long_t>                                           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_double_t, user_long_t>                                          END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, std::complex<user_int_t>, user_long_t>                               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, SomeBadType<user_int_t,user_long_t, user_int_t>, user_long_t>        END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, size_t, user_long_t>        END_CHECK

  // global ordinal (third slot) must be int, long, long long, ssize_t,
  // unsigned int, unsigned long, unsigned long long, size_t
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_int_t, user_int_t, std::complex<user_int_t>>                                  END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_int_t, user_int_t, SomeBadType<user_int_t,user_long_t, user_int_t>>           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_int_t, user_int_t, user_float_t>                                              END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_int_t, user_int_t, user_double_t>                                             END_CHECK


  // some more checks that should all fail - this doesn't check all
  // combinations but just tries a bunch of things on different class types
  BEGIN_CHECK    Zoltan2::BasicUserTypes<unsigned int, unsigned int, unsigned long>     END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<unsigned long, unsigned long, long>            END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<double, float, long long>                            END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, long, float>                                  END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, double, unsigned long>                        END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, unsigned long, double>                        END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<SomeBadType<int, int, int>, int, long long>          END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<int, SomeBadType<int, double, float>>                 END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<SomeBadType<int, long long, int>, unsigned int>       END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<int, double>                                          END_CHECK

  // set the PASS keyword
  std::cout << "Validated bad InputTraits - The test PASSED because it "
    "compiled with the static_assert checks inverted." << std::endl;
  return 0;
}

