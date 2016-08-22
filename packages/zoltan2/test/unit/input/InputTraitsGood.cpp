// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Test for Zoltan2::InputTraitsGood

#include <Zoltan2_InputTraits.hpp>

using Zoltan2::InputTraits;

#define BEGIN_CHECK   if( InputTraits<
#define END_CHECK     >::name() == "unused string" ) \
  { std::cout << "Unused - Just make sure it compiles." << std::endl; }

int main(int argc, char *argv[])
{
  // scalar ordinal (first slot) must be float, double, or int
  // this validates they all work
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, signed int, long>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<double, long>                            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<int, int, long>                          END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<int32_t, int, long>                      END_CHECK

  // local ordinal (second slot) must always be signed
  // this validates they all work
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, signed int, long>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, signed long, long>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, signed long long, long>           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, ssize_t, long>                    END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int32_t, long>                    END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int64_t, long>                    END_CHECK

  // global ordinal (third slot) must be int, long, long long, ssize_t,
  // unsigned int, unsigned long, unsigned long long, size_t
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, int>                         END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, long>                        END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, long long>                   END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, ssize_t>                     END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, unsigned int>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, unsigned long>               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, unsigned long long>          END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, size_t>                      END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, int32_t>                     END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, int64_t>                     END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, uint32_t>                    END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, int, uint64_t>                    END_CHECK

  // we also want to make sure that user typedefs will work just fine
  // lots of redundancy here  to make sure std::is_same is working as we think
  typedef signed int user_int_t;
  typedef unsigned int user_unsigned_int_t;
  typedef float user_float_t;
  typedef double user_double_t;
  typedef int32_t user_int32_t;
  typedef uint32_t user_uint32_t;
  typedef int64_t user_int64_t;
  typedef uint64_t user_uint64_t;
  typedef signed long user_long_t;
  typedef signed long long user_long_long_t;
  typedef unsigned long user_unsigned_long_t;
  typedef unsigned long long user_unsigned_long_long_t;
  typedef size_t user_size_t;
  typedef ssize_t user_ssize_t;

  // scalar ordinal (first slot) must be float, double, or int
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_long_t>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_double_t, user_long_t>                            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_int_t, user_int_t, user_long_t>                   END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_int32_t, user_int_t, user_long_t>                 END_CHECK

  // local ordinal (second slot) must always be signed
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_long_t>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_long_t, user_long_t>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_long_long_t, user_long_t>           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_ssize_t, user_long_t>               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int32_t, user_long_t>               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int64_t, user_long_t>               END_CHECK

  // global ordinal (third slot) must be int, long, long long, ssize_t,
  // unsigned int, unsigned long, unsigned long long, size_t
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_int_t>                  END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_long_t>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_long_long_t>            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_ssize_t>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_unsigned_int_t>         END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_unsigned_long_t>        END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_unsigned_long_long_t>   END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_size_t>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_int32_t>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_int64_t>                END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_uint32_t>               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<user_float_t, user_int_t, user_uint64_t>               END_CHECK

  // this section takes each specific type and validates a few combinations
  // note that some types may pass only two template parameters, not three,
  // so it can be scalar + local + global or it can be just local + global.

  // validate Zoltan2::BasicUserTypes
  BEGIN_CHECK    Zoltan2::BasicUserTypes<double, signed int, long long>           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<double, int, long long>                  END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, signed long, unsigned int>        END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, signed long long, long>           END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<double, ssize_t, unsigned int>           END_CHECK

  // validate Xpetra::CrsMatrix
  BEGIN_CHECK    Xpetra::CrsMatrix<double, signed int, long long>                 END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, signed long, unsigned int>              END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, signed long long, long>                 END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<double, ssize_t, unsigned int>                 END_CHECK

  // validate Tpetra::CrsMatrix
  BEGIN_CHECK    Tpetra::CrsMatrix<double, signed int, long long>                 END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<float, signed long, unsigned int>              END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<float, signed long long, long>                 END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<double, ssize_t, unsigned int>                 END_CHECK

  // validate Xpetra::RowMatrix
  BEGIN_CHECK    Xpetra::RowMatrix<double, signed int, long long>                 END_CHECK
  BEGIN_CHECK    Xpetra::RowMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Xpetra::RowMatrix<float, signed long, unsigned int>              END_CHECK
  BEGIN_CHECK    Xpetra::RowMatrix<float, signed long long, long>                 END_CHECK
  BEGIN_CHECK    Xpetra::RowMatrix<double, ssize_t, unsigned int>                 END_CHECK

  // validate Tpetra::RowMatrix
  BEGIN_CHECK    Tpetra::RowMatrix<double, signed int, long long>                 END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<float, signed long, unsigned int>              END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<float, signed long long, long>                 END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<double, ssize_t, unsigned int>                 END_CHECK

  // validate Xpetra::Vector
  BEGIN_CHECK    Xpetra::Vector<double, signed int, long long>                    END_CHECK
  BEGIN_CHECK    Xpetra::Vector<double, int, long long>                           END_CHECK
  BEGIN_CHECK    Xpetra::Vector<float, signed long, unsigned int>                 END_CHECK
  BEGIN_CHECK    Xpetra::Vector<float, signed long long, long>                    END_CHECK
  BEGIN_CHECK    Xpetra::Vector<double, ssize_t, unsigned int>                    END_CHECK

  // validate Tpetra::RowMatrix
  BEGIN_CHECK    Tpetra::RowMatrix<double, signed int, long long>                 END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<float, signed long, unsigned int>              END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<float, signed long long, long>                 END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<double, ssize_t, unsigned int>                 END_CHECK

  // validate Xpetra::MultiVector
  BEGIN_CHECK    Xpetra::MultiVector<double, signed int, long long>               END_CHECK
  BEGIN_CHECK    Xpetra::MultiVector<double, int, long long>                      END_CHECK
  BEGIN_CHECK    Xpetra::MultiVector<float, signed long, unsigned int>            END_CHECK
  BEGIN_CHECK    Xpetra::MultiVector<float, signed long long, long>               END_CHECK
  BEGIN_CHECK    Xpetra::MultiVector<double, ssize_t, unsigned int>               END_CHECK

  // validate Tpetra::MultiVector
  BEGIN_CHECK    Tpetra::MultiVector<double, signed int, long long>               END_CHECK
  BEGIN_CHECK    Tpetra::MultiVector<double, int, long long>                      END_CHECK
  BEGIN_CHECK    Tpetra::MultiVector<float, signed long, unsigned int>            END_CHECK
  BEGIN_CHECK    Tpetra::MultiVector<float, signed long long, long>               END_CHECK
  BEGIN_CHECK    Tpetra::MultiVector<double, ssize_t, unsigned int>               END_CHECK

  // validate Xpetra::CrsGraph
  BEGIN_CHECK    Xpetra::CrsGraph<signed int, long long>                          END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<int, long long>                                 END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<signed long, unsigned int>                      END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<signed long long, long>                         END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<ssize_t, unsigned int>                          END_CHECK

  // validate Tpetra::CrsGraph
  BEGIN_CHECK    Tpetra::CrsGraph<signed int, long long>                          END_CHECK
  BEGIN_CHECK    Tpetra::CrsGraph<int, long long>                                 END_CHECK
  BEGIN_CHECK    Tpetra::CrsGraph<signed long, unsigned int>                      END_CHECK
  BEGIN_CHECK    Tpetra::CrsGraph<signed long long, long>                         END_CHECK
  BEGIN_CHECK    Tpetra::CrsGraph<ssize_t, unsigned int>                          END_CHECK

  // set the PASS keyword
  std::cout << "Validated InputTraits - The test PASSED "
    "because it compiled" << std::endl;
  return 0;
}

