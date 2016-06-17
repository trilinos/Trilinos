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
// Test for Zoltan2::InputTraitsBad

#define INVERT_STATIC_ASSERT_FOR_UNIT_TESTING  // This should only be used for this test - varifies that bad forms would not compile

#include <Zoltan2_InputTraits.hpp>

using Zoltan2::InputTraits;

#define BEGIN_CHECK   if( InputTraits<
#define END_CHECK     >::name() == "unused string" ) { std::cout << "Unused - Just make sure it compiles." << std::endl; }

template< class A, class B, class C >
class SomeBadType {
public: SomeBadType() {}
};

int main(int argc, char *argv[])
{
  // makes a general list of 'bad' types - generally just check these all would fail - the #define above inverts the status of the static_asserts
  // that means these should all compile but only in this test and nowhere else
  BEGIN_CHECK    Zoltan2::BasicUserTypes<long, int, long long>                          END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<long long, long, unsigned int>                 END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<unsigned int, unsigned int, unsigned long>     END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<unsigned long, unsigned long, long>            END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<long long, long long, unsigned int>            END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<double, float, long long>                            END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, long, float>                                  END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, double, unsigned long>                        END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, unsigned long, double>                        END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<long long, long long, unsigned int>                  END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<SomeBadType<int, int, int>, int, long long>          END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<int, SomeBadType<int, double, float>>                 END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<SomeBadType<int, long long, int>, unsigned int>       END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<int, double>                                          END_CHECK

  // set the PASS keyword
  std::cout << "Validated bad InputTraits - The test PASSED because it compiled with the static_assert checks inverted." << std::endl;
  return 0;
}

