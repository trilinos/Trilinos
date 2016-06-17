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
#define END_CHECK     >::name() == "unused string" ) { std::cout << "Unused - Just make sure it compiles." << std::endl; }

int main(int argc, char *argv[])
{
  // validate Zoltan2::BasicUserTypes
  BEGIN_CHECK    Zoltan2::BasicUserTypes<double, int, long long>                  END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, long, unsigned int>               END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, unsigned int, unsigned long>      END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<float, unsigned long, long>              END_CHECK
  BEGIN_CHECK    Zoltan2::BasicUserTypes<double, long long, unsigned int>         END_CHECK

  // validate Xpetra::CrsMatrix
  BEGIN_CHECK    Xpetra::CrsMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, long, unsigned int>                     END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, unsigned int, unsigned long>            END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<float, unsigned long, long>                    END_CHECK
  BEGIN_CHECK    Xpetra::CrsMatrix<double, long long, unsigned int>               END_CHECK

  // validate Tpetra::CrsMatrix
  BEGIN_CHECK    Tpetra::CrsMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<float, long, unsigned int>                     END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<float, unsigned int, unsigned long>            END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<float, unsigned long, long>                    END_CHECK
  BEGIN_CHECK    Tpetra::CrsMatrix<double, long long, unsigned int>               END_CHECK

  // validate Xpetra::RowMatrix
  BEGIN_CHECK    Xpetra::RowMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Xpetra::RowMatrix<float, long, unsigned int>                     END_CHECK
  BEGIN_CHECK    Xpetra::RowMatrix<float, unsigned int, unsigned long>            END_CHECK
  BEGIN_CHECK    Xpetra::RowMatrix<float, unsigned long, long>                    END_CHECK
  BEGIN_CHECK    Xpetra::RowMatrix<double, long long, unsigned int>               END_CHECK

  // validate Tpetra::RowMatrix
  BEGIN_CHECK    Tpetra::RowMatrix<double, int, long long>                        END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<float, long, unsigned int>                     END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<float, unsigned int, unsigned long>            END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<float, unsigned long, long>                    END_CHECK
  BEGIN_CHECK    Tpetra::RowMatrix<double, long long, unsigned int>               END_CHECK

  // validate Xpetra::Vector
  BEGIN_CHECK    Xpetra::Vector<double, int, long long>                           END_CHECK
  BEGIN_CHECK    Xpetra::Vector<float, long, unsigned int>                        END_CHECK
  BEGIN_CHECK    Xpetra::Vector<float, unsigned int, unsigned long>               END_CHECK
  BEGIN_CHECK    Xpetra::Vector<float, unsigned long, long>                       END_CHECK
  BEGIN_CHECK    Xpetra::Vector<double, long long, unsigned int>                  END_CHECK

  // validate Tpetra::RowMatrix
  BEGIN_CHECK    Tpetra::Vector<double, int, long long>                           END_CHECK
  BEGIN_CHECK    Tpetra::Vector<float, long, unsigned int>                        END_CHECK
  BEGIN_CHECK    Tpetra::Vector<float, unsigned int, unsigned long>               END_CHECK
  BEGIN_CHECK    Tpetra::Vector<float, unsigned long, long>                       END_CHECK
  BEGIN_CHECK    Tpetra::Vector<double, long long, unsigned int>                  END_CHECK

  // validate Xpetra::MultiVector
  BEGIN_CHECK    Xpetra::MultiVector<double, int, long long>                      END_CHECK
  BEGIN_CHECK    Xpetra::MultiVector<float, long, unsigned int>                   END_CHECK
  BEGIN_CHECK    Xpetra::MultiVector<float, unsigned int, unsigned long>          END_CHECK
  BEGIN_CHECK    Xpetra::MultiVector<float, unsigned long, long>                  END_CHECK
  BEGIN_CHECK    Xpetra::MultiVector<double, long long, unsigned int>             END_CHECK

  // validate Tpetra::MultiVector
  BEGIN_CHECK    Tpetra::MultiVector<double, int, long long>                      END_CHECK
  BEGIN_CHECK    Tpetra::MultiVector<float, long, unsigned int>                   END_CHECK
  BEGIN_CHECK    Tpetra::MultiVector<float, unsigned int, unsigned long>          END_CHECK
  BEGIN_CHECK    Tpetra::MultiVector<float, unsigned long, long>                  END_CHECK
  BEGIN_CHECK    Tpetra::MultiVector<double, long long, unsigned int>             END_CHECK

  // validate Xpetra::CrsGraph
  BEGIN_CHECK    Xpetra::CrsGraph<int, long long>                                 END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<long, unsigned int>                             END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<unsigned int, unsigned long>                    END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<unsigned long, long>                            END_CHECK
  BEGIN_CHECK    Xpetra::CrsGraph<long long, unsigned int>                        END_CHECK

  // validate Tpetra::CrsGraph
  BEGIN_CHECK    Tpetra::CrsGraph<int, long long>                                 END_CHECK
  BEGIN_CHECK    Tpetra::CrsGraph<long, unsigned int>                             END_CHECK
  BEGIN_CHECK    Tpetra::CrsGraph<unsigned int, unsigned long>                    END_CHECK
  BEGIN_CHECK    Tpetra::CrsGraph<unsigned long, long>                            END_CHECK
  BEGIN_CHECK    Tpetra::CrsGraph<long long, unsigned int>                        END_CHECK

  // set the PASS keyword
  std::cout << "Validated InputTraits - The test PASSED because it compiled" << std::endl;
  return 0;
}

