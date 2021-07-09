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
// This is another test where you need to look at the output to
// know if it's right.  This should be fixed.

/*! \file StridedData.cpp
 *  \brief Tests the StridedData class.
 *
 *   \todo Some of the tests require that you look at the output
 *          to know if they did the right thing.  Enhance this so
 *          the test itself determines correctness.
 */

#include <Zoltan2_StridedData.hpp>   
#include <Zoltan2_TestHelpers.hpp>   

using Zoltan2::StridedData;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::Array;

void StridedDataTest(const Teuchos::SerialComm<int> &comm)
{
  // StridedData template arguments

  typedef int     index_t;
  typedef double  value_t;
  typedef float   different_value_t;

  typedef StridedData<index_t, value_t> stridedInput_t;
  bool aok = true;

  /*! \test strided input with stride 1
   */

  ArrayRCP<value_t> input1(new value_t [12], 0, 12, true);
  for (int i=0; i < 12; i++)
    input1[i] = (i+1) * 5; 

  RCP<stridedInput_t> s1;

  try{
    s1 = rcp<stridedInput_t>(new stridedInput_t(input1, 1));
  }
  catch (std::exception &e){
    aok = false;
  }
  TEST_FAIL_AND_EXIT(comm, aok, "Error in constructor 1", 1);

  std::cout << std::endl;
  std::cout << "Test 1, input: " << input1 << std::endl;
  std::cout << "[] test: ";
  for (int i=0; i < 12; i++)
    std::cout << (*s1)[i] << " ";
  std::cout << std::endl;

  ArrayRCP<const value_t> fromS1;
  s1->getInputArray(fromS1);
  std::cout << "getInputArray test: ";
  for (int i=0; i < 12; i++)
    std::cout << fromS1[i] << " ";
  std::cout << std::endl;

  stridedInput_t s1Copy;
  s1Copy = *s1;

  std::cout << "assignment operator test: ";
  for (int i=0; i < 12; i++)
    std::cout << s1Copy[i] << " ";
  std::cout << std::endl;

  // test a different type
  ArrayRCP<const different_value_t> fromS1too;
  s1->getInputArray(fromS1too);
  std::cout << "getInputArray test -- different type: ";
  for (int i=0; i < 12; i++)
    std::cout << fromS1too[i] << " ";
  std::cout << std::endl;

  /*! \test strided input with stride 3
   */

  ArrayRCP<value_t> input2(new value_t [12], 0, 12, true);
  for (int i=0; i < 12; i+=3)
    input2[i] = (i+1) * -5.0; 

  RCP<stridedInput_t> s2;

  try{
    s2 = rcp<stridedInput_t>(new stridedInput_t(input2, 3));
  }
  catch (std::exception &e){
    aok = false;
  }
  TEST_FAIL_AND_EXIT(comm, aok, "Error in constructor 2", 2);

  std::cout << std::endl;
  std::cout << "Test 2, input: " << input2 << std::endl;
  std::cout << "[] test: ";
  for (int i=0; i < 4; i++)
    std::cout << (*s2)[i] << " ";
  std::cout << std::endl;

  ArrayRCP<const value_t> fromS2;
  s2->getInputArray(fromS2);
  std::cout << "getInputArray test: ";
  for (int i=0; i < 4; i++)
    std::cout << fromS2[i] << " ";
  std::cout << std::endl;

  stridedInput_t s2Copy;
  s2Copy = *s2;

  std::cout << "assignment operator test: ";
  for (int i=0; i < 4; i++)
    std::cout << s2Copy[i] << " ";
  std::cout << std::endl;

  // test a different type
  ArrayRCP<const different_value_t> fromS2too;
  s1->getInputArray(fromS2too);
  std::cout << "getInputArray test -- different type: ";
  for (int i=0; i < 4; i++)
    std::cout << fromS2too[i] << " ";
  std::cout << std::endl;
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > tcomm = Tpetra::getDefaultComm();

  // Run the test on only one rank. 
  // There's no parallelism involved in StridedData, 
  // and the output is neater on only one proc.
  if (tcomm->getRank() > 0)
    return 0;

  Teuchos::SerialComm<int> comm;

  StridedDataTest(comm);

  std::cout << "PASS" << std::endl;
}
