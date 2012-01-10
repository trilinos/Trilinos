/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_iostream.hpp>
#include <fei_Pattern.hpp>

#include <vector>
#include <cmath>

TEUCHOS_UNIT_TEST(Pattern, Pattern_test1)
{
  int numIDs = 6;
  std::vector<int> idTypes(numIDs);
  std::vector<snl_fei::RecordCollection*> recColls(numIDs,(snl_fei::RecordCollection*)NULL);
  std::vector<int> fieldsPerID(numIDs);
  std::vector<int> fieldIDs(3);
  std::vector<int> fieldSizes(3, 1);

  idTypes[0] = 0;
  idTypes[1] = 0;
  idTypes[2] = 0;
  idTypes[3] = 1;
  idTypes[4] = 1;
  idTypes[5] = 1;

  fieldsPerID[0] = 0;
  fieldsPerID[1] = 0;
  fieldsPerID[2] = 0;
  fieldsPerID[3] = 1;
  fieldsPerID[4] = 1;
  fieldsPerID[5] = 1;

  fieldIDs[0] = 0;
  fieldIDs[1] = 0;
  fieldIDs[2] = 0;

  fei::Pattern pattern1(numIDs, 0, recColls[0], &fieldsPerID[0], &fieldIDs[0], &fieldSizes[0]);

  TEUCHOS_TEST_EQUALITY(pattern1.getTotalNumFields(), 3, out, success);
  TEUCHOS_TEST_EQUALITY(pattern1.getNumIndices(), 3, out, success);

  fei::Pattern pattern2(numIDs, &idTypes[0], &recColls[0], &fieldsPerID[0], &fieldIDs[0], &fieldSizes[0]);

  TEUCHOS_TEST_EQUALITY(pattern2.getTotalNumFields(), 3, out, success);
  TEUCHOS_TEST_EQUALITY(pattern2.getNumIndices(), 3, out, success);
}

