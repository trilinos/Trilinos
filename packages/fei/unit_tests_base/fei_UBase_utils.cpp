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
#include <fei_ArrayUtils.hpp>

#include <vector>
#include <cmath>

TEUCHOS_UNIT_TEST(fei_utils, insertion_sort_with_companions)
{
  int len = 5;
  std::vector<int> iarray(len);
  std::vector<double> darray(len);

  iarray[0] = 2;
  iarray[1] = 3;
  iarray[2] = 0;
  iarray[3] = 4;
  iarray[4] = 1;

  darray[0] = 2.0;
  darray[1] = 3.0;
  darray[2] = 0.0;
  darray[3] = 4.0;
  darray[4] = 1.0;

  fei::insertion_sort_with_companions(len, &iarray[0], &darray[0]);

  for(int i=0; i<len; ++i) {
    TEUCHOS_TEST_EQUALITY(iarray[i], i, out, success);

    TEUCHOS_TEST_EQUALITY(std::abs(darray[i] - 1.0*i) < 1.e-49, true, out, success);
  }

  iarray.resize(20);

  len = 4;

  iarray[10] = 91;
  iarray[11] = 2225;
  iarray[12] = 214;
  iarray[13] = 3;

  fei::insertion_sort_with_companions(len, &iarray[10], &darray[0]);
}

