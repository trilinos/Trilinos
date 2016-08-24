// @HEADER                                                                                                                                    
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <unistd.h>
#include <vector>
#include "Panzer_MemUtils.hpp"

namespace panzer
{
  // Test that memory usage increases after creating a bunch of variables.
  TEUCHOS_UNIT_TEST(memUtils, beforeAndAfter)
  {
    static const size_t NUM(4096);
    MemUsage before, after, diff;
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
      Teuchos::DefaultComm<int>::getComm();

    // Get the initial memory usage.
    before = getMemoryUsage(*comm);

    // Create a bunch of variables to use up some memory.
    char    plainOldChar[NUM];
    int     plainOldInt[NUM];
    double  plainOldDouble[NUM];
    char*   newChar   = new char[NUM];
    int*    newInt    = new int[NUM];
    double* newDouble = new double[NUM];

    // Get the new memory usage and compute the difference.
    after = getMemoryUsage(*comm);
    diff  = after - before;

    // Test that memory usage has increased.
    TEST_COMPARE(after.currMin, >=, before.currMin);
    TEST_COMPARE(after.currMax, >=, before.currMax);
    TEST_COMPARE(after.currTot, >=, before.currTot);
    TEST_COMPARE(after.peakMin, >=, before.peakMin);
    TEST_COMPARE(after.peakMax, >=, before.peakMax);
    TEST_COMPARE(after.peakTot, >=, before.peakTot);

    // Clean up.
    (void)(plainOldChar);   (void)(newChar);
    (void)(plainOldInt);    (void)(newInt);
    (void)(plainOldDouble); (void)(newDouble);
    delete[] newChar;
    delete[] newInt;
    delete[] newDouble;
  } // end of TEUCHOS_UNIT_TEST()

  // Test that the peak memory usage is always greater than the current memory
  // usage.
  TEUCHOS_UNIT_TEST(memUtils, currentVsPeak)
  {
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
      Teuchos::DefaultComm<int>::getComm();

    // Get the current and peak memory usage.
    MemUsage mem = getMemoryUsage(*comm);

    // Test that the peak usage is greater than the current usage.
    TEST_COMPARE(mem.peakMin, >=, mem.currMin);
    TEST_COMPARE(mem.peakMax, >=, mem.currMax);
    TEST_COMPARE(mem.peakTot, >=, mem.currTot);
  } // end of TEUCHOS_UNIT_TEST()
} // end namespace panzer
