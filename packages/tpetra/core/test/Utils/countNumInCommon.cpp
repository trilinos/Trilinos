/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_Util.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <algorithm>
#include <vector>

namespace { // (anonymous)

  void
  testLists(Teuchos::FancyOStream& out,
            bool& success,
            std::vector<int>& list1,
            std::vector<int>& list2,
            const size_t expectedNumInCommon,
            const size_t expectedUnionSize)
  {
    std::sort(list1.begin(), list1.end());
    std::sort(list2.begin(), list2.end());
    auto iter1 = std::unique(list1.begin(), list1.end());
    auto iter2 = std::unique(list2.begin(), list2.end());

    using Tpetra::Details::countNumInCommon;
    const size_t numInCommon =
      countNumInCommon(list1.begin(), iter1,
                       list2.begin(), iter2);
    TEST_EQUALITY( numInCommon, expectedNumInCommon );

    const size_t unionSize =
      list1.size() + list2.size() - numInCommon;
    TEST_EQUALITY( unionSize, expectedUnionSize );

    std::vector<int> unionInds(list1.size() + list2.size());
    TEST_EQUALITY( unionInds.size(), list1.size() + list2.size() );

    auto unionIter = std::set_union(list1.begin(), iter1,
                                    list2.begin(), iter2,
                                    unionInds.begin());
    const size_t unionSize2 = size_t(unionIter - unionInds.begin());
    TEST_EQUALITY( unionSize2, expectedUnionSize );
  }

  TEUCHOS_UNIT_TEST( Utils, CountNumInCommon_short )
  {
    {
      std::vector<int> list1;
      std::vector<int> list2;
      testLists(out, success, list1, list2, 0, 0);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{666}};
      std::vector<int> list2;
      testLists(out, success, list1, list2, 0, 1);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1;
      std::vector<int> list2{{666}};
      testLists(out, success, list1, list2, 0, 1);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{418}};
      std::vector<int> list2{{418}};
      testLists(out, success, list1, list2, 1, 1);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{418, 419}};
      std::vector<int> list2{{418}};
      testLists(out, success, list1, list2, 1, 2);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{418}};
      std::vector<int> list2{{418, 419}};
      testLists(out, success, list1, list2, 1, 2);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{417, 418}};
      std::vector<int> list2{{418}};
      testLists(out, success, list1, list2, 1, 2);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{418}};
      std::vector<int> list2{{417, 418}};
      testLists(out, success, list1, list2, 1, 2);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{417, 418, 419}};
      std::vector<int> list2{{418}};
      testLists(out, success, list1, list2, 1, 3);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{418}};
      std::vector<int> list2{{417, 418, 419}};
      testLists(out, success, list1, list2, 1, 3);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{415, 418, 421}};
      std::vector<int> list2{{418}};
      testLists(out, success, list1, list2, 1, 3);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{418}};
      std::vector<int> list2{{415, 418, 421}};
      testLists(out, success, list1, list2, 1, 3);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{415, 419, 421}};
      std::vector<int> list2{{418}};
      testLists(out, success, list1, list2, 0, 4);
      if (! success) {
        return;
      }
    }
    {
      std::vector<int> list1{{418}};
      std::vector<int> list2{{415, 419, 421}};
      testLists(out, success, list1, list2, 0, 4);
      if (! success) {
        return;
      }
    }
  }

  // This came from an application test problem.
  TEUCHOS_UNIT_TEST( Utils, CountNumInCommon_app )
  {
    std::vector<int> newGblColInds {{142944, 142945, 142946, 142947, 142948, 142949, 142950, 142951, 142952, 142953, 142954, 142955, 142959, 142960, 142961, 142965, 142966, 142967, 142968, 142969, 142970, 143142, 143143, 143144, 198279, 198280, 198281, 198282, 198283, 198284, 198291, 198292, 198293, 198303, 198304, 198305, 198309, 198310, 198311, 198333, 198334, 198335, 198336, 198337, 198338, 198339, 198340, 198341, 198342, 198343, 198344, 198345, 198346, 198347, 198348, 198349, 198350, 198351, 198352, 198353, 198354, 198355, 198356, 198699, 198700, 198701, 198702, 198703, 198704, 198705, 198706, 198707, 198708, 198709, 198710, 198711, 198712, 198713, 198729, 198730, 198731, 198732, 198733, 198734, 198735, 198736, 198737, 198738, 198739, 198740, 198741, 198742, 198743, 198744, 198745, 198746}};

    std::vector<int> curGblColInds {{166215, 166216, 166217, 166218, 166219, 166220, 166221, 166222, 166223, 166224, 166225, 166226, 166227, 166228, 166229, 166230, 166231, 166232, 166233, 166234, 166235, 166236, 166237, 166238, 166239, 166240, 166241, 166242, 166243, 166244, 166245, 166246, 166247, 198279, 198280, 198281, 198282, 198283, 198284, 198285, 198286, 198287, 198288, 198289, 198290, 198291, 198292, 198293, 198294, 198295, 198296, 198297, 198298, 198299, 198300, 198301, 198302, 198303, 198304, 198305, 198306, 198307, 198308, 198309, 198310, 198311, 198312, 198313, 198314, 198315, 198316, 198317, 198333, 198334, 198335, 198336, 198337, 198338, 198339, 198340, 198341, 198342, 198343, 198344, 198345, 198346, 198347, 198348, 198349, 198350, 198351, 198352, 198353, 198354, 198355, 198356}};

    constexpr size_t newGblColIndsSize (96);
    constexpr size_t curGblColIndsSize (96);

    TEST_EQUALITY( newGblColInds.size(), newGblColIndsSize );
    TEST_EQUALITY( curGblColInds.size(), curGblColIndsSize );

    std::sort(newGblColInds.begin(), newGblColInds.end());
    std::sort(curGblColInds.begin(), curGblColInds.end());

    TEST_EQUALITY( newGblColInds.size(), newGblColIndsSize );
    TEST_EQUALITY( curGblColInds.size(), curGblColIndsSize );

    auto newIter = std::unique(newGblColInds.begin(), newGblColInds.end());
    auto curIter = std::unique(curGblColInds.begin(), curGblColInds.end());

    TEST_EQUALITY( size_t(newIter - newGblColInds.begin()), newGblColIndsSize );
    TEST_EQUALITY( size_t(curIter - curGblColInds.begin()), curGblColIndsSize );

    constexpr size_t expectedMergeSize = 153;
    constexpr size_t expectedNumInCommon = newGblColIndsSize +
      curGblColIndsSize - expectedMergeSize;
    testLists(out, success, newGblColInds, curGblColInds,
              expectedNumInCommon, expectedMergeSize);
  }

} // namespace (anonymous)
