// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Details_Merge.hpp"
#include <vector>

namespace { // (anonymous)

using Tpetra::Details::countMergeUnsortedIndices;
using Tpetra::Details::countMergeSortedIndices;
using Tpetra::Details::mergeUnsortedIndices;
using Tpetra::Details::mergeSortedIndices;

TEUCHOS_UNIT_TEST( Merge, CountUnsorted )
{
  typedef long long ordinal_type;
  typedef int index_type;

  // Repeats are allowed in the input, but NOT in the current indices.
  // This is because the whole point of these methods is to keep the
  // current indices without repeats ("merged in").

  // Thoroughly test the case of zero current indices.  This is a
  // common case, where the row has no entries currently and the user
  // wants to insert some entries.
  {
    const index_type numCurInds = 0;
    const ordinal_type* curInds = NULL;

    {
      const index_type numInputInds = 0;
      const ordinal_type* inputInds = NULL;
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {42};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {42, 42};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 3;
      const ordinal_type inputInds[] = {42, 44, 42};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 4;
      const ordinal_type inputInds[] = {42, 44, 42, 42};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 3;
      const ordinal_type inputInds[] = {6, 2, 4};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 6;
      const ordinal_type inputInds[] = {6, 6, 6, 2, 4, 4};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 7;
      const ordinal_type inputInds[] = {6, 6, 6, 2, 2, 4, 4};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
  }

  // Thoroughly test the case of zero input indices.
  {
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;

    {
      const index_type numCurInds = 0;
      const ordinal_type* curInds = NULL;
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {42};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 3;
      const ordinal_type curInds[] = {6, 2, 4};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
  }

  // Thoroughly test the case of one current index, and one input index.
  {
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {42};
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {44};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {44};
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {42};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {6};
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {6};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }
  }

  // Thoroughly test the case of more current indices than input
  // indices.  Be sure to vary the orders and test edge cases, like
  // zero or one inputs.  Test cases where there are no merges, all
  // merges, or a mix of input that can and can't be merged.
  {
    const index_type numCurInds = 6;
    const ordinal_type curInds[] = {5, 0, 4, 1, 3, 2};

    {
      const index_type numInputInds = 0;
      const ordinal_type* inputInds = NULL;
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {1};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {1, 1};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }
    {
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {7};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {7, 7};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {7, 1};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }
    {
      const index_type numInputInds = 3;
      const ordinal_type inputInds[] = {7, 1, 1};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }
    {
      const index_type numInputInds = 4;
      const ordinal_type inputInds[] = {7, 1, 7, 1};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {9, 7};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {3, 0};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }
    {
      const index_type numInputInds = 3;
      const ordinal_type inputInds[] = {3, 0, 2};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (3) );
    }
    {
      const index_type numInputInds = 4;
      const ordinal_type inputInds[] = {3, 7, 0, 2};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (3) );
    }
    {
      const index_type numInputInds = 6;
      const ordinal_type inputInds[] = {2, 5, 3, 0, 1, 4};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (6) );
    }
  }

  // Test the case of more input indices than current indices.
  {
    const index_type numInputInds = 6;
    const ordinal_type inputInds[] = {5, 0, 4, 1, 3, 2};

    {
      const index_type numCurInds = 0;
      const ordinal_type* curInds = NULL;
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {1};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }

    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {7};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }

    {
      const index_type numCurInds = 2;
      const ordinal_type curInds[] = {7, 1};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }

    {
      const index_type numCurInds = 2;
      const ordinal_type curInds[] = {9, 7};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }

    {
      const index_type numCurInds = 2;
      const ordinal_type curInds[] = {3, 0};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }

    {
      const index_type numCurInds = 3;
      const ordinal_type curInds[] = {3, 0, 2};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (3) );
    }

    {
      const index_type numCurInds = 4;
      const ordinal_type curInds[] = {3, 7, 0, 2};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (3) );
    }

    {
      const index_type numCurInds = 6;
      const ordinal_type curInds[] = {2, 5, 3, 0, 1, 4};
      const index_type count =
        countMergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (6) );
    }
  }
}


TEUCHOS_UNIT_TEST( Merge, MergeUnsorted )
{
  typedef long long ordinal_type;
  typedef int index_type;

  // Test 0 current indices and 0 input indices,
  // with 0, 1, or 3 extra spaces for insertions.
  {
    // No current indices, no extra slots.
    const index_type numCurInds = 0;
    const index_type numSpace = 0;
    ordinal_type* curInds = NULL;
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (0) );
  }
  {
    // No current indices, but one extra slot.
    const index_type numCurInds = 0;
    const index_type numSpace = 1;
    ordinal_type curInds[1];
    curInds[0] = -1; // initialize to some flag value
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (0) );
  }
  {
    // No current indices, but 3 extra slots.
    const index_type numCurInds = 0;
    const index_type numSpace = 3;
    ordinal_type curInds[4];
    curInds[0] = -1; // initialize to some flag value
    curInds[1] = -1; // initialize to some flag value
    curInds[2] = -1; // initialize to some flag value
    curInds[3] = -1; // initialize to some flag value
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (0) );
  }

  // Test 0 current indices and 1 input index,
  // with 0, 1, or 3 extra spaces for insertions.
  {
    // No current indices, no extra slots.
    const index_type numCurInds = 0;
    const index_type numSpace = 0;
    ordinal_type* curInds = NULL;
    const index_type numInputInds = 1;
    const ordinal_type inputInds[] = {42};
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, false );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }
  {
    // No current indices, 1 extra slot.
    const index_type numCurInds = 0;
    const index_type numSpace = 1;
    ordinal_type curInds[1];
    curInds[0] = -1; // flag
    const index_type numInputInds = 1;
    const ordinal_type inputInds[] = {42};
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }
  {
    // No current indices, 3 extra slots.
    const index_type numCurInds = 0;
    const index_type numSpace = 3;
    ordinal_type curInds[3];
    curInds[0] = -1; // flag
    curInds[1] = -1; // flag
    curInds[2] = -1; // flag
    const index_type numInputInds = 1;
    const ordinal_type inputInds[] = {42};
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }

  // Test 1 current index and 0 input indices,
  // with 0, 1, or 3 extra spaces for insertions.
  {
    // No extra slots.
    const index_type numCurInds = 1;
    const index_type numSpace = 1;
    ordinal_type curInds[1];
    curInds[0] = 42;
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }
  {
    // 1 extra slot.
    const index_type numCurInds = 1;
    const index_type numSpace = 2;
    ordinal_type curInds[2];
    curInds[0] = 42;
    curInds[1] = -1; // flag
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }
  {
    // 3 extra slots.
    const index_type numCurInds = 1;
    const index_type numSpace = 4;
    ordinal_type curInds[4];
    curInds[0] = 42;
    curInds[1] = -1; // flag
    curInds[2] = -1; // flag
    curInds[3] = -1; // flag
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }

  // Test several current and input indices.
  // More current indices than input indices.
  // No merges.
  {
    // No extra slots.
    const index_type numCurInds = 3;
    const index_type numSpace = 3;
    ordinal_type curInds[3];
    curInds[0] = 46;
    curInds[1] = 42;
    curInds[2] = 44;
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 47;
    inputInds[1] = 45;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, false );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (5) );
  }
  {
    // 1 extra slot.
    const index_type numCurInds = 3;
    const index_type numSpace = 4;
    ordinal_type curInds[4];
    curInds[0] = 46;
    curInds[1] = 42;
    curInds[2] = 44;
    curInds[3] = -1; // flag
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 47;
    inputInds[1] = 45;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, false );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (5) );
  }
  {
    // 3 extra slots.
    const index_type numCurInds = 3;
    const index_type numSpace = 6;
    ordinal_type curInds[6];
    curInds[0] = 46;
    curInds[1] = 42;
    curInds[2] = 44;
    curInds[3] = -1; // flag
    curInds[4] = -1; // flag
    curInds[5] = -1; // flag
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 47;
    inputInds[1] = 45;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (5) );
  }

  // Test several current and input indices.
  // More current indices than input indices.
  // Some merges possible, but not all inputs mergeable.
  {
    // No extra slots.
    const index_type numCurInds = 3;
    const index_type numSpace = 3;
    ordinal_type curInds[3];
    curInds[0] = 46;
    curInds[1] = 42;
    curInds[2] = 44;
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 42;
    inputInds[1] = 40;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, false );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (4) );
  }
  {
    // 1 extra slot.
    const index_type numCurInds = 3;
    const index_type numSpace = 4;
    ordinal_type curInds[4];
    curInds[0] = 46;
    curInds[1] = 42;
    curInds[2] = 44;
    curInds[3] = -1; // flag
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 42;
    inputInds[1] = 40;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (4) );
  }
  {
    // 3 extra slots.
    const index_type numCurInds = 3;
    const index_type numSpace = 6;
    ordinal_type curInds[6];
    curInds[0] = 46;
    curInds[1] = 42;
    curInds[2] = 44;
    curInds[3] = -1; // flag
    curInds[4] = -1; // flag
    curInds[5] = -1; // flag
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 42;
    inputInds[1] = 40;
    const std::pair<bool, index_type> result =
      mergeUnsortedIndices<ordinal_type, index_type> (curInds,
                                                      numCurInds,
                                                      numSpace,
                                                      inputInds,
                                                      numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (4) );
  }
}


TEUCHOS_UNIT_TEST( Merge, CountSorted )
{
  typedef long long ordinal_type;
  typedef int index_type;

  // Thoroughly test the case of zero current indices.  This is a
  // common case, where the row has no entries currently and the user
  // wants to insert some entries.
  {
    const index_type numCurInds = 0;
    const ordinal_type* curInds = NULL;

    {
      const index_type numInputInds = 0;
      const ordinal_type* inputInds = NULL;
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {42};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 3;
      const ordinal_type inputInds[] = {2, 4, 6};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 4;
      const ordinal_type inputInds[] = {2, 4, 6, 6};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 5;
      const ordinal_type inputInds[] = {2, 2, 4, 6, 6};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
  }

  // Thoroughly test the case of zero input indices.
  {
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;

    {
      const index_type numCurInds = 0;
      const ordinal_type* curInds = NULL;
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {42};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 3;
      const ordinal_type curInds[] = {2, 4, 6};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
  }

  // Thoroughly test the case of one current index, and one input index.
  {
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {42};
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {44};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {44};
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {42};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {6};
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {6};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }
  }

  // Thoroughly test the case of more current indices than input
  // indices.  Be sure to vary the orders and test edge cases, like
  // zero or one inputs.  Test cases where there are no merges, all
  // merges, or a mix of input that can and can't be merged.
  {
    const index_type numCurInds = 6;
    const ordinal_type curInds[] = {0, 1, 2, 3, 4, 5};

    {
      const index_type numInputInds = 0;
      const ordinal_type* inputInds = NULL;
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {1};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {1, 1};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }
    {
      const index_type numInputInds = 1;
      const ordinal_type inputInds[] = {7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {7, 7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {1, 7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }
    {
      const index_type numInputInds = 3;
      const ordinal_type inputInds[] = {1, 1, 7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }
    {
      const index_type numInputInds = 3;
      const ordinal_type inputInds[] = {1, 7, 7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }
    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {7, 9};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }

    {
      const index_type numInputInds = 2;
      const ordinal_type inputInds[] = {0, 3};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }

    {
      const index_type numInputInds = 3;
      const ordinal_type inputInds[] = {0, 2, 3};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (3) );
    }

    {
      const index_type numInputInds = 4;
      const ordinal_type inputInds[] = {0, 2, 3, 7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (3) );
    }

    {
      const index_type numInputInds = 6;
      const ordinal_type inputInds[] = {0, 1, 2, 3, 4, 5};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (6) );
    }
  }

  // Test the case of more input indices than current indices.
  {
    const index_type numInputInds = 6;
    const ordinal_type inputInds[] = {0, 1, 2, 3, 4, 5};

    {
      const index_type numCurInds = 0;
      const ordinal_type* curInds = NULL;
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }
    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {1};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }

    {
      const index_type numCurInds = 1;
      const ordinal_type curInds[] = {7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }

    {
      const index_type numCurInds = 2;
      const ordinal_type curInds[] = {1, 7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (1) );
    }

    {
      const index_type numCurInds = 2;
      const ordinal_type curInds[] = {7, 9};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    }

    {
      const index_type numCurInds = 2;
      const ordinal_type curInds[] = {0, 3};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (2) );
    }

    {
      const index_type numCurInds = 3;
      const ordinal_type curInds[] = {0, 2, 3};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (3) );
    }

    {
      const index_type numCurInds = 4;
      const ordinal_type curInds[] = {0, 2, 3, 7};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (3) );
    }

    {
      const index_type numCurInds = 6;
      const ordinal_type curInds[] = {0, 1, 2, 3, 4, 5};
      const index_type count =
        countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                             numCurInds,
                                                             inputInds,
                                                             numInputInds);
      TEST_EQUALITY_CONST( count, static_cast<index_type> (6) );
    }
  }
}


TEUCHOS_UNIT_TEST( Merge, MergeSorted )
{
  typedef long long ordinal_type;
  typedef int index_type;

  // Test 0 current indices and 0 input indices,
  // with 0, 1, or 3 extra spaces for insertions.
  {
    // No current indices, no extra slots.
    const index_type numCurInds = 0;
    const index_type numSpace = 0;
    ordinal_type* curInds = NULL;
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (0) );
  }
  {
    // No current indices, but one extra slot.
    const index_type numCurInds = 0;
    const index_type numSpace = 1;
    ordinal_type curInds[1];
    curInds[0] = -1; // initialize to some flag value
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (0) );
  }
  {
    // No current indices, but 3 extra slots.
    const index_type numCurInds = 0;
    const index_type numSpace = 3;
    ordinal_type curInds[4];
    curInds[0] = -1; // initialize to some flag value
    curInds[1] = -1; // initialize to some flag value
    curInds[2] = -1; // initialize to some flag value
    curInds[3] = -1; // initialize to some flag value
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (0) );
  }

  // Test 0 current indices and 1 input index,
  // with 0, 1, or 3 extra spaces for insertions.
  {
    // No current indices, no extra slots.
    const index_type numCurInds = 0;
    const index_type numSpace = 0;
    ordinal_type* curInds = NULL;
    const index_type numInputInds = 1;
    const ordinal_type inputInds[] = {42};
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, false );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }
  {
    // No current indices, 1 extra slot.
    const index_type numCurInds = 0;
    const index_type numSpace = 1;
    ordinal_type curInds[1];
    curInds[0] = -1; // flag
    const index_type numInputInds = 1;
    const ordinal_type inputInds[] = {42};
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }
  {
    // No current indices, 3 extra slots.
    const index_type numCurInds = 0;
    const index_type numSpace = 3;
    ordinal_type curInds[3];
    curInds[0] = -1; // flag
    curInds[1] = -1; // flag
    curInds[2] = -1; // flag
    const index_type numInputInds = 1;
    const ordinal_type inputInds[] = {42};
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }

  // Test 1 current index and 0 input indices,
  // with 0, 1, or 3 extra spaces for insertions.
  {
    // No extra slots.
    const index_type numCurInds = 1;
    const index_type numSpace = 1;
    ordinal_type curInds[1];
    curInds[0] = 42;
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }
  {
    // 1 extra slot.
    const index_type numCurInds = 1;
    const index_type numSpace = 2;
    ordinal_type curInds[2];
    curInds[0] = 42;
    curInds[1] = -1; // flag
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }
  {
    // 3 extra slots.
    const index_type numCurInds = 1;
    const index_type numSpace = 4;
    ordinal_type curInds[4];
    curInds[0] = 42;
    curInds[1] = -1; // flag
    curInds[2] = -1; // flag
    curInds[3] = -1; // flag
    const index_type numInputInds = 0;
    const ordinal_type* inputInds = NULL;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (1) );
  }

  // Test several current and input indices.
  // More current indices than input indices.
  // No merges.
  {
    // No extra slots.
    const index_type numCurInds = 3;
    const index_type numSpace = 3;
    ordinal_type curInds[3];
    curInds[0] = 42;
    curInds[1] = 44;
    curInds[2] = 46;
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 45;
    inputInds[1] = 47;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, false );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (5) );
  }
  {
    // 1 extra slot.
    const index_type numCurInds = 3;
    const index_type numSpace = 4;
    ordinal_type curInds[4];
    curInds[0] = 42;
    curInds[1] = 44;
    curInds[2] = 46;
    curInds[3] = -1; // flag
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 45;
    inputInds[1] = 47;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, false );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (5) );
  }
  {
    // 3 extra slots.
    const index_type numCurInds = 3;
    const index_type numSpace = 6;
    ordinal_type curInds[6];
    curInds[0] = 42;
    curInds[1] = 44;
    curInds[2] = 46;
    curInds[3] = -1; // flag
    curInds[4] = -1; // flag
    curInds[5] = -1; // flag
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 45;
    inputInds[1] = 47;

    const index_type count =
      countMergeSortedIndices<ordinal_type, index_type> (curInds,
                                                         numCurInds,
                                                         inputInds,
                                                         numInputInds);
    TEST_EQUALITY_CONST( count, static_cast<index_type> (0) );
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (5) );
  }

  // Test several current and input indices.
  // More current indices than input indices.
  // Some merges possible, but not all inputs mergeable.
  {
    // No extra slots.
    const index_type numCurInds = 3;
    const index_type numSpace = 3;
    ordinal_type curInds[3];
    curInds[0] = 42;
    curInds[1] = 44;
    curInds[2] = 46;
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 40;
    inputInds[1] = 42;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, false );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (4) );
  }
  {
    // 1 extra slot.
    const index_type numCurInds = 3;
    const index_type numSpace = 4;
    ordinal_type curInds[4];
    curInds[0] = 42;
    curInds[1] = 44;
    curInds[2] = 46;
    curInds[3] = -1; // flag
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 40;
    inputInds[1] = 42;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (4) );
  }
  {
    // 3 extra slots.
    const index_type numCurInds = 3;
    const index_type numSpace = 6;
    ordinal_type curInds[6];
    curInds[0] = 42;
    curInds[1] = 44;
    curInds[2] = 46;
    curInds[3] = -1; // flag
    curInds[4] = -1; // flag
    curInds[5] = -1; // flag
    const index_type numInputInds = 2;
    ordinal_type inputInds[2];
    inputInds[0] = 40;
    inputInds[1] = 42;
    const std::pair<bool, index_type> result =
      mergeSortedIndices<ordinal_type, index_type> (curInds,
                                                    numCurInds,
                                                    numSpace,
                                                    inputInds,
                                                    numInputInds);
    TEST_EQUALITY_CONST( result.first, true );
    TEST_EQUALITY_CONST( result.second, static_cast<index_type> (4) );
  }
}


} // namespace (anonymous)

