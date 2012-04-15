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
#include <fei_CommMap.hpp>
#include <iostream>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(CommMap, test0, T)
{
  typename fei::CommMap<T>::Type comm_map;

  std::vector<T> input_items(4);
  input_items[0] = 2;
  input_items[1] = 3;
  input_items[2] = 1;
  input_items[3] = 0;

  fei::addItemsToCommMap(0, input_items.size(), &input_items[0], comm_map);
  fei::addItemsToCommMap(1, input_items.size(), &input_items[0], comm_map, false);

  std::vector<T>& sorted_items = comm_map[0];
  std::vector<T>& unsorted_items = comm_map[1];

  std::vector<T> expected_sorted_items(4);
  expected_sorted_items[0] = 0;
  expected_sorted_items[1] = 1;
  expected_sorted_items[2] = 2;
  expected_sorted_items[3] = 3;

  bool sorted_correct = sorted_items == expected_sorted_items;
  bool unsorted_correct = unsorted_items == input_items;

  TEUCHOS_TEST_EQUALITY(sorted_correct, true, out, success);
  TEUCHOS_TEST_EQUALITY(unsorted_correct, true, out, success);
}

#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(CommMap,test0,TYPE)

typedef long int longint;
UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(longint)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)

}//namespace <anonymous>

