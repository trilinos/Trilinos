// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/*! \file Ifpack2_UnitTestHeap.cpp

\brief Ifpack2 Unit test for the Heap templates.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Ifpack2_Version.hpp>
#include <iostream>

#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_Heap.hpp>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2Heap, Test1, Scalar, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  Teuchos::Array<GlobalOrdinal> heap1;
  typename Teuchos::Array<GlobalOrdinal>::size_type heap1_len = 0;

  GlobalOrdinal idx = 1;

  Ifpack2::add_to_heap(idx, heap1, heap1_len);
  ++idx;
  Ifpack2::add_to_heap(idx, heap1, heap1_len);
  ++idx;
  Ifpack2::add_to_heap(idx, heap1, heap1_len);

  TEST_EQUALITY(heap1.size() == 3, true);
  TEST_EQUALITY(heap1_len == 3, true);
  TEST_EQUALITY(heap1.front() == 1, true);

  Ifpack2::rm_heap_root(heap1, heap1_len);
  TEST_EQUALITY(heap1.size() == 3, true);
  TEST_EQUALITY(heap1_len == 2, true);
  TEST_EQUALITY(heap1.front() == 2, true);

  Ifpack2::rm_heap_root(heap1, heap1_len);
  TEST_EQUALITY(heap1.size() == 3, true);
  TEST_EQUALITY(heap1_len == 1, true);
  TEST_EQUALITY(heap1.front() == 3, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(Ifpack2Heap, Test3, Scalar, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Ifpack2::Version();
  out << "Ifpack2::Version(): " << version << std::endl;

  Teuchos::Array<Scalar> vals(4);
  vals[0] = 3;
  vals[1] = 2;
  vals[2] = 1;
  vals[3] = 0;

  Ifpack2::greater_indirect<Scalar,GlobalOrdinal> vals_comp(vals);

  Teuchos::Array<GlobalOrdinal> heap1;
  typename Teuchos::Array<GlobalOrdinal>::size_type heap1_len = 0;

  GlobalOrdinal idx = 0;

  Ifpack2::add_to_heap(idx, heap1, heap1_len, vals_comp);
  ++idx;
  Ifpack2::add_to_heap(idx, heap1, heap1_len, vals_comp);
  ++idx;
  Ifpack2::add_to_heap(idx, heap1, heap1_len, vals_comp);
  ++idx;
  Ifpack2::add_to_heap(idx, heap1, heap1_len, vals_comp);

  TEST_EQUALITY(heap1.size(), 4);
  TEST_EQUALITY(heap1_len, 4);
  TEST_EQUALITY(heap1.front(), 3);

  Ifpack2::rm_heap_root(heap1, heap1_len, vals_comp);
  TEST_EQUALITY(heap1.size(), 4);
  TEST_EQUALITY(heap1_len, 3);
  TEST_EQUALITY(heap1.front(), 2);

  Ifpack2::rm_heap_root(heap1, heap1_len, vals_comp);
  TEST_EQUALITY(heap1.size(), 4);
  TEST_EQUALITY(heap1_len, 2);
  TEST_EQUALITY(heap1.front(), 1);

  Ifpack2::rm_heap_root(heap1, heap1_len, vals_comp);
  TEST_EQUALITY(heap1_len, 1);
  TEST_EQUALITY(heap1.front(), 0);

  Ifpack2::rm_heap_root(heap1, heap1_len, vals_comp);
  TEST_EQUALITY(heap1_len, 0);
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2Heap, Test1, Scalar, GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Ifpack2Heap, Test3, Scalar, GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int)
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, int)

}//namespace <anonymous>

