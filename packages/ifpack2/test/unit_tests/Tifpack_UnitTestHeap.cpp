// ***********************************************************************
// 
//      Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


/*! \file Tifpack_UnitTestHeap.cpp

\brief Tifpack Unit test for the Heap templates.
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Tifpack_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Tifpack_Version.hpp>
#include <iostream>

#include <Tifpack_UnitTestHelpers.hpp>
#include <Tifpack_Heap.hpp>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(TifpackHeap, Test1, Scalar, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Tifpack::Version();
  out << "Tifpack::Version(): " << version << std::endl;

  Teuchos::Array<GlobalOrdinal> heap1;
  typename Teuchos::Array<GlobalOrdinal>::size_type heap1_len = 0;

  GlobalOrdinal idx = 1;

  Tifpack::add_to_heap(idx, heap1, heap1_len);
  ++idx;
  Tifpack::add_to_heap(idx, heap1, heap1_len);
  ++idx;
  Tifpack::add_to_heap(idx, heap1, heap1_len);

  TEUCHOS_TEST_EQUALITY(heap1.size() == 3, true, out, success);
  TEUCHOS_TEST_EQUALITY(heap1_len == 3, true, out, success);
  TEUCHOS_TEST_EQUALITY(heap1.front() == 1, true, out, success);

  Tifpack::rm_heap_root(heap1, heap1_len);
  TEUCHOS_TEST_EQUALITY(heap1.size() == 3, true, out, success);
  TEUCHOS_TEST_EQUALITY(heap1_len == 2, true, out, success);
  TEUCHOS_TEST_EQUALITY(heap1.front() == 2, true, out, success);

  Tifpack::rm_heap_root(heap1, heap1_len);
  TEUCHOS_TEST_EQUALITY(heap1.size() == 3, true, out, success);
  TEUCHOS_TEST_EQUALITY(heap1_len == 1, true, out, success);
  TEUCHOS_TEST_EQUALITY(heap1.front() == 3, true, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(TifpackHeap, Test3, Scalar, GlobalOrdinal)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  std::string version = Tifpack::Version();
  out << "Tifpack::Version(): " << version << std::endl;

  Teuchos::Array<Scalar> vals(4);
  vals[0] = 3;
  vals[1] = 2;
  vals[2] = 1;
  vals[3] = 0;

  Tifpack::greater_indirect<Scalar,GlobalOrdinal> vals_comp(vals);

  Teuchos::Array<GlobalOrdinal> heap1;
  typename Teuchos::Array<GlobalOrdinal>::size_type heap1_len = 0;

  GlobalOrdinal idx = 0;

  Tifpack::add_to_heap(idx, heap1, heap1_len, vals_comp);
  ++idx;
  Tifpack::add_to_heap(idx, heap1, heap1_len, vals_comp);
  ++idx;
  Tifpack::add_to_heap(idx, heap1, heap1_len, vals_comp);
  ++idx;
  Tifpack::add_to_heap(idx, heap1, heap1_len, vals_comp);

  TEUCHOS_TEST_EQUALITY(heap1.size(), 4, out, success);
  TEUCHOS_TEST_EQUALITY(heap1_len, 4, out, success);
  TEUCHOS_TEST_EQUALITY(heap1.front(), 3, out, success);

  Tifpack::rm_heap_root(heap1, heap1_len, vals_comp);
  TEUCHOS_TEST_EQUALITY(heap1.size(), 4, out, success);
  TEUCHOS_TEST_EQUALITY(heap1_len, 3, out, success);
  TEUCHOS_TEST_EQUALITY(heap1.front(), 2, out, success);

  Tifpack::rm_heap_root(heap1, heap1_len, vals_comp);
  TEUCHOS_TEST_EQUALITY(heap1.size(), 4, out, success);
  TEUCHOS_TEST_EQUALITY(heap1_len, 2, out, success);
  TEUCHOS_TEST_EQUALITY(heap1.front(), 1, out, success);

  Tifpack::rm_heap_root(heap1, heap1_len, vals_comp);
  TEUCHOS_TEST_EQUALITY(heap1_len, 1, out, success);
  TEUCHOS_TEST_EQUALITY(heap1.front(), 0, out, success);

  Tifpack::rm_heap_root(heap1, heap1_len, vals_comp);
  TEUCHOS_TEST_EQUALITY(heap1_len, 0, out, success);
}

#define UNIT_TEST_GROUP_SCALAR_ORDINAL(Scalar,GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( TifpackHeap, Test1, Scalar, GlobalOrdinal) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( TifpackHeap, Test3, Scalar, GlobalOrdinal)

UNIT_TEST_GROUP_SCALAR_ORDINAL(double, int)
UNIT_TEST_GROUP_SCALAR_ORDINAL(float, int)

}//namespace <anonymous>

