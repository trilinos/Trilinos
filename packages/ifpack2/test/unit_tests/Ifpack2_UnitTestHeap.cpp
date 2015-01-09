/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

// ***********************************************************************
// 
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************


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

