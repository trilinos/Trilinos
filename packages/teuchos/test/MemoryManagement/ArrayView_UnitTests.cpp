/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER
*/

#include "Array_UnitTest_helpers.hpp"

#include "Teuchos_implicit_cast.hpp"


namespace {


using ArrayUnitTestHelpers::n;
using ArrayUnitTestHelpers::generateArray;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::arcp;
using Teuchos::ArrayView;
using Teuchos::arrayView;
using Teuchos::av_const_cast;
using Teuchos::av_reinterpret_cast;
using Teuchos::DanglingReferenceError;
using Teuchos::as;
using Teuchos::null;
using Teuchos::implicit_ptr_cast;


TEUCHOS_UNIT_TEST( ArrayView, assignSelf )
{
  ArrayView<int> av;
  av = av;
  TEST_ASSERT(is_null(av));
  TEST_ASSERT(!nonnull(av));
}


TEUCHOS_UNIT_TEST( ArrayView, assignFuncSelf )
{
  Array<int> a = generateArray<int>(n);
  ArrayView<int> av = a;
  av.assign(av);
}


TEUCHOS_UNIT_TEST( ArrayView, av_const_cast_null )
{
  ArrayView<const int> av_int1 = null;
  ArrayView<int> av_int2 = av_const_cast<int>(av_int1);
  TEST_ASSERT(is_null(av_int2));
}


TEUCHOS_UNIT_TEST( ArrayView, av_const_cast )
{
  ArrayRCP<const int> arcp_int = arcp<int>(n);
  ArrayView<const int> av_int1 = arcp_int();
  ArrayView<int> av_int2 = av_const_cast<int>(av_int1);
  TEST_ASSERT(nonnull(av_int2));
  TEST_EQUALITY(av_int2.getRawPtr(), av_int1.getRawPtr());
  TEST_EQUALITY(av_int2.getRawPtr(), arcp_int.getRawPtr());
}


TEUCHOS_UNIT_TEST( ArrayView, av_reinterpret_cast_null )
{
  ArrayView<char> av_char = null;
  ArrayView<int> av_int = av_reinterpret_cast<int>(av_char);
  TEST_ASSERT(is_null(av_int));
}


TEUCHOS_UNIT_TEST( ArrayView, av_reinterpret_cast_char_to_int )
{

  const int sizeOfInt = sizeof(int);
  const int sizeOfChar = sizeof(char);
  const int num_ints = n;
  const int num_chars = (num_ints*sizeOfInt)/sizeOfChar;
  out << "num_ints = " << num_ints << "\n";
  out << "num_chars = " << num_chars << "\n";

  ArrayRCP<char> arcp_char = arcp<char>(num_chars);
  ArrayView<int> av_int = av_reinterpret_cast<int>(arcp_char());
  TEST_EQUALITY(av_int.size(), num_ints);
  TEST_EQUALITY(implicit_ptr_cast<void>(&av_int[0]),
    implicit_ptr_cast<void>(&arcp_char[0]));
  TEST_EQUALITY(implicit_ptr_cast<void>((&av_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&arcp_char[num_chars-1])+1));

}


TEUCHOS_UNIT_TEST( ArrayView, av_reinterpret_cast_int_to_char )
{

  const int sizeOfInt = sizeof(int);
  const int sizeOfChar = sizeof(char);
  const int num_ints = n;
  const int num_chars = (num_ints*sizeOfInt)/sizeOfChar;
  out << "num_ints = " << num_ints << "\n";
  out << "num_chars = " << num_chars << "\n";

  ArrayRCP<int> arcp_int = arcp<int>(num_ints);
  ArrayView<char> av_char = av_reinterpret_cast<char>(arcp_int());
  TEST_EQUALITY(av_char.size(), num_chars);
  TEST_EQUALITY(implicit_ptr_cast<void>(&arcp_int[0]),
    implicit_ptr_cast<void>(&av_char[0]));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&av_char[num_chars-1])+1));
  TEST_EQUALITY(implicit_ptr_cast<void>((&arcp_int[num_ints-1])+1),
    implicit_ptr_cast<void>((&av_char[num_chars-1])+1));

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, arrayView_construct_zero_size, T )
{
  Array<T> a;
  const ArrayView<T> av = arrayView(a.getRawPtr(), a.size());
  TEST_EQUALITY_CONST(av.size(), 0);
  TEST_ASSERT(is_null(av));
  TEST_ASSERT(!nonnull(av));
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, arrayView, T )
{
  Array<T> a = generateArray<T>(n);
  const ArrayView<T> av = arrayView(&a[0], a.size());
  TEST_COMPARE_ARRAYS( a, av );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, view_size_zero, T )
{
  Array<T> a = generateArray<T>(n);
  const ArrayView<T> av = a();
  TEST_THROW(const ArrayView<T> av_zero1 = av(0,0), Teuchos::RangeError);
  TEST_THROW(const ArrayView<T> av_zero2 = av(1,0), Teuchos::RangeError);
  TEST_THROW(const ArrayView<T> av_zero2 = av(n-1,0), Teuchos::RangeError);
  try {
    const ArrayView<T> av_zero1 = av(0,0);
    TEST_ASSERT(0); // If you get here you failed!
  }
  catch(const Teuchos::RangeError& except) {
    const std::string exceptMsg = except.what(); 
    out << "Exception message: " << exceptMsg << "\n";
    const size_t found_substr_i = exceptMsg.find("Error, size=0 is not allowed");
    TEST_INEQUALITY(found_substr_i, std::string::npos);
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, assignmentOperator, T )
{
  Array<T> a = generateArray<T>(n);
  ArrayView<T> av1;
  av1 = a;
  ArrayView<T> av2;
  av2 = av1;
  TEST_EQUALITY( av1.getRawPtr(), a.getRawPtr() );
  TEST_EQUALITY( av1.size(), as<int>(a.size()) );
  TEST_EQUALITY( av1.getRawPtr(), av2.getRawPtr() );
  TEST_EQUALITY( av1.size(), av2.size() );
  TEST_COMPARE_ARRAYS( av1, a );
  TEST_COMPARE_ARRAYS( av1, av2 );
  av1 = null;
  TEST_EQUALITY_CONST( av1.getRawPtr(), 0 );
  TEST_EQUALITY_CONST( av1.size(), 0 );
  av2 = null;
  TEST_EQUALITY_CONST( av2.getRawPtr(), 0 );
  TEST_EQUALITY_CONST( av2.size(), 0 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, iterators, T )
{
  typedef typename ArrayView<T>::iterator iter_t;
  typedef Teuchos::ScalarTraits<T> ST;
  ECHO(Array<T> a = generateArray<T>(n));
  ECHO(ArrayView<T> av = a);
  ECHO(const iter_t av_begin = av.begin());
  ECHO(const iter_t av_end = av.end());
#ifdef TEUCHOS_DEBUG
  TEST_ASSERT(av_begin.shares_resource(av_end));
#endif
  ECHO(std::fill(av_begin, av_end, ST::random()));
  ECHO(Array<T> a2 = generateArray<T>(n));
  ECHO(ArrayView<T> av2 = a2);
  ECHO(std::copy(av.begin(), av.end(), av2.begin()));
  TEST_COMPARE_ARRAYS(a, a2);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, danglingView_std_vector, T )
{
  ArrayView<T> av;
  T* badPtr = 0;
  {
    std::vector<T> v(n);
    av = v;
    badPtr = &v[0];
  }
  // Access the raw pointer but it now points to invalid memory!
  TEST_EQUALITY(av.getRawPtr(), badPtr);
  // Above, we have no way to detect that the underlying std::vector object
  // has gone away.  This is the whole point of needing Teuchos::Array and
  // having an integrated set of utility classes that all work together!
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayView, danglingView_rcp_std_vector, T )
{
  ArrayView<T> av;
  {
    ArrayRCP<T> ap = arcp(rcp(new std::vector<T>(n)));
    av = ap();
  }
#ifdef TEUCHOS_DEBUG
  TEST_THROW(av.getRawPtr(), DanglingReferenceError);
#endif
  // Above, because we wrapped the initial std::vector in an RCP object, we
  // can sucessfully detect when the object goes away in debug mode!
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


//
// Instantiations
//



#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, view_size_zero, T )

#else // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, arrayView_construct_zero_size, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, arrayView, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, assignmentOperator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, iterators, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, danglingView_std_vector, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayView, danglingView_rcp_std_vector, T ) \
  DEBUG_UNIT_TEST_GROUP( T )


UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)


} // namespace
