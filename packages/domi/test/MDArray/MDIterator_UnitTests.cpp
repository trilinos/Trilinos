// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "MDArray_UnitTest_helpers.hpp"

namespace
{

using Teuchos::tuple;
using Domi::MDArray;
using Domi::MDArrayView;
using Domi::MDArrayRCP;
using Domi::Slice;
using Teuchos::RangeError;
using MDArrayUnitTestHelpers::nrows;
using MDArrayUnitTestHelpers::ncols;
using MDArrayUnitTestHelpers::generateMDArray;

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayIteratorCtor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  typename MDArray< T >::iterator it(a);
  TEST_EQUALITY_CONST(it.index(0), 0);
  TEST_EQUALITY_CONST(it.index(1), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayViewIteratorCtor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  typename MDArrayView< T >::iterator it(av);
  TEST_EQUALITY_CONST(it.index(0), 0);
  TEST_EQUALITY_CONST(it.index(1), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayRCPIteratorCtor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  typename MDArrayRCP< T >::iterator it(a);
  TEST_EQUALITY_CONST(it.index(0), 0);
  TEST_EQUALITY_CONST(it.index(1), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayIteratorCtorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  typename MDArray< T >::iterator it(a,true);
  TEST_EQUALITY_CONST(it.index(0), 2);
  TEST_EQUALITY_CONST(it.index(1), 5);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayViewIteratorCtorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  typename MDArrayView< T >::iterator it(av,true);
  TEST_EQUALITY_CONST(it.index(0), 3);
  TEST_EQUALITY_CONST(it.index(1), 4);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayRCPIteratorCtorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  typename MDArrayRCP< T >::iterator it(a,true);
  TEST_EQUALITY_CONST(it.index(0), 4);
  TEST_EQUALITY_CONST(it.index(1), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayIteratorBegin, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  typename MDArray< T >::iterator it = a.begin();
  TEST_EQUALITY_CONST(it.index(0), 0);
  TEST_EQUALITY_CONST(it.index(1), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayViewIteratorBegin, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  typename MDArrayView< T >::iterator it = av.begin();
  TEST_EQUALITY_CONST(it.index(0), 0);
  TEST_EQUALITY_CONST(it.index(1), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayRCPIteratorBegin, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  typename MDArrayRCP< T >::iterator it = a.begin();
  TEST_EQUALITY_CONST(it.index(0), 0);
  TEST_EQUALITY_CONST(it.index(1), 0);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayIteratorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  typename MDArray< T >::iterator it = a.end();
  TEST_EQUALITY_CONST(it.index(0), 2);
  TEST_EQUALITY_CONST(it.index(1), 5);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayViewIteratorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  typename MDArrayView< T >::iterator it = av.end();
  TEST_EQUALITY_CONST(it.index(0), 3);
  TEST_EQUALITY_CONST(it.index(1), 4);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayRCPIteratorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  typename MDArrayRCP< T >::iterator it = a.end();
  TEST_EQUALITY_CONST(it.index(0), 4);
  TEST_EQUALITY_CONST(it.index(1), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayIteratorIndex, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  Teuchos::ArrayView< dim_type > index = tuple< dim_type >(1,4);
  typename MDArray< T >::iterator it(a, index);
  TEST_EQUALITY_CONST(it.index(0), 1);
  TEST_EQUALITY_CONST(it.index(1), 4);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayViewIteratorIndex, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  Teuchos::ArrayView< dim_type > index = tuple< dim_type >(2,3);
  typename MDArrayView< T >::iterator it(av, index);
  TEST_EQUALITY_CONST(it.index(0), 2);
  TEST_EQUALITY_CONST(it.index(1), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayRCPIteratorIndex, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  Teuchos::ArrayView< dim_type > index = tuple< dim_type >(3,2);
  typename MDArrayRCP< T >::iterator it(a, index);
  TEST_EQUALITY_CONST(it.index(0), 3);
  TEST_EQUALITY_CONST(it.index(1), 2);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, MDArrayViewAssign, T )
{
  typedef typename Domi::dim_type ord;
  MDArray< T > a(tuple< ord >(8,5));
  // Initialize all of a to the value of 5.  This uses the
  // MDArray<T>::iterator
  a.assign(5);
  // Take a view of the left half of the array
  MDArrayView< T > av = a[Slice(4)][Slice()];
  // Set all of av (the left hand side of a) to the value of -5.  This
  // uses the MDArrayView<T>::iterator
  av.assign(-5);
  // Check that the array a has the correct values.
  for (ord j = 0; j < 5; ++j)
    for (ord i = 0; i < 8; ++i)
    {
      if (i < 4)
      {
        TEST_EQUALITY_CONST(a(i,j), -5);
      }
      else
      {
        TEST_EQUALITY_CONST(a(i,j),  5);
      }
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, FirstIndexFastest, T )
{
  typedef typename Domi::dim_type ord;
  MDArrayRCP< T > a(tuple< ord >(4,7), Domi::FIRST_INDEX_FASTEST);
  typename MDArrayRCP< T >::iterator it = a.begin();
  for (ord j = 0; j < 7; ++j)
    for (ord i = 0; i < 4; ++i)
    {
      TEST_EQUALITY(it.index(0), i);
      TEST_EQUALITY(it.index(1), j);
      ++it;
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDIterator, LastIndexFastest, T )
{
  typedef typename Domi::dim_type ord;
  MDArrayRCP< T > a(tuple< ord >(5,6), Domi::LAST_INDEX_FASTEST);
  typename MDArrayRCP< T >::iterator it = a.begin();
  for (ord i = 0; i < 5; ++i)
    for (ord j = 0; j < 6; ++j)
    {
      TEST_EQUALITY(it.index(0), i);
      TEST_EQUALITY(it.index(1), j);
      ++it;
    }
}

//
// Instantiations
//
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
#  define DEBUG_UNIT_TEST_GROUP( T )
#else // HAVE_DOMI_ARRAY_BOUNDSCHECK
#  define DEBUG_UNIT_TEST_GROUP( T )
#endif // HAVE_DOMI_ARRAY_BOUNDSCHECK

#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayIteratorCtor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayViewIteratorCtor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayRCPIteratorCtor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayIteratorCtorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayViewIteratorCtorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayRCPIteratorCtorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayIteratorBegin, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayViewIteratorBegin, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayRCPIteratorBegin, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayIteratorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayViewIteratorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayRCPIteratorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayIteratorIndex, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayViewIteratorIndex, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayRCPIteratorIndex, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, MDArrayViewAssign, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, FirstIndexFastest, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDIterator, LastIndexFastest, T ) \
  DEBUG_UNIT_TEST_GROUP( T )

UNIT_TEST_GROUP(int)
#if 1
UNIT_TEST_GROUP(long)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)
#endif

}  // namespace
