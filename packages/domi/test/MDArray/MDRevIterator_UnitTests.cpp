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

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRevIteratorCtor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  typename MDArray< T >::reverse_iterator it(a);
  TEST_EQUALITY_CONST(it.index(0), 1);
  TEST_EQUALITY_CONST(it.index(1), 4);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayViewRevIteratorCtor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  typename MDArrayView< T >::reverse_iterator it(av);
  TEST_EQUALITY_CONST(it.index(0), 2);
  TEST_EQUALITY_CONST(it.index(1), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRCPRevIteratorCtor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  typename MDArrayRCP< T >::reverse_iterator it(a);
  TEST_EQUALITY_CONST(it.index(0), 3);
  TEST_EQUALITY_CONST(it.index(1), 2);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRevIteratorCtorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  typename MDArray< T >::reverse_iterator it(a,true);
  TEST_EQUALITY_CONST(it.index(0), -1);
  TEST_EQUALITY_CONST(it.index(1), -1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayViewRevIteratorCtorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  typename MDArrayView< T >::reverse_iterator it(av,true);
  TEST_EQUALITY_CONST(it.index(0), -1);
  TEST_EQUALITY_CONST(it.index(1), -1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRCPRevIteratorCtorEnd, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  typename MDArrayRCP< T >::reverse_iterator it(a,true);
  TEST_EQUALITY_CONST(it.index(0), -1);
  TEST_EQUALITY_CONST(it.index(1), -1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRevIteratorRbegin, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  typename MDArray< T >::reverse_iterator it = a.rbegin();
  TEST_EQUALITY_CONST(it.index(0), 1);
  TEST_EQUALITY_CONST(it.index(1), 4);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayViewRevIteratorRbegin, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  typename MDArrayView< T >::reverse_iterator it = av.rbegin();
  TEST_EQUALITY_CONST(it.index(0), 2);
  TEST_EQUALITY_CONST(it.index(1), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRCPRevIteratorRbegin, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  typename MDArrayRCP< T >::reverse_iterator it = a.rbegin();
  TEST_EQUALITY_CONST(it.index(0), 3);
  TEST_EQUALITY_CONST(it.index(1), 2);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRevIteratorRend, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  typename MDArray< T >::reverse_iterator it = a.rend();
  TEST_EQUALITY_CONST(it.index(0), -1);
  TEST_EQUALITY_CONST(it.index(1), -1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayViewRevIteratorRend, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  typename MDArrayView< T >::reverse_iterator it = av.rend();
  TEST_EQUALITY_CONST(it.index(0), -1);
  TEST_EQUALITY_CONST(it.index(1), -1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRCPRevIteratorRend, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  typename MDArrayRCP< T >::reverse_iterator it = a.rend();
  TEST_EQUALITY_CONST(it.index(0), -1);
  TEST_EQUALITY_CONST(it.index(1), -1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRevIteratorIndex, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,5));
  Teuchos::ArrayView< dim_type > index = tuple< dim_type >(1,4);
  typename MDArray< T >::reverse_iterator it(a, index);
  TEST_EQUALITY_CONST(it.index(0), 1);
  TEST_EQUALITY_CONST(it.index(1), 4);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayViewRevIteratorIndex, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,4));
  MDArrayView< T > av = a();
  Teuchos::ArrayView< dim_type > index = tuple< dim_type >(2,3);
  typename MDArrayView< T >::reverse_iterator it(av, index);
  TEST_EQUALITY_CONST(it.index(0), 2);
  TEST_EQUALITY_CONST(it.index(1), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, MDArrayRCPRevIteratorIndex, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a(tuple< dim_type >(4,3));
  Teuchos::ArrayView< dim_type > index = tuple< dim_type >(3,2);
  typename MDArrayRCP< T >::reverse_iterator it(a, index);
  TEST_EQUALITY_CONST(it.index(0), 3);
  TEST_EQUALITY_CONST(it.index(1), 2);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, FirstIndexFastest, T )
{
  typedef typename Domi::dim_type ord;
  MDArrayRCP< T > a(tuple< ord >(4,7), Domi::FIRST_INDEX_FASTEST);
  typename MDArrayRCP< T >::reverse_iterator it = a.rbegin();
  for (ord j = 6; j >= 0; --j)
    for (ord i = 3; i >= 0; --i)
    {
      TEST_EQUALITY(it.index(0), i);
      TEST_EQUALITY(it.index(1), j);
      ++it;
    }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDRevIterator, LastIndexFastest, T )
{
  typedef typename Domi::dim_type ord;
  MDArrayRCP< T > a(tuple< ord >(5,6), Domi::LAST_INDEX_FASTEST);
  typename MDArrayRCP< T >::reverse_iterator it = a.rbegin();
  for (ord i = 4; i >= 0; --i)
    for (ord j = 5; j >= 0; --j)
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
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRevIteratorCtor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayViewRevIteratorCtor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRCPRevIteratorCtor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRevIteratorCtorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayViewRevIteratorCtorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRCPRevIteratorCtorEnd, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRevIteratorRbegin, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayViewRevIteratorRbegin, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRCPRevIteratorRbegin, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRevIteratorRend, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayViewRevIteratorRend, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRCPRevIteratorRend, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRevIteratorIndex, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayViewRevIteratorIndex, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, MDArrayRCPRevIteratorIndex, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, FirstIndexFastest, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDRevIterator, LastIndexFastest, T ) \
  DEBUG_UNIT_TEST_GROUP( T )

UNIT_TEST_GROUP(int)
#if 1
UNIT_TEST_GROUP(long)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)
#endif

}  // namespace
