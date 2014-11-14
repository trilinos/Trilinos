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
using Domi::Slice;
using Domi::RangeError;
using MDArrayUnitTestHelpers::nrows;
using MDArrayUnitTestHelpers::ncols;
using MDArrayUnitTestHelpers::generateMDArray;

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, defaultConstructor, T )
{
  MDArray< T > a;
  TEST_EQUALITY_CONST(a.numDims()   , 1);
  TEST_EQUALITY_CONST(a.dimension(0), 0);
  TEST_EQUALITY_CONST(a.size()      , 0);
  TEST_EQUALITY_CONST(a.strides()[0], 1);
  TEUCHOS_ASSERT(a.empty());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, simpleConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(2,3));
  TEST_EQUALITY_CONST(a.numDims()   , 2);
  TEST_EQUALITY_CONST(a.dimension(0), 2);
  TEST_EQUALITY_CONST(a.dimension(1), 3);
  TEST_EQUALITY_CONST(a.size()      , 6);
  TEST_EQUALITY_CONST(a.strides()[0], 1);
  TEST_EQUALITY_CONST(a.strides()[1], 2);
  TEST_EQUALITY_CONST(a.layout(), Domi::DEFAULT_ORDER);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, dimsAndValConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  T val = 1;
  MDArray< T > a(tuple< dim_type >(2,2), val);
  TEST_EQUALITY_CONST(a.numDims()   , 2);
  TEST_EQUALITY_CONST(a.dimension(0), 2);
  TEST_EQUALITY_CONST(a.dimension(1), 2);
  TEST_EQUALITY_CONST(a.size()      , 4);
  TEST_EQUALITY_CONST(a.strides()[0], 1);
  TEST_EQUALITY_CONST(a.strides()[1], 2);
  TEST_EQUALITY_CONST(a.layout(), Domi::DEFAULT_ORDER);
  TEST_EQUALITY_CONST(a(0,0), val);
  TEST_EQUALITY_CONST(a(0,1), val);
  TEST_EQUALITY_CONST(a(1,0), val);
  TEST_EQUALITY_CONST(a(1,1), val);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, dimsAndOrderConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a(tuple< dim_type >(3,3), Domi::C_ORDER);
  TEST_EQUALITY_CONST(a.numDims()   , 2);
  TEST_EQUALITY_CONST(a.dimension(0), 3);
  TEST_EQUALITY_CONST(a.dimension(1), 3);
  TEST_EQUALITY_CONST(a.size()      , 9);
  TEST_EQUALITY_CONST(a.strides()[0], 3);
  TEST_EQUALITY_CONST(a.strides()[1], 1);
  TEST_EQUALITY_CONST(a.layout(), Domi::C_ORDER);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, dimsValAndOrderConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  T val = 2;
  MDArray< T > a(tuple< dim_type >(2,2,2), val, Domi::ROW_MAJOR);
  TEST_EQUALITY_CONST(a.numDims()   , 3);
  TEST_EQUALITY_CONST(a.dimension(0), 2);
  TEST_EQUALITY_CONST(a.dimension(1), 2);
  TEST_EQUALITY_CONST(a.dimension(2), 2);
  TEST_EQUALITY_CONST(a.size()      , 8);
  TEST_EQUALITY_CONST(a.strides()[0], 4);
  TEST_EQUALITY_CONST(a.strides()[1], 2);
  TEST_EQUALITY_CONST(a.strides()[2], 1);
  TEST_EQUALITY_CONST(a.layout(), Domi::ROW_MAJOR);
  TEST_EQUALITY_CONST(a(0,0,0), val);
  TEST_EQUALITY_CONST(a(0,0,1), val);
  TEST_EQUALITY_CONST(a(0,1,0), val);
  TEST_EQUALITY_CONST(a(0,1,1), val);
  TEST_EQUALITY_CONST(a(1,0,0), val);
  TEST_EQUALITY_CONST(a(1,0,1), val);
  TEST_EQUALITY_CONST(a(1,1,0), val);
  TEST_EQUALITY_CONST(a(1,1,1), val);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, copyConstructor, T )
{
  MDArray< T > a = generateMDArray< T >(3,4);
  MDArray< T > b(a);
  TEUCHOS_ASSERT( a == b );
  TEUCHOS_ASSERT( b == a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, mdArrayViewConstructor0, T )
{
  MDArray< T >     a = generateMDArray< T >(4,2);
  MDArrayView< T > b = a.mdArrayView();
  MDArray< T >     c(b);
  TEUCHOS_ASSERT( a == c );
  TEUCHOS_ASSERT( c == a );
  TEUCHOS_ASSERT( b == c );
  TEUCHOS_ASSERT( c == b );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, mdArrayViewConstructor1, T )
{
  MDArray< T >     a = generateMDArray< T >(5,6);
  MDArrayView< T > b = a[Slice(1,4)][Slice(2,5)];
  MDArray< T >     c(b);
  TEUCHOS_ASSERT( b == c );
  TEUCHOS_ASSERT( c == b );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, equalOperatorDiffLayout, T )
{
  MDArray< T > a = generateMDArray< T >(3,4,Domi::FORTRAN_ORDER);
  MDArray< T > b = generateMDArray< T >(3,4,Domi::C_ORDER);
  TEST_EQUALITY((a == b), false);
  TEUCHOS_ASSERT( a != b );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, equalOperatorMDArrayView, T )
{
  MDArray< T > a = generateMDArray< T >(3,4);
  MDArrayView< T > b(a);
  TEUCHOS_ASSERT( a == b );
  TEUCHOS_ASSERT( b == a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, inequalityOperator, T )
{
  //typedef typename Domi::dim_type dim_type;
  MDArray< T > a = generateMDArray< T >(4,4);
  MDArray< T > b(a);
  b(1,2) = b(1,2) * b(1,2);
  TEUCHOS_ASSERT( b != a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, indexing4D, T )
{
  typedef Domi::Ordinal ord;
  using Teuchos::tuple;
  ord ni = 2;
  ord nj = 2;
  ord nk = 3;
  ord nm = 3;
  MDArray< T > mda(tuple(ni,nj,nk,nm));
  for (ord m=0; m < nm; m++)
    for (ord k=0; k < nk; k++)
      for (ord j=0; j < nj; j++)
        for (ord i=0; i < ni; i++)
        {
          mda(i,j,k,m) = 3;
          TEST_EQUALITY(mda(i,j,k,m), 3);
        }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, indexing5D, T )
{
  typedef Domi::Ordinal ord;
  using Teuchos::tuple;
  ord ni = 2;
  ord nj = 2;
  ord nk = 3;
  ord nm = 3;
  ord nn = 4;
  MDArray< T > mda(tuple(ni,nj,nk,nm,nn));
  for (ord n=0; n < nn; n++)
    for (ord m=0; m < nm; m++)
      for (ord k=0; k < nk; k++)
        for (ord j=0; j < nj; j++)
          for (ord i=0; i < ni; i++)
          {
            mda(i,j,k,m,n) = 4;
            TEST_EQUALITY(mda(i,j,k,m,n), 4);
          }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, indexing6D, T )
{
  typedef Domi::Ordinal ord;
  using Teuchos::tuple;
  ord ni = 2;
  ord nj = 2;
  ord nk = 3;
  ord nm = 3;
  ord nn = 4;
  ord np = 4;
  MDArray< T > mda(tuple(ni,nj,nk,nm,nn,np));
  for (ord p=0; p < np; p++)
    for (ord n=0; n < nn; n++)
      for (ord m=0; m < nm; m++)
        for (ord k=0; k < nk; k++)
          for (ord j=0; j < nj; j++)
            for (ord i=0; i < ni; i++)
            {
              mda(i,j,k,m,n,p) = 5;
              TEST_EQUALITY(mda(i,j,k,m,n,p), 5);
            }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, indexing7D, T )
{
  // Note: this unit test exercises the variadic argument of the last
  // overloaded operator()
  typedef Domi::Ordinal ord;
  using Teuchos::tuple;
  ord ni = 2;
  ord nj = 2;
  ord nk = 3;
  ord nm = 3;
  ord nn = 4;
  ord np = 4;
  ord nq = 3;
  MDArray< T > mda(tuple(ni,nj,nk,nm,nn,np,nq));
  for (ord q=0; q < nq; q++)
    for (ord p=0; p < np; p++)
      for (ord n=0; n < nn; n++)
        for (ord m=0; m < nm; m++)
          for (ord k=0; k < nk; k++)
            for (ord j=0; j < nj; j++)
              for (ord i=0; i < ni; i++)
              {
                mda(i,j,k,m,n,p,q) = 6;
                TEST_EQUALITY(mda(i,j,k,m,n,p,q), 6);
              }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, squareBracketOrdinal, T )
{
  MDArray< T > a = generateMDArray< T >(3,4);

  MDArrayView< T > view = a[0][Slice()];
  TEST_EQUALITY_CONST(view(0), 0);
  TEST_EQUALITY_CONST(view(1), 3);
  TEST_EQUALITY_CONST(view(2), 6);
  TEST_EQUALITY_CONST(view(3), 9);

  view = a[1][Slice()];
  TEST_EQUALITY_CONST(view(0),  1);
  TEST_EQUALITY_CONST(view(1),  4);
  TEST_EQUALITY_CONST(view(2),  7);
  TEST_EQUALITY_CONST(view(3), 10);

  view = a[2][Slice()];
  TEST_EQUALITY_CONST(view(0),  2);
  TEST_EQUALITY_CONST(view(1),  5);
  TEST_EQUALITY_CONST(view(2),  8);
  TEST_EQUALITY_CONST(view(3), 11);

  // This tests a bug in computing the private _next_axis data member.
  // If incorrect, we'll get an exception here.
  view = a[Slice()][0];
  MDArrayView< T > view2 = view[1];

  view = a[1][1];
  TEST_EQUALITY_CONST(view(0), 4);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, squareBracketOrdinalCOrder, T )
{
  MDArray< T > a = generateMDArray< T >(3,4, Domi::C_ORDER);
  MDArrayView< T > view = a[0][Slice()];
  TEST_EQUALITY_CONST(view(0), 0);
  TEST_EQUALITY_CONST(view(1), 3);
  TEST_EQUALITY_CONST(view(2), 6);
  TEST_EQUALITY_CONST(view(3), 9);

  view = a[1][Slice()];
  TEST_EQUALITY_CONST(view(0),  1);
  TEST_EQUALITY_CONST(view(1),  4);
  TEST_EQUALITY_CONST(view(2),  7);
  TEST_EQUALITY_CONST(view(3), 10);

  view = a[2][Slice()];
  TEST_EQUALITY_CONST(view(0),  2);
  TEST_EQUALITY_CONST(view(1),  5);
  TEST_EQUALITY_CONST(view(2),  8);
  TEST_EQUALITY_CONST(view(3), 11);

  view = a[1][1];
  TEST_EQUALITY_CONST(view(0), 4);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, squareBracketSlice1, T )
{
  MDArray< T > a = generateMDArray< T >(4,5);
  MDArrayView< T > view = a[Slice(1,-1)][Slice(1,-1)];
  TEST_EQUALITY_CONST(view(0,0),  5);
  TEST_EQUALITY_CONST(view(0,1),  9);
  TEST_EQUALITY_CONST(view(0,2), 13);
  TEST_EQUALITY_CONST(view(1,0),  6);
  TEST_EQUALITY_CONST(view(1,1), 10);
  TEST_EQUALITY_CONST(view(1,2), 14);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, squareBracketSlice2, T )
{
  MDArray< T > a = generateMDArray< T >(4,3);
  MDArrayView< T > view = a[Slice(1,-1)][Slice()];
  TEST_EQUALITY_CONST(view(0,0),  1);
  TEST_EQUALITY_CONST(view(0,1),  5);
  TEST_EQUALITY_CONST(view(0,2),  9);
  TEST_EQUALITY_CONST(view(1,0),  2);
  TEST_EQUALITY_CONST(view(1,1),  6);
  TEST_EQUALITY_CONST(view(1,2), 10);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, rangeError, T )
{
  MDArray< T > a = generateMDArray< T >(3,4);
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEST_THROW(a(3,3), RangeError);
  TEST_THROW(a(0,4), RangeError);
#else
  a(0,0);  // Prevent unused variable warning
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, rangeErrorCOrder, T )
{
  MDArray< T > a = generateMDArray< T >(3,4,Domi::C_ORDER);
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEST_THROW(a(3,3), RangeError);
  TEST_THROW(a(0,4), RangeError);
#else
  a(0,0);  // Prevent unused variable warning
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, assign, T )
{
  typedef Domi::Ordinal ord;
  MDArray< T > a = generateMDArray< T >(3,4,5);
  a.assign(-1);
  for (ord k = 0; k < 5; ++k)
    for (ord j = 0; j < 4; ++j)
      for (ord i = 0; i < 3; ++i)
        TEST_EQUALITY_CONST(a(i,j,k), -1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, legalAt, T )
{
  MDArray< T > a = generateMDArray< T >(2,3);
  TEST_EQUALITY_CONST(a.at(0,0), 0);
  TEST_EQUALITY_CONST(a.at(1,0), 1);
  TEST_EQUALITY_CONST(a.at(0,1), 2);
  TEST_EQUALITY_CONST(a.at(1,1), 3);
  TEST_EQUALITY_CONST(a.at(0,2), 4);
  TEST_EQUALITY_CONST(a.at(1,2), 5);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, illegalAt, T )
{
  MDArray< T > a = generateMDArray< T >(2,3);
  TEST_THROW(a.at(2,0), RangeError);
  TEST_THROW(a.at(0,3), RangeError);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, clearEmpty, T )
{
  MDArray< T > a = generateMDArray< T >(10,10);
  TEST_EQUALITY_CONST(a.numDims()   ,  2);
  TEST_EQUALITY_CONST(a.dimension(0), 10);
  TEST_EQUALITY_CONST(a.dimension(1), 10);
  TEUCHOS_ASSERT(!a.empty());
  a.clear();
  TEST_EQUALITY_CONST(a.numDims()   , 1);
  TEST_EQUALITY_CONST(a.dimension(0), 0);
  TEUCHOS_ASSERT(a.empty());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, max_size, T )
{
  MDArray< T > a = generateMDArray< T >(5, 12);
  TEST_COMPARE( a.size(), <=, a.max_size() );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, resize, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArray< T > a = generateMDArray< T >(11,4);
  TEST_EQUALITY_CONST(a.numDims()   ,  2);
  TEST_EQUALITY_CONST(a.dimension(0), 11);
  TEST_EQUALITY_CONST(a.dimension(1),  4);
  TEUCHOS_ASSERT(!a.empty());
  a.resize(tuple< dim_type >(5,9,3));
  TEST_EQUALITY_CONST(a.numDims()   , 3);
  TEST_EQUALITY_CONST(a.dimension(0), 5);
  TEST_EQUALITY_CONST(a.dimension(1), 9);
  TEST_EQUALITY_CONST(a.dimension(2), 3);
  TEUCHOS_ASSERT(!a.empty());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, equality, T )
{
  MDArray< T > a = generateMDArray< T >(2,3,4);
  MDArray< T > b = generateMDArray< T >(2,3,4);
  TEST_EQUALITY(a, b)
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, inequality, T )
{
  MDArray< T > a = generateMDArray< T >(2,3,4);
  MDArray< T > b = generateMDArray< T >(2,3,4);
  b(1,1,1) = -1;
  TEST_INEQUALITY(a, b)
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, swap, T )
{
  MDArray< T > a = generateMDArray< T >(2,3);
  MDArray< T > b = generateMDArray< T >(4,5);
  TEST_EQUALITY_CONST(a.numDims()   , 2);
  TEST_EQUALITY_CONST(a.dimension(0), 2);
  TEST_EQUALITY_CONST(a.dimension(1), 3);
  TEST_EQUALITY_CONST(b.numDims()   , 2);
  TEST_EQUALITY_CONST(b.dimension(0), 4);
  TEST_EQUALITY_CONST(b.dimension(1), 5);
  a.swap(b);
  TEST_EQUALITY_CONST(a.numDims()   , 2);
  TEST_EQUALITY_CONST(a.dimension(0), 4);
  TEST_EQUALITY_CONST(a.dimension(1), 5);
  TEST_EQUALITY_CONST(b.numDims()   , 2);
  TEST_EQUALITY_CONST(b.dimension(0), 2);
  TEST_EQUALITY_CONST(b.dimension(1), 3);
  swap(a,b);
  TEST_EQUALITY_CONST(a.numDims()   , 2);
  TEST_EQUALITY_CONST(a.dimension(0), 2);
  TEST_EQUALITY_CONST(a.dimension(1), 3);
  TEST_EQUALITY_CONST(b.numDims()   , 2);
  TEST_EQUALITY_CONST(b.dimension(0), 4);
  TEST_EQUALITY_CONST(b.dimension(1), 5);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, toStringNull, T )
{
  MDArray< T > a;
  TEST_EQUALITY_CONST(a.toString(), "[]");
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, toString1D, T )
{
  typedef typename Domi::dim_type dim_type;
  T val = 3;
  MDArray< T > a(tuple< dim_type >(3), val);
  TEST_EQUALITY_CONST(a.toString(), "[3, 3, 3]");
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, toString2D, T )
{
  MDArray< T > a = generateMDArray< T >(2,3);
  TEST_EQUALITY_CONST(a.toString(), "[[0, 2, 4],\n [1, 3, 5]]");
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArray, toString3D, T )
{
  MDArray< T > a = generateMDArray< T >(2,3,2);
  TEST_EQUALITY_CONST(a.toString(), "[[[0, 1],\n  [2, 3],\n  [4, 5]],\n [[6, 7],\n  [8, 9],\n  [10, 11]]]");
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
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, defaultConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, simpleConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, dimsAndValConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, dimsAndOrderConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, dimsValAndOrderConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, copyConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, mdArrayViewConstructor0, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, mdArrayViewConstructor1, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, equalOperatorDiffLayout, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, equalOperatorMDArrayView, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, inequalityOperator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, indexing4D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, indexing5D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, indexing6D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, indexing7D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, squareBracketOrdinal, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, squareBracketOrdinalCOrder, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, squareBracketSlice1, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, squareBracketSlice2, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, rangeError, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, rangeErrorCOrder, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, assign, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, legalAt, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, illegalAt, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, clearEmpty, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, max_size, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, resize, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, equality, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, inequality, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, swap, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, toStringNull, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, toString1D, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, toString2D, T) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArray, toString3D, T) \
  DEBUG_UNIT_TEST_GROUP( T )

UNIT_TEST_GROUP(int)
#if 1
UNIT_TEST_GROUP(long)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)
#endif

}  // namespace
