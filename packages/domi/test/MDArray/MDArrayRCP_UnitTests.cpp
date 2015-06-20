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

#include "Domi_MDArrayRCP.hpp"
#include "Domi_MDArray.hpp"
#include "MDArray_UnitTest_helpers.hpp"

namespace
{

using Teuchos::tuple;
using Teuchos::Array;
using Domi::MDArray;
using Domi::MDArrayRCP;
using Domi::MDArrayView;
using Domi::Slice;
using Domi::RangeError;
using MDArrayUnitTestHelpers::generateMDArray;
using MDArrayUnitTestHelpers::generateMDArrayRCP;

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, defaultConstructor, T )
{
  MDArrayRCP< T > mdar;
  TEST_EQUALITY_CONST(mdar.numDims()   , 1);
  TEST_EQUALITY_CONST(mdar.dimension(0), 0);
  TEST_EQUALITY_CONST(mdar.size()      , 0);
  TEST_EQUALITY_CONST(mdar.strides()[0], 1);
  TEUCHOS_ASSERT(mdar.empty());
  TEUCHOS_ASSERT(mdar.is_null());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, arrayViewDimsConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  Array< T > a(60);
  MDArrayRCP< T > mdar(a,tuple< dim_type >(3,4,5));
  TEST_EQUALITY(mdar.numDims()   ,  3);
  TEST_EQUALITY(mdar.dimension(0),  3);
  TEST_EQUALITY(mdar.dimension(1),  4);
  TEST_EQUALITY(mdar.dimension(2),  5);
  TEST_EQUALITY(mdar.size()      , 60);
  TEST_EQUALITY(mdar.strides()[0],  1);
  TEST_EQUALITY(mdar.strides()[1],  3);
  TEST_EQUALITY(mdar.strides()[2], 12);
  TEST_EQUALITY(mdar.layout(), Domi::DEFAULT_ORDER);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, arrayViewDimsConstructorBad, T )
{
  typedef typename Domi::dim_type dim_type;
  Array< T > a(50);
  TEST_THROW(MDArrayRCP< T > mdar(a,tuple< dim_type >(8,8)), RangeError);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, arrayViewDimsOrderConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  Array< T > a(60);
  MDArrayRCP< T > mdar(a,tuple< dim_type >(3,4,5),Domi::C_ORDER);
  TEST_EQUALITY(mdar.numDims()   ,  3);
  TEST_EQUALITY(mdar.dimension(0),  3);
  TEST_EQUALITY(mdar.dimension(1),  4);
  TEST_EQUALITY(mdar.dimension(2),  5);
  TEST_EQUALITY(mdar.size()      , 60);
  TEST_EQUALITY(mdar.strides()[0], 20);
  TEST_EQUALITY(mdar.strides()[1],  5);
  TEST_EQUALITY(mdar.strides()[2],  1);
  TEST_EQUALITY(mdar.layout(), Domi::C_ORDER);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, dimsConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > mdar(tuple< dim_type >(3,4,5));
  TEST_EQUALITY(mdar.numDims()   ,  3);
  TEST_EQUALITY(mdar.dimension(0),  3);
  TEST_EQUALITY(mdar.dimension(1),  4);
  TEST_EQUALITY(mdar.dimension(2),  5);
  TEST_EQUALITY(mdar.size()      , 60);
  TEST_EQUALITY(mdar.strides()[0],  1);
  TEST_EQUALITY(mdar.strides()[1],  3);
  TEST_EQUALITY(mdar.strides()[2], 12);
  TEST_EQUALITY(mdar.layout(), Domi::DEFAULT_ORDER);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, dimsValConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > mdar(tuple< dim_type >(3,4), 12);
  TEST_EQUALITY(mdar.numDims()   ,  2);
  TEST_EQUALITY(mdar.dimension(0),  3);
  TEST_EQUALITY(mdar.dimension(1),  4);
  TEST_EQUALITY(mdar.size()      , 12);
  TEST_EQUALITY(mdar.strides()[0],  1);
  TEST_EQUALITY(mdar.strides()[1],  3);
  TEST_EQUALITY(mdar.layout(), Domi::DEFAULT_ORDER);
  TEST_EQUALITY(mdar(0,0), 12);
  TEST_EQUALITY(mdar(1,0), 12);
  TEST_EQUALITY(mdar(2,0), 12);
  TEST_EQUALITY(mdar(0,1), 12);
  TEST_EQUALITY(mdar(1,1), 12);
  TEST_EQUALITY(mdar(2,1), 12);
  TEST_EQUALITY(mdar(0,2), 12);
  TEST_EQUALITY(mdar(1,2), 12);
  TEST_EQUALITY(mdar(2,2), 12);
  TEST_EQUALITY(mdar(0,3), 12);
  TEST_EQUALITY(mdar(1,3), 12);
  TEST_EQUALITY(mdar(2,3), 12);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, copyConstructor, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > mdar1(tuple< dim_type >(2,3), 2);
  MDArrayRCP< T > mdar2(mdar1);
  TEST_EQUALITY(mdar1.numDims(),    mdar2.numDims()   );
  TEST_EQUALITY(mdar1.dimension(0), mdar2.dimension(0));
  TEST_EQUALITY(mdar1.dimension(1), mdar2.dimension(1));
  TEST_EQUALITY(mdar1.size(),       mdar2.size()      );
  TEST_EQUALITY(mdar1.strides()[0], mdar2.strides()[0]);
  TEST_EQUALITY(mdar1.strides()[1], mdar2.strides()[1]);
  TEST_EQUALITY(mdar1.layout(),     mdar2.layout()    );
  TEST_EQUALITY(mdar1(0,0),         mdar2(0,0)        );
  TEST_EQUALITY(mdar1(1,0),         mdar2(1,0)        );
  TEST_EQUALITY(mdar1(0,1),         mdar2(0,1)        );
  TEST_EQUALITY(mdar1(1,1),         mdar2(1,1)        );
  TEST_EQUALITY(mdar1(0,2),         mdar2(0,2)        );
  TEST_EQUALITY(mdar1(1,2),         mdar2(1,2)        );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, equalOperator, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(3,4);
  MDArrayRCP< T > b = generateMDArrayRCP< T >(3,4);
  TEUCHOS_ASSERT( a == b );
  TEUCHOS_ASSERT( b == a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, equalOperatorMDArray, T )
{
  MDArray< T >    a = generateMDArray< T >(4,5);
  MDArrayRCP< T > b = generateMDArrayRCP< T >(4,5);
  TEUCHOS_ASSERT( a == b );
  TEUCHOS_ASSERT( b == a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, equalOperatorMDArrayView, T )
{
  MDArrayRCP< T >  a = generateMDArrayRCP< T >(2,6);
  MDArrayView< T > b = a();
  TEUCHOS_ASSERT( a == b );
  TEUCHOS_ASSERT( b == a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, inequalityOperator, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(3,4);
  MDArrayRCP< T > b = generateMDArrayRCP< T >(3,4);
  b(2,2) = -1;
  TEUCHOS_ASSERT( a != b );
  TEUCHOS_ASSERT( b != a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, inequalityOperatorMDArray, T )
{
  MDArray< T >    a = generateMDArray< T >(5,3);
  MDArrayRCP< T > b = generateMDArrayRCP< T >(5,3);
  a(4,2) = -1;
  TEUCHOS_ASSERT( a != b );
  TEUCHOS_ASSERT( b != a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, inequalityOperatorMDArrayView, T )
{
  MDArrayRCP< T >  a = generateMDArrayRCP< T >(2,4);
  MDArrayView< T > b = a[Slice()][Slice(1,3)];
  TEUCHOS_ASSERT( a != b );
  TEUCHOS_ASSERT( b != a );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, assignmentOperator, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > mdar1(tuple< dim_type >(2,2,2), 8);
  MDArrayRCP< T > mdar2 = mdar1;
  TEST_EQUALITY(mdar1.numDims(),    mdar2.numDims()   );
  TEST_EQUALITY(mdar1.dimension(0), mdar2.dimension(0));
  TEST_EQUALITY(mdar1.dimension(1), mdar2.dimension(1));
  TEST_EQUALITY(mdar1.dimension(2), mdar2.dimension(2));
  TEST_EQUALITY(mdar1.size(),       mdar2.size()      );
  TEST_EQUALITY(mdar1.strides()[0], mdar2.strides()[0]);
  TEST_EQUALITY(mdar1.strides()[1], mdar2.strides()[1]);
  TEST_EQUALITY(mdar1.strides()[2], mdar2.strides()[2]);
  TEST_EQUALITY(mdar1.layout(),     mdar2.layout()    );
  TEST_EQUALITY(mdar1(0,0,0),       mdar2(0,0,0)      );
  TEST_EQUALITY(mdar1(1,0,0),       mdar2(1,0,0)      );
  TEST_EQUALITY(mdar1(0,1,0),       mdar2(0,1,0)      );
  TEST_EQUALITY(mdar1(1,1,0),       mdar2(1,1,0)      );
  TEST_EQUALITY(mdar1(0,0,1),       mdar2(0,0,1)      );
  TEST_EQUALITY(mdar1(1,0,1),       mdar2(1,0,1)      );
  TEST_EQUALITY(mdar1(0,1,1),       mdar2(0,1,1)      );
  TEST_EQUALITY(mdar1(1,1,1),       mdar2(1,1,1)      );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, indexing4D, T )
{
  typedef Domi::Ordinal ord;
  ord ni = 2;
  ord nj = 2;
  ord nk = 3;
  ord nm = 3;
  MDArrayRCP< T > mdarcp(tuple(ni,nj,nk,nm));
  MDArrayView< T > mdav = mdarcp();
  for (ord m=0; m < nm; m++)
    for (ord k=0; k < nk; k++)
      for (ord j=0; j < nj; j++)
        for (ord i=0; i < ni; i++)
        {
          mdarcp(i,j,k,m) = 3;
          TEST_EQUALITY(mdarcp(i,j,k,m), 3);
        }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, indexing5D, T )
{
  typedef Domi::Ordinal ord;
  ord ni = 2;
  ord nj = 2;
  ord nk = 3;
  ord nm = 3;
  ord nn = 4;
  MDArrayRCP< T > mdarcp(tuple(ni,nj,nk,nm,nn));
  for (ord n=0; n < nn; n++)
    for (ord m=0; m < nm; m++)
      for (ord k=0; k < nk; k++)
        for (ord j=0; j < nj; j++)
          for (ord i=0; i < ni; i++)
          {
            mdarcp(i,j,k,m,n) = 4;
            TEST_EQUALITY(mdarcp(i,j,k,m,n), 4);
          }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, indexing6D, T )
{
  typedef Domi::Ordinal ord;
  ord ni = 2;
  ord nj = 2;
  ord nk = 3;
  ord nm = 3;
  ord nn = 4;
  ord np = 4;
  MDArrayRCP< T > mdarcp(tuple(ni,nj,nk,nm,nn,np));
  for (ord p=0; p < np; p++)
    for (ord n=0; n < nn; n++)
      for (ord m=0; m < nm; m++)
        for (ord k=0; k < nk; k++)
          for (ord j=0; j < nj; j++)
            for (ord i=0; i < ni; i++)
            {
              mdarcp(i,j,k,m,n,p) = 5;
              TEST_EQUALITY(mdarcp(i,j,k,m,n,p), 5);
            }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, indexing7D, T )
{
  // Note: this unit test exercises the variadic argument of the last
  // overloaded operator()
  typedef Domi::Ordinal ord;
  ord ni = 2;
  ord nj = 2;
  ord nk = 3;
  ord nm = 3;
  ord nn = 4;
  ord np = 4;
  ord nq = 3;
  MDArrayRCP< T > mdarcp(tuple(ni,nj,nk,nm,nn,np,nq));
  for (ord q=0; q < nq; q++)
    for (ord p=0; p < np; p++)
      for (ord n=0; n < nn; n++)
        for (ord m=0; m < nm; m++)
          for (ord k=0; k < nk; k++)
            for (ord j=0; j < nj; j++)
              for (ord i=0; i < ni; i++)
              {
                mdarcp(i,j,k,m,n,p,q) = 6;
                TEST_EQUALITY(mdarcp(i,j,k,m,n,p,q), 6);
              }
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, squareBracketOrdinal, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(3,4);
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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, squareBracketOrdinalCOrder, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(3,4, Domi::C_ORDER);
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

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, squareBracketSlice1, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(4,5);
  MDArrayView< T > view = a[Slice(1,-1)][Slice(1,-1)];
  TEST_EQUALITY_CONST(view(0,0),  5);
  TEST_EQUALITY_CONST(view(0,1),  9);
  TEST_EQUALITY_CONST(view(0,2), 13);
  TEST_EQUALITY_CONST(view(1,0),  6);
  TEST_EQUALITY_CONST(view(1,1), 10);
  TEST_EQUALITY_CONST(view(1,2), 14);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, squareBracketSlice2, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(4,3);
  MDArrayView< T > view = a[Slice(1,-1)][Slice()];
  TEST_EQUALITY_CONST(view(0,0),  1);
  TEST_EQUALITY_CONST(view(0,1),  5);
  TEST_EQUALITY_CONST(view(0,2),  9);
  TEST_EQUALITY_CONST(view(1,0),  2);
  TEST_EQUALITY_CONST(view(1,1),  6);
  TEST_EQUALITY_CONST(view(1,2), 10);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, rangeError, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(3,4);
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEST_THROW(a(3,3), RangeError);
  TEST_THROW(a(0,4), RangeError);
#else
  a(0,0);  // Prevent unused variable warning
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, rangeErrorCOrder, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(3,4,Domi::C_ORDER);
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEST_THROW(a(3,3), RangeError);
  TEST_THROW(a(0,4), RangeError);
#else
  a(0,0);  // Prevent unused variable warning
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, assign, T )
{
  typedef Domi::Ordinal ord;
  MDArrayRCP< T > a = generateMDArrayRCP< T >(3,4,5);
  a.assign(-1);
  for (ord k = 0; k < 5; ++k)
    for (ord j = 0; j < 4; ++j)
      for (ord i = 0; i < 3; ++i)
        TEST_EQUALITY_CONST(a(i,j,k), -1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, legalAt, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(2,3);
  TEST_EQUALITY_CONST(a.at(0,0), 0);
  TEST_EQUALITY_CONST(a.at(1,0), 1);
  TEST_EQUALITY_CONST(a.at(0,1), 2);
  TEST_EQUALITY_CONST(a.at(1,1), 3);
  TEST_EQUALITY_CONST(a.at(0,2), 4);
  TEST_EQUALITY_CONST(a.at(1,2), 5);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, illegalAt, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(2,3);
  TEST_THROW(a.at(2,0), RangeError);
  TEST_THROW(a.at(0,3), RangeError);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, clearEmpty, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(10,10);
  TEST_EQUALITY_CONST(a.numDims()   ,  2);
  TEST_EQUALITY_CONST(a.dimension(0), 10);
  TEST_EQUALITY_CONST(a.dimension(1), 10);
  TEUCHOS_ASSERT(!a.empty());
  a.clear();
  TEST_EQUALITY_CONST(a.numDims()   , 1);
  TEST_EQUALITY_CONST(a.dimension(0), 0);
  TEUCHOS_ASSERT(a.empty());
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, resize, T )
{
  typedef typename Domi::dim_type dim_type;
  MDArrayRCP< T > a = generateMDArrayRCP< T >(11,4);
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

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, equality, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(2,3,4);
  MDArrayRCP< T > b = generateMDArrayRCP< T >(2,3,4);
  TEST_EQUALITY(a, b)
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, inequality, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(2,3,4);
  MDArrayRCP< T > b = generateMDArrayRCP< T >(2,3,4);
  b(1,1,1) = -1;
  TEST_INEQUALITY(a, b)
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, toStringNull, T )
{
  MDArrayRCP< T > a;
  TEST_EQUALITY_CONST(a.toString(), "[]");
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, toString1D, T )
{
  typedef typename Domi::dim_type dim_type;
  T val = 3;
  MDArrayRCP< T > a(tuple< dim_type >(3), val);
  TEST_EQUALITY_CONST(a.toString(), "[3, 3, 3]");
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, toString2D, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(2,3);
  TEST_EQUALITY_CONST(a.toString(), "[[0, 2, 4],\n [1, 3, 5]]");
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDArrayRCP, toString3D, T )
{
  MDArrayRCP< T > a = generateMDArrayRCP< T >(2,3,2);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, defaultConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, arrayViewDimsConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, arrayViewDimsConstructorBad, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, arrayViewDimsOrderConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, dimsConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, dimsValConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, copyConstructor, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, equalOperator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, equalOperatorMDArray, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, equalOperatorMDArrayView, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, inequalityOperator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, inequalityOperatorMDArray, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, inequalityOperatorMDArrayView, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, assignmentOperator, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, indexing4D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, indexing5D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, indexing6D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, indexing7D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, squareBracketOrdinal, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, squareBracketOrdinalCOrder, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, squareBracketSlice1, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, squareBracketSlice2, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, rangeError, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, rangeErrorCOrder, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, assign, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, legalAt, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, illegalAt, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, clearEmpty, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, resize, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, equality, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, inequality, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, toStringNull, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, toString1D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, toString2D, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDArrayRCP, toString3D, T ) \
  DEBUG_UNIT_TEST_GROUP( T )

UNIT_TEST_GROUP(int)
#if 1
UNIT_TEST_GROUP(long)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)
#endif

}  // End anonymous namespace
