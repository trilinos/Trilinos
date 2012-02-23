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


#include "Teuchos_FilteredIterator.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


template<typename T>
class SelectAll {
public:
  bool operator()(const T& x) const
    { return true; }
};


template<typename IntegralType>
class SelectEven {
public:
  bool operator()(const IntegralType& x) const
    { return ( (x % 2) == 0 ); }
};


template<typename IntegralType>
class SelectOdd {
public:
  bool operator()(const IntegralType& x) const
    { return ( (x % 2) != 0 ); }
};


} // namespace


namespace Teuchos {


TEUCHOS_UNIT_TEST( FilteredIterator, default_construct )
{
  FilteredIterator<int*, SelectAll<int> > itr;
  // There is just nothing that we can check for an uninitialized iterator!
}


TEUCHOS_UNIT_TEST( FilteredIterator, construct )
{
  typedef Array<int>::iterator itr_t;
  Array<int> a;
  a.push_back(1);
  FilteredIterator<itr_t,SelectAll<int> > itr(a.begin(), a.begin(), a.end());
  TEST_ITER_EQUALITY(itr.current(), a.begin());
  TEST_ITER_EQUALITY(itr.begin(), a.begin());
  TEST_ITER_EQUALITY(itr.end(), a.end());
  FilteredIterator<itr_t,SelectAll<int> > itr_end(a.end(), a.begin(), a.end());
  TEST_ITER_EQUALITY(itr_end.current(), a.end());
  TEST_ITER_EQUALITY(itr_end.begin(), a.begin());
  TEST_ITER_EQUALITY(itr_end.end(), a.end());
}


TEUCHOS_UNIT_TEST( FilteredIterator, deref )
{
  typedef Array<int>::iterator itr_t;
  Array<int> a;
  a.push_back(2);
  FilteredIterator<itr_t,SelectAll<int> > itr(a.begin(), a.begin(), a.end());
  TEST_EQUALITY_CONST(*itr, 2);
}


TEUCHOS_UNIT_TEST( FilteredIterator, member_access )
{
  typedef std::pair<int,int> value_t;
  typedef Array<value_t>::iterator itr_t;
  Array<value_t> a;
  a.push_back(std::make_pair(2, 4));
  FilteredIterator<itr_t,SelectAll<value_t> > itr(a.begin(), a.begin(), a.end());
  TEST_EQUALITY_CONST(itr->first, 2);
  TEST_EQUALITY_CONST(itr->second, 4);
}


TEUCHOS_UNIT_TEST( FilteredIterator, copy_construct_same )
{
  typedef Array<int>::iterator itr_t;
  Array<int> a;
  a.push_back(1);
  FilteredIterator<itr_t,SelectAll<int> > itr1(a.begin(), a.begin(), a.end());
  FilteredIterator<itr_t,SelectAll<int> > itr2(itr1);
  TEST_ITER_EQUALITY(itr2.current(), a.begin());
  TEST_ITER_EQUALITY(itr2.begin(), a.begin());
  TEST_ITER_EQUALITY(itr2.end(), a.end());
  TEST_EQUALITY(*itr1, *itr2);
}


TEUCHOS_UNIT_TEST( FilteredIterator, copy_construct_different )
{
  typedef Array<int>::iterator itr_t;
  typedef Array<int>::const_iterator citr_t;
  Array<int> a;
  a.push_back(1);
  FilteredIterator<itr_t,SelectAll<int> > itr1(a.begin(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectAll<int> > itr2(itr1);
  TEST_ITER_EQUALITY(itr2.current(), a.begin());
  TEST_ITER_EQUALITY(itr2.begin(), a.begin());
  TEST_ITER_EQUALITY(itr2.end(), a.end());
  TEST_EQUALITY(*itr1, *itr2);
}


TEUCHOS_UNIT_TEST( FilteredIterator, assign_same )
{
  typedef Array<int>::iterator itr_t;
  Array<int> a;
  a.push_back(1);
  FilteredIterator<itr_t,SelectAll<int> > itr1(a.begin(), a.begin(), a.end());
  FilteredIterator<itr_t,SelectAll<int> > itr2;
  itr2 = itr1;
  TEST_ITER_EQUALITY(itr2.current(), a.begin());
  TEST_ITER_EQUALITY(itr2.begin(), a.begin());
  TEST_ITER_EQUALITY(itr2.end(), a.end());
  TEST_EQUALITY(*itr1, *itr2);
}


TEUCHOS_UNIT_TEST( FilteredIterator, assign_different )
{
  typedef Array<int>::iterator itr_t;
  typedef Array<int>::const_iterator citr_t;
  Array<int> a;
  a.push_back(1);
  FilteredIterator<itr_t,SelectAll<int> > itr1(a.begin(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectAll<int> > itr2;
  itr2 = itr1;
  TEST_ITER_EQUALITY(itr2.current(), a.begin());
  TEST_ITER_EQUALITY(itr2.begin(), a.begin());
  TEST_ITER_EQUALITY(itr2.end(), a.end());
  TEST_EQUALITY(*itr1, *itr2);
}


TEUCHOS_UNIT_TEST( FilteredIterator, equality_operators_same )
{
  typedef Array<int>::iterator itr_t;
  Array<int> a;
  a.push_back(1);
  FilteredIterator<itr_t,SelectAll<int> > itr1(a.begin(), a.begin(), a.end());
  FilteredIterator<itr_t,SelectAll<int> > itr2(itr1);
  TEST_EQUALITY_CONST(itr2 == itr1, true);
  TEST_EQUALITY_CONST(itr2 != itr1, false);
}


TEUCHOS_UNIT_TEST( FilteredIterator, equality_operators_different )
{
  typedef Array<int>::iterator itr_t;
  Array<int> a;
  a.push_back(1);
  FilteredIterator<itr_t,SelectAll<int> > itr1(a.begin(), a.begin(), a.end());
  FilteredIterator<itr_t,SelectAll<int> > itr2(a.end(), a.begin(), a.end());
  TEST_EQUALITY_CONST(itr2 == itr1, false);
  TEST_EQUALITY_CONST(itr2 != itr1, true);
}


TEUCHOS_UNIT_TEST( FilteredIterator, pre_iterate_forward_no_filtering )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3);
  FilteredIterator<citr_t,SelectAll<int> > itr(a.begin(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectAll<int> > itr_end(a.end(), a.begin(), a.end());
  TEST_ITER_INEQUALITY(itr, itr_end);
  TEST_EQUALITY_CONST(*itr, 1);
  ECHO(++itr);
  TEST_EQUALITY_CONST(*itr, 2);
  ECHO(++itr);
  TEST_EQUALITY_CONST(*itr, 3);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, itr_end);
}


TEUCHOS_UNIT_TEST( FilteredIterator, post_iterate_forward_no_filtering )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3);
  FilteredIterator<citr_t,SelectAll<int> > itr(a.begin(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectAll<int> > itr_end(a.end(), a.begin(), a.end());
  TEST_ITER_INEQUALITY(itr, itr_end);
  ECHO(const int v0 = *itr++);
  TEST_EQUALITY_CONST(v0, 1);
  ECHO(const int v1 = *itr++);
  TEST_EQUALITY_CONST(v1, 2);
  ECHO(const int v2 = *itr++);
  TEST_EQUALITY_CONST(v2, 3);
  TEST_ITER_EQUALITY(itr, itr_end);
}


TEUCHOS_UNIT_TEST( FilteredIterator, pre_iterate_backward_no_filtering )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3);
  FilteredIterator<citr_t,SelectAll<int> > itr(a.end(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectAll<int> > itr_begin(a.begin(), a.begin(), a.end());
  ECHO(--itr);
  TEST_EQUALITY_CONST(*itr, 3);
  ECHO(--itr);
  TEST_EQUALITY_CONST(*itr, 2);
  ECHO(--itr);
  TEST_EQUALITY_CONST(*itr, 1);
  TEST_ITER_EQUALITY(itr, itr_begin);
}


TEUCHOS_UNIT_TEST( FilteredIterator, post_iterate_backward_no_filtering )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3);
  FilteredIterator<citr_t,SelectAll<int> > itr(a.end(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectAll<int> > itr_begin(a.begin(), a.begin(), a.end());
  ECHO(--itr);
  ECHO(const int v0 = *itr--);
  TEST_EQUALITY_CONST(v0, 3);
  ECHO(const int v1 = *itr--);
  TEST_EQUALITY_CONST(v1, 2);
  ECHO(const int v2 = *itr);
  TEST_EQUALITY_CONST(v2, 1);
  TEST_ITER_EQUALITY(itr, itr_begin);
}


TEUCHOS_UNIT_TEST( FilteredIterator, pre_iterate_forward_filter_even )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3, 4);
  FilteredIterator<citr_t,SelectEven<int> > itr(a.begin(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectEven<int> > itr_end(a.end(), a.begin(), a.end());
  TEST_ITER_INEQUALITY(itr, itr_end);
  TEST_EQUALITY_CONST(*itr, 2);
  ECHO(++itr);
  TEST_EQUALITY_CONST(*itr, 4);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, itr_end);
}


TEUCHOS_UNIT_TEST( FilteredIterator, pre_iterate_forward_filter_odd )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3, 4);
  FilteredIterator<citr_t,SelectOdd<int> > itr(a.begin(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectOdd<int> > itr_end(a.end(), a.begin(), a.end());
  TEST_ITER_INEQUALITY(itr, itr_end);
  TEST_EQUALITY_CONST(*itr, 1);
  ECHO(++itr);
  TEST_EQUALITY_CONST(*itr, 3);
  ECHO(++itr);
  TEST_ITER_EQUALITY(itr, itr_end);
}


TEUCHOS_UNIT_TEST( FilteredIterator, post_iterate_forward_filter_even )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3, 4);
  FilteredIterator<citr_t,SelectEven<int> > itr(a.begin(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectEven<int> > itr_end(a.end(), a.begin(), a.end());
  TEST_ITER_INEQUALITY(itr, itr_end);
  ECHO(const int v0 = *itr++);
  TEST_EQUALITY_CONST(v0, 2);
  ECHO(const int v1 = *itr++);
  TEST_EQUALITY_CONST(v1, 4);
  TEST_ITER_EQUALITY(itr, itr_end);
}


TEUCHOS_UNIT_TEST( FilteredIterator, post_iterate_forward_filter_odd )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3, 4);
  FilteredIterator<citr_t,SelectOdd<int> > itr(a.begin(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectOdd<int> > itr_end(a.end(), a.begin(), a.end());
  TEST_ITER_INEQUALITY(itr, itr_end);
  ECHO(const int v0 = *itr++);
  TEST_EQUALITY_CONST(v0, 1);
  ECHO(const int v1 = *itr++);
  TEST_EQUALITY_CONST(v1, 3);
  TEST_ITER_EQUALITY(itr, itr_end);
}


TEUCHOS_UNIT_TEST( FilteredIterator, pre_iterate_backward_filter_even )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3, 4);
  FilteredIterator<citr_t,SelectEven<int> > itr(a.end(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectEven<int> > itr_begin(a.begin(), a.begin(), a.end());
  ECHO(--itr);
  TEST_EQUALITY_CONST(*itr, 4);
  ECHO(--itr);
  TEST_EQUALITY_CONST(*itr, 2);
  TEST_ITER_EQUALITY(itr, itr_begin);
}


TEUCHOS_UNIT_TEST( FilteredIterator, pre_iterate_backward_filter_odd )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3, 4);
  FilteredIterator<citr_t,SelectOdd<int> > itr(a.end(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectOdd<int> > itr_begin(a.begin(), a.begin(), a.end());
  ECHO(--itr);
  TEST_EQUALITY_CONST(*itr, 3);
  ECHO(--itr);
  TEST_EQUALITY_CONST(*itr, 1);
  TEST_ITER_EQUALITY(itr, itr_begin);
}


TEUCHOS_UNIT_TEST( FilteredIterator, post_iterate_backward_filter_even )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3, 4);
  FilteredIterator<citr_t,SelectEven<int> > itr(a.end(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectEven<int> > itr_begin(a.begin(), a.begin(), a.end());
  ECHO(--itr);
  ECHO(const int v0 = *itr--);
  TEST_EQUALITY_CONST(v0, 4);
  ECHO(const int v1 = *itr);
  TEST_EQUALITY_CONST(v1, 2);
  TEST_ITER_EQUALITY(itr, itr_begin);
}


TEUCHOS_UNIT_TEST( FilteredIterator, post_iterate_backward_filter_odd )
{
  typedef Array<int>::const_iterator citr_t;
  Array<int> a = tuple<int>(1, 2, 3, 4);
  FilteredIterator<citr_t,SelectOdd<int> > itr(a.end(), a.begin(), a.end());
  FilteredIterator<citr_t,SelectOdd<int> > itr_begin(a.begin(), a.begin(), a.end());
  ECHO(--itr);
  ECHO(const int v0 = *itr--);
  TEST_EQUALITY_CONST(v0, 3);
  ECHO(const int v1 = *itr);
  TEST_EQUALITY_CONST(v1, 1);
  TEST_ITER_EQUALITY(itr, itr_begin);
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


TEUCHOS_UNIT_TEST( FilteredIterator, iterate_forward_past_end_throw )
{
  // Need to use an unchecked iterator to make sure class thows
  int a_raw = 1;
  FilteredIterator<int*,SelectAll<int> > itr_end((&a_raw)+1, &a_raw, (&a_raw)+1);
  FilteredIterator<int*,SelectAll<int> > itr = itr_end;
  TEST_THROW(++itr, RangeError);
  TEST_ITER_EQUALITY(itr, itr_end);
  TEST_THROW(itr++, RangeError);
  TEST_ITER_EQUALITY(itr, itr_end);
}


TEUCHOS_UNIT_TEST( FilteredIterator, iterate_backward_past_begin_throw )
{
  // Need to use an unchecked iterator to make sure class thows
  int a_raw = 1;
  FilteredIterator<int*,SelectAll<int> > itr_begin(&a_raw, &a_raw, (&a_raw)+1);
  FilteredIterator<int*,SelectAll<int> > itr = itr_begin;
  TEST_THROW(--itr, RangeError);
  TEST_ITER_EQUALITY(itr, itr_begin);
  TEST_THROW(itr--, RangeError);
  TEST_ITER_EQUALITY(itr, itr_begin);
}


#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


} // namespace Teuchos



