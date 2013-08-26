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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Range1D.hpp"


namespace {


using Teuchos::Range1D;
typedef Teuchos::Ordinal Ordinal;


TEUCHOS_UNIT_TEST( Range1D, range_0_0 )
{
  ECHO(Range1D rng(0,0));
  TEST_EQUALITY_CONST(rng.lbound(), 0);
  TEST_EQUALITY_CONST(rng.ubound(), 0);;
  TEST_EQUALITY_CONST(rng.size(), 1);
  TEST_ASSERT(!rng.full_range());
  TEST_ASSERT(!rng.in_range(-1));
  TEST_ASSERT(rng.in_range(0));
  TEST_ASSERT(!rng.in_range(1));
  TEST_ASSERT(!rng.in_range(2));
}


TEUCHOS_UNIT_TEST( Range1D, range_1_2 )
{
  ECHO(Range1D rng(1,2));
  TEST_EQUALITY_CONST(rng.lbound(), 1);
  TEST_EQUALITY_CONST(rng.ubound(), 2);;
  TEST_EQUALITY_CONST(rng.size(), 2);
  TEST_ASSERT(!rng.full_range());
  TEST_ASSERT(!rng.in_range(-1));
  TEST_ASSERT(!rng.in_range(0));
  TEST_ASSERT(rng.in_range(1));
  TEST_ASSERT(rng.in_range(2));
  TEST_ASSERT(!rng.in_range(3));
}


TEUCHOS_UNIT_TEST( Range1D, range_full )
{
  ECHO(Range1D rng);
  TEST_EQUALITY_CONST(rng.lbound(), 0);
  TEST_EQUALITY_CONST(rng.ubound(), std::numeric_limits<Ordinal>::max()-1);
  TEST_EQUALITY_CONST(rng.size(), std::numeric_limits<Ordinal>::max());
  TEST_ASSERT(rng.full_range());
  TEST_ASSERT(!rng.in_range(-1));
  TEST_ASSERT(rng.in_range(0));
  TEST_ASSERT(rng.in_range(1));
  TEST_ASSERT(rng.in_range(2));
  TEST_ASSERT(rng.in_range(std::numeric_limits<Ordinal>::max()-1));
  TEST_ASSERT(!rng.in_range(std::numeric_limits<Ordinal>::max()));
}


TEUCHOS_UNIT_TEST( Range1D, range_invalid )
{
  ECHO(Range1D rng(Range1D::INVALID));
  TEST_EQUALITY_CONST(rng.lbound(), 0);
  TEST_EQUALITY_CONST(rng.ubound(), -2);
  TEST_EQUALITY_CONST(rng.size(), -1);
  TEST_ASSERT(!rng.full_range());
  TEST_ASSERT(!rng.in_range(-1));
  TEST_ASSERT(!rng.in_range(0));
  TEST_ASSERT(!rng.in_range(1));
}


TEUCHOS_UNIT_TEST( Range1D, range_0_m1 )
{
  ECHO(Range1D rng(0,-1));
  TEST_EQUALITY_CONST(rng.lbound(), 0);
  TEST_EQUALITY_CONST(rng.ubound(), -1);;
  TEST_EQUALITY_CONST(rng.size(), 0);
  TEST_ASSERT(!rng.full_range());
  TEST_ASSERT(!rng.in_range(-1));
  TEST_ASSERT(!rng.in_range(0));
  TEST_ASSERT(!rng.in_range(1));
}


TEUCHOS_UNIT_TEST( Range1D, range_1_0 )
{
  ECHO(Range1D rng(1,0));
  TEST_EQUALITY_CONST(rng.lbound(), 1);
  TEST_EQUALITY_CONST(rng.ubound(), 0);;
  TEST_EQUALITY_CONST(rng.size(), 0);
  TEST_ASSERT(!rng.full_range());
  TEST_ASSERT(!rng.in_range(0));
  TEST_ASSERT(!rng.in_range(1));
  TEST_ASSERT(!rng.in_range(2));
}


TEUCHOS_UNIT_TEST( Range1D, range_4_3 )
{
  ECHO(Range1D rng(4,3));
  TEST_EQUALITY_CONST(rng.lbound(), 4);
  TEST_EQUALITY_CONST(rng.ubound(), 3);
  TEST_EQUALITY_CONST(rng.size(), 0);
  TEST_ASSERT(!rng.full_range());
  TEST_ASSERT(!rng.in_range(2));
  TEST_ASSERT(!rng.in_range(3));
  TEST_ASSERT(!rng.in_range(4));
}


TEUCHOS_UNIT_TEST( Range1D, equalityOp )
{

  ECHO(Range1D rng1(0,1));

  TEST_EQUALITY(rng1, rng1);

  ECHO(Range1D rng2(0,1));
  TEST_EQUALITY(rng2, rng1);

  ECHO(Range1D rng3(0,2));
  TEST_INEQUALITY(rng3, rng1);

  ECHO(Range1D rng4(1,2));
  TEST_INEQUALITY(rng3, rng1);

  ECHO(Range1D rng5 = rng4 - 1);
  TEST_EQUALITY(rng5, rng1);

}


TEUCHOS_UNIT_TEST( Range1D, increment )
{
  ECHO(Range1D rng1(4,6));
  TEST_EQUALITY_CONST(rng1.lbound(), 4);
  TEST_EQUALITY_CONST(rng1.ubound(), 6);

  ECHO(rng1 += 3);
  TEST_EQUALITY_CONST(rng1.lbound(), 7);
  TEST_EQUALITY_CONST(rng1.ubound(), 9);

  ECHO(rng1 -= 1);
  TEST_EQUALITY_CONST(rng1.lbound(), 6);
  TEST_EQUALITY_CONST(rng1.ubound(), 8);

  ECHO(rng1 -= 6);
  TEST_EQUALITY_CONST(rng1.lbound(), 0);
  TEST_EQUALITY_CONST(rng1.ubound(), 2);

  ECHO(Range1D rng2 = Range1D(2,3) + 4);
  TEST_EQUALITY_CONST(rng2.lbound(), 6);
  TEST_EQUALITY_CONST(rng2.ubound(), 7);

  ECHO(Range1D rng2b = 4 + Range1D(2,3));
  TEST_EQUALITY_CONST(rng2b.lbound(), 6);
  TEST_EQUALITY_CONST(rng2b.ubound(), 7);

  ECHO(Range1D rng3 = Range1D(4,5) - 2);
  TEST_EQUALITY_CONST(rng3.lbound(), 2);
  TEST_EQUALITY_CONST(rng3.ubound(), 3);

  ECHO(Range1D rng4 = Range1D(4,4) - 4);
  TEST_EQUALITY_CONST(rng4.lbound(), 0);
  TEST_EQUALITY_CONST(rng4.ubound(), 0);

  ECHO(Range1D rng5 = Range1D(4,4) + (-4));
  TEST_EQUALITY_CONST(rng5.lbound(), 0);
  TEST_EQUALITY_CONST(rng5.ubound(), 0);
}


#ifdef TEUCHOS_DEBUG

TEUCHOS_UNIT_TEST( Range1D, outOfRange )
{
  TEST_THROW(Range1D(-1,-1), std::out_of_range);
  TEST_THROW(Range1D(-1,1), std::out_of_range);
  TEST_THROW(Range1D(2,0), std::out_of_range);
  TEST_THROW(Range1D(3,0), std::out_of_range);
  TEST_THROW(Range1D(5,3), std::out_of_range);
  TEST_THROW(Range1D(0,0)-1, std::out_of_range);
  TEST_THROW(Range1D(0,0)+(-1), std::out_of_range);
}

#endif // TEUCHOS_DEBUG

// ToDo: Test creating invalid ranges

// ToDo: Test invalid lower increment.



} // namespace
