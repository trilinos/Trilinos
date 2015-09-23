/*
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
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Domi_Slice.hpp"

namespace
{

using Domi::Slice;
const Domi::Ordinal & Default = Domi::Slice::Default;

TEUCHOS_UNIT_TEST( Slice, defaultConstructor )
{
  Slice slc;
  TEST_EQUALITY_CONST(slc.start(), 0      );
  TEST_EQUALITY_CONST(slc.stop(),  Default);
  TEST_EQUALITY_CONST(slc.step(),  1      );
}

TEUCHOS_UNIT_TEST( Slice, stopConstructor )
{
  Slice slc(5);
  TEST_EQUALITY_CONST(slc.start(), 0);
  TEST_EQUALITY_CONST(slc.stop(),  5);
  TEST_EQUALITY_CONST(slc.step(),  1);
}

TEUCHOS_UNIT_TEST( Slice, startStopConstructor1 )
{
  Slice slc(2,6);
  TEST_EQUALITY_CONST(slc.start(), 2);
  TEST_EQUALITY_CONST(slc.stop(),  6);
  TEST_EQUALITY_CONST(slc.step(),  1);
}

TEUCHOS_UNIT_TEST( Slice, startStopConstructor2 )
{
  Slice slc(Default,9);
  TEST_EQUALITY_CONST(slc.start(), 0);
  TEST_EQUALITY_CONST(slc.stop(),  9);
  TEST_EQUALITY_CONST(slc.step(),  1);
}

TEUCHOS_UNIT_TEST( Slice, startStopConstructor3 )
{
  Slice slc(1,Default);
  TEST_EQUALITY_CONST(slc.start(), 1      );
  TEST_EQUALITY_CONST(slc.stop(),  Default);
  TEST_EQUALITY_CONST(slc.step(),  1      );
}

TEUCHOS_UNIT_TEST( Slice, startStopStepConstructor1 )
{
  Slice slc(1,10,2);
  TEST_EQUALITY_CONST(slc.start(),  1);
  TEST_EQUALITY_CONST(slc.stop(),  10);
  TEST_EQUALITY_CONST(slc.step(),   2);
}

TEUCHOS_UNIT_TEST( Slice, startStopStepConstructor2 )
{
  Slice slc(5,10,Default);
  TEST_EQUALITY_CONST(slc.start(),  5);
  TEST_EQUALITY_CONST(slc.stop(),  10);
  TEST_EQUALITY_CONST(slc.step(),   1);
}

TEUCHOS_UNIT_TEST( Slice, startStopStepConstructor3 )
{
  Slice slc(1,Default,-1);
  TEST_EQUALITY_CONST(slc.start(), 1);
  TEST_EQUALITY_CONST(slc.stop(),  0);
  TEST_EQUALITY_CONST(slc.step(), -1);
}

TEUCHOS_UNIT_TEST( Slice, startStopStepConstructor4 )
{
  Slice slc(Default,10,-2);
  TEST_EQUALITY_CONST(slc.start(), Default);
  TEST_EQUALITY_CONST(slc.stop(),  10     );
  TEST_EQUALITY_CONST(slc.step(),  -2     );
}

TEUCHOS_UNIT_TEST( Slice, startStopStepConstructor5 )
{
  Slice slc(3,Default,Default);
  TEST_EQUALITY_CONST(slc.start(), 3      );
  TEST_EQUALITY_CONST(slc.stop(),  Default);
  TEST_EQUALITY_CONST(slc.step(),  1      );
}

TEUCHOS_UNIT_TEST( Slice, startStopStepConstructor6 )
{
  Slice slc(Default,Default,-1);
  TEST_EQUALITY_CONST(slc.start(), Default);
  TEST_EQUALITY_CONST(slc.stop(),   0     );
  TEST_EQUALITY_CONST(slc.step(),  -1     );
}

TEUCHOS_UNIT_TEST( Slice, startStopStepConstructor7 )
{
  Slice slc(Default,6,Default);
  TEST_EQUALITY_CONST(slc.start(), 0);
  TEST_EQUALITY_CONST(slc.stop(),  6);
  TEST_EQUALITY_CONST(slc.step(),  1);
}

TEUCHOS_UNIT_TEST( Slice, step0Constructor )
{
  TEST_THROW(Slice(1,2,0), std::invalid_argument);
}

TEUCHOS_UNIT_TEST( Slice, bounds1 )
{
  Slice slc(3, 14);
  Slice bounds = slc.bounds(20);
  TEST_EQUALITY_CONST(bounds.start(),  3);
  TEST_EQUALITY_CONST(bounds.stop(),  14);
  TEST_EQUALITY_CONST(bounds.step(),   1);
}

TEUCHOS_UNIT_TEST( Slice, bounds2 )
{
  Slice slc(3, 20);
  Slice bounds = slc.bounds(14);
  TEST_EQUALITY_CONST(bounds.start(),  3);
  TEST_EQUALITY_CONST(bounds.stop(),  14);
  TEST_EQUALITY_CONST(bounds.step(),   1);
}

TEUCHOS_UNIT_TEST( Slice, bounds3 )
{
  Slice slc(Default, 14, -1);
  Slice bounds = slc.bounds(20);
  TEST_EQUALITY_CONST(bounds.start(), 19);
  TEST_EQUALITY_CONST(bounds.stop(),  14);
  TEST_EQUALITY_CONST(bounds.step(),  -1);
}

TEUCHOS_UNIT_TEST( Slice, bounds4 )
{
  Slice slc(3, Default);
  Slice bounds = slc.bounds(20);
  TEST_EQUALITY_CONST(bounds.start(),  3);
  TEST_EQUALITY_CONST(bounds.stop(),  20);
  TEST_EQUALITY_CONST(bounds.step(),   1);
}

TEUCHOS_UNIT_TEST( Slice, bounds5 )
{
  Slice slc(3, -1);
  Slice bounds = slc.bounds(20);
  TEST_EQUALITY_CONST(bounds.start(),  3);
  TEST_EQUALITY_CONST(bounds.stop(),  19);
  TEST_EQUALITY_CONST(bounds.step(),   1);
}

TEUCHOS_UNIT_TEST( Slice, bounds6 )
{
  Slice slc(-2, 2, -1);
  Slice bounds = slc.bounds(20);
  TEST_EQUALITY_CONST(bounds.start(), 18);
  TEST_EQUALITY_CONST(bounds.stop(),   2);
  TEST_EQUALITY_CONST(bounds.step(),  -1);
}

TEUCHOS_UNIT_TEST( Slice, bounds7 )
{
  // This one is a little tricky.  I set up a slice from [0:20) with a
  // step() of 2.  When I ask for bounds, I cut it down to [0:11) with a
  // step() of 2.  This corresponds to iterates {0,2,4,6,8,10}.  By
  // returning a stop index of 12 (instead of 11), I can code the
  // following:
  //
  //     for (i=bounds.start(); i != bounds.stop(); i += bounds.step()) ...
  //
  // and it will work, because i == 12 at the stopping condition.
  // Note that I want to use (i != bounds.stop()) as my stopping
  // condition so that I can handle both positive and negative step()
  // sizes.

  Slice slc(0, 20, 2);
  Slice bounds = slc.bounds(11);
  TEST_EQUALITY_CONST(bounds.start(),  0);
  TEST_EQUALITY_CONST(bounds.stop(),  12);
  TEST_EQUALITY_CONST(bounds.step(),   2);
}

TEUCHOS_UNIT_TEST( Slice, toString0 )
{
  Slice slc;
  TEST_EQUALITY_CONST(slc.toString(), "[:]");
}

TEUCHOS_UNIT_TEST( Slice, toString1 )
{
  Slice slc(Default);
  TEST_EQUALITY_CONST(slc.toString(), "[:]");
}

TEUCHOS_UNIT_TEST( Slice, toString2 )
{
  Slice slc(10);
  TEST_EQUALITY_CONST(slc.toString(), "[:10]");
}

TEUCHOS_UNIT_TEST( Slice, toString3 )
{
  Slice slc(-1);
  TEST_EQUALITY_CONST(slc.toString(), "[:-1]");
}

TEUCHOS_UNIT_TEST( Slice, toString4 )
{
  Slice slc(1,-1);
  TEST_EQUALITY_CONST(slc.toString(), "[1:-1]");
}

TEUCHOS_UNIT_TEST( Slice, toString5 )
{
  Slice slc(2,Default);
  TEST_EQUALITY_CONST(slc.toString(), "[2:]");
}

TEUCHOS_UNIT_TEST( Slice, toString6 )
{
  Slice slc(3,14);
  TEST_EQUALITY_CONST(slc.toString(), "[3:14]");
}

TEUCHOS_UNIT_TEST( Slice, toString7 )
{
  Slice slc(Default,Default,-1);
  TEST_EQUALITY_CONST(slc.toString(), "[::-1]");
}

TEUCHOS_UNIT_TEST( Slice, toString8 )
{
  Slice slc(0,Default,2);
  TEST_EQUALITY_CONST(slc.toString(), "[::2]");
}

TEUCHOS_UNIT_TEST( Slice, toString9 )
{
  Slice slc(1,-1,2);
  TEST_EQUALITY_CONST(slc.toString(), "[1:-1:2]");
}

TEUCHOS_UNIT_TEST( Slice, streamOut )
{
  Slice slc(10,100,2);
  std::stringstream ss;
  ss << slc;
  TEST_EQUALITY_CONST(ss.str(), "[10:100:2]");
}

}  // namespace
