/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <fei_LinearDecomposition.hpp>
#include <iostream>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(LinearDecomposition, test0, T)
{
  int myproc = 3;
  int numprocs = 4;
  T lowest_global = 1;
  T highest_global = 100;

  fei::LinearDecomposition<T> ld1(myproc, numprocs, lowest_global, highest_global);
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(0), -1, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(5), 0, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(25), 0, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(26), 1, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(50), 1, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(51), 2, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(76), 3, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(100), 3, out, success)
  TEUCHOS_TEST_EQUALITY(ld1.which_proc(101), -1, out, success)

  lowest_global = 1;
  highest_global = 12;
  numprocs = 10;
  fei::LinearDecomposition<T> ld2(myproc, numprocs, lowest_global, highest_global);
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(0), -1, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(1), 0, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(2), 0, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(3), 1, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(5), 2, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(11), 8, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(12), 9, out, success)
  TEUCHOS_TEST_EQUALITY(ld2.which_proc(13), -1, out, success)
}

#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(LinearDecomposition,test0,TYPE)

typedef long int longint;
UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(longint)

}//namespace <anonymous>

