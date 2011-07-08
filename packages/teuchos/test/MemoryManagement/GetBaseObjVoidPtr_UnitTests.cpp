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

#include "Teuchos_getBaseObjVoidPtr.hpp"
#include "Teuchos_RCP.hpp"

#include "TestClasses.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace Teuchos {


TEUCHOS_UNIT_TEST( GetBaseObjVoidPtr, polymorphicClasses )
{
  RCP<C> c_ptr(new C);
  RCP<B1> b1_ptr = c_ptr;
  RCP<B2> b2_ptr = c_ptr;
  RCP<A> a_ptr = c_ptr;
  TEST_EQUALITY( getBaseObjVoidPtr(&*c_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_EQUALITY( getBaseObjVoidPtr(&*b1_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_INEQUALITY( static_cast<const void*>(&*b1_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_EQUALITY( getBaseObjVoidPtr(&*b2_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_INEQUALITY( static_cast<const void*>(&*b2_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_EQUALITY( getBaseObjVoidPtr(&*a_ptr), static_cast<const void*>(&*c_ptr) );
  TEST_INEQUALITY( static_cast<const void*>(&*a_ptr), static_cast<const void*>(&*c_ptr) );
}


TEUCHOS_UNIT_TEST( GetBaseObjVoidPtr, nonPolymorphicClasses )
{
  RCP<E> e_ptr(new E);
  RCP<D> d_ptr = e_ptr;
  TEST_EQUALITY( getBaseObjVoidPtr(&*e_ptr), static_cast<const void*>(&*e_ptr) );
  TEST_EQUALITY( getBaseObjVoidPtr(&*d_ptr), static_cast<const void*>(&*e_ptr) );
}


TEUCHOS_UNIT_TEST( GetBaseObjVoidPtr, nonPolymorphicBuiltInTypes )
{
  RCP<int> i_ptr(new int);
  TEST_EQUALITY( getBaseObjVoidPtr(&*i_ptr), static_cast<const void*>(&*i_ptr) );
}


} // namespace Teuchos namespace
