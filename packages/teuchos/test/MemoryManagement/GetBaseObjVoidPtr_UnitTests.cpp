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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
