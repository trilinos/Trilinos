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

#include "TestClasses.hpp"


int C::A_g_on_delete_ = -2;


void deallocA(A* ptr)
{
  std::cout << "\nCalled deallocA(...)!\n";
  delete ptr;
}


void deallocHandleA(A** handle)
{
  std::cout << "\nCalled deallocHandleA(...)!\n";
  A *ptr = *handle;
  delete ptr;
  *handle = 0;
}


struct UndefinedType {
  int val_;
};


Opaque_handle createOpaque()
{
  Opaque_handle obj = new UndefinedType;
  obj->val_ = getOpaqueValue_return;
  return obj;
}


int getOpaqueValue( Opaque_handle opaque )
{
  return opaque->val_;
}


void destroyOpaque( Opaque_handle * opaque )
{
  std::cout << "\nCalled destroyOpaque(...)!\n";
  delete *opaque;
  *opaque = 0;
}


struct UndefinedType2 {
  int val_;
};


Opaque2_handle createOpaque2()
{
  Opaque2_handle obj = new UndefinedType2;
  obj->val_ = getOpaque2Value_return;
  return obj;
}


int getOpaque2Value( Opaque2_handle opaque )
{
  return opaque->val_;
}


void destroyOpaque2( Opaque2_handle * opaque )
{
  std::cout << "\nCalled destroyOpaque2(...)!\n";
  delete *opaque;
  *opaque = 0;
}
