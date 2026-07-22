// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
