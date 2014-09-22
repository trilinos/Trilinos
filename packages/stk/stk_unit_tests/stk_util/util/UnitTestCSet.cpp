// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>
#include <stk_util/util/CSet.hpp>       // for CSet


namespace stk {
namespace cset_unit {

class A {
public:
  virtual int id() const = 0 ;
  virtual ~A(){}
};

class B {
public:
  virtual int id() const = 0 ;
  virtual ~B(){}
};

class U : public A {
public:
  int id() const 
  {
    return static_cast<int>(ID);
  }
  ~U() {}
  enum {ID = 0};
};

class V : public B {
public:
  int id() const 
  {
    return static_cast<int>(ID);
  }
  ~V() {}
  enum {ID = 1};
};

class W : public B {
public:
  int id() const 
  {
    return static_cast<int>(ID);
  }
  ~W() {}
  enum {ID = 2};
};

class X : public A , public B {
public:
  int id() const
  {
    return static_cast<int>(ID);
  }
  ~X() {}
  enum {ID = 3};
};

class Y : public A , public B {
public:
  int id() const
  {
    return static_cast<int>(ID);
  }
  ~Y() {}
  enum {ID = 4};
};

class Z {
public:
  int id() const
  {
    return static_cast<int>(ID);
  }
  ~Z() {}
  enum {ID = 5};
};

}//namespace cset_unit
}//namespace stk

using namespace stk;
using namespace stk::cset_unit;

TEST( UnitTestCSet, UnitTest)
{
//This unit-test imported from its previous home in the bottom of
//the CSet implementation file.
  const A * sa = 0 ;
  const B * sb = 0 ;
  bool flag = false ;

  U  u;
  V  v;
  W  w;
  X  x;
  Y  y;

  {
    CSet cs ;

    sa = cs.insert_no_delete<A>(&u);
    ASSERT_TRUE(sa->id() == static_cast<int>(U::ID));

    sb = cs.insert_no_delete<B>(&v);
    ASSERT_TRUE(sb->id() == static_cast<int>(V::ID));

    // Should not replace:
    sb = cs.insert_no_delete<B>(&w);
    ASSERT_TRUE(sb->id() == static_cast<int>(V::ID));

    flag = cs.remove<A>( &u );
    ASSERT_TRUE(flag);

    flag = cs.remove<B>( &v );
    ASSERT_TRUE(flag);

    sa = cs.insert_no_delete<A>(&x);
    sb = cs.insert_no_delete<B>(&x);
    ASSERT_TRUE(sa->id() == static_cast<int>(X::ID));
    ASSERT_TRUE(sb->id() == static_cast<int>(X::ID));

    sa = cs.insert_no_delete<A>(&y);
    sb = cs.insert_no_delete<B>(&y);
    ASSERT_TRUE(sa->id() == static_cast<int>(X::ID));
    ASSERT_TRUE(sb->id() == static_cast<int>(X::ID));
  }
}

