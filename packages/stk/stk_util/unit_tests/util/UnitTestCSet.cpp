/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

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

