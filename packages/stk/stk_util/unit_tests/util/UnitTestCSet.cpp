/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/util/CSet.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

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
    return (int)ID;
  }
  ~U() {}
  enum {ID = 0};
};

class V : public B {
public:
  int id() const 
  {
    return (int)ID;
  }
  ~V() {}
  enum {ID = 1};
};

class W : public B {
public:
  int id() const 
  {
    return (int)ID;
  }
  ~W() {}
  enum {ID = 2};
};

class X : public A , public B {
public:
  int id() const
  {
    return (int)ID;
  }
  ~X() {}
  enum {ID = 3};
};

class Y : public A , public B {
public:
  int id() const
  {
    return (int)ID;
  }
  ~Y() {}
  enum {ID = 4};
};

class Z {
public:
  int id() const
  {
    return (int)ID;
  }
  ~Z() {}
  enum {ID = 5};
};

}//namespace cset_unit
}//namespace stk

using namespace stk;
using namespace stk::cset_unit;

STKUNIT_UNIT_TEST( UnitTestCSet, UnitTest)
{
//This unit-test imported from its previous home in the bottom of
//the CSet implementation file.
  const A * sa ;
  const B * sb ;
  bool flag ;

  U * u = new U();
  V * v = new V();
  W * w = new W();
  X * x = new X();
  Y * y = new Y();

  {
    CSet cs ;

    sa = cs.insert_no_delete<A>(u);
    STKUNIT_ASSERT(sa->id() == (int)U::ID);

    sb = cs.insert_no_delete<B>(v);
    STKUNIT_ASSERT(sb->id() == (int)V::ID);

    // Should not replace:
    sb = cs.insert_no_delete<B>(w);
    STKUNIT_ASSERT(sb->id() == (int)V::ID);

    flag = cs.remove<A>( u );
    STKUNIT_ASSERT(flag);

    flag = cs.remove<B>( v );
    STKUNIT_ASSERT(flag);

    sa = cs.insert_no_delete<A>(x);
    sb = cs.insert_no_delete<B>(x);
    STKUNIT_ASSERT(sa->id() == (int)X::ID);
    STKUNIT_ASSERT(sb->id() == (int)X::ID);

    sa = cs.insert_no_delete<A>(y);
    sb = cs.insert_no_delete<B>(y);
    STKUNIT_ASSERT(sa->id() == (int)X::ID);
    STKUNIT_ASSERT(sb->id() == (int)X::ID);
  }

  delete x ; x = NULL ;
  delete y ; y = NULL ;
  delete w ; w = NULL ;
  delete v ; v = NULL ;
  delete u ; u = NULL ;
}

