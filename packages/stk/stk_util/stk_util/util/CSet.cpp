#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

#include <stk_util/util/CSet.hpp>

namespace stk {

namespace {

typedef void (* DeleteFunction )( void * );

typedef std::pair< const std::type_info * , DeleteFunction > Manager ;

// Comparison for sorted vector

struct less_cset {
  bool operator()( const Manager        & lhs ,
                   const std::type_info & rhs ) const ;
  bool operator()( const std::type_info & lhs ,
                   const Manager        & rhs ) const ;
};

// On some systems, namely AIX, std::type_info::before(...)
// has a bug where it returns true instead of false for equality.
// Thus we pay a small price on all systems to specifically
// test for and eliminate equality.

bool less_cset::operator()( const Manager        & lhs ,
                            const std::type_info & rhs ) const
{ return lhs.first->before( rhs ) && * lhs.first != rhs ; }

bool less_cset::operator()( const std::type_info & lhs ,
                            const Manager        & rhs ) const
{ return lhs.before( *rhs.first ) && lhs != *rhs.first ; }


std::vector< Manager >::const_iterator
lower_bound( const std::vector< Manager > & v , const std::type_info & t )
{
  std::vector< Manager >::const_iterator i = v.begin();
  std::vector< Manager >::const_iterator j = v.end();

  return std::lower_bound( i , j , t , less_cset() );
}

std::vector< Manager >::iterator
lower_bound( std::vector< Manager > & v , const std::type_info & t )
{
  std::vector< Manager >::iterator i = v.begin();
  std::vector< Manager >::iterator j = v.end();

  return std::lower_bound( i , j , t , less_cset() );
}

}

//----------------------------------------------------------------------

const void * CSet::p_get( const std::type_info & t ) const
{
  const void * result = NULL ;

  const std::vector< Manager >::const_iterator im = lower_bound(m_manager,t);

  if ( im < m_manager.end() && t == * im->first ) {
    const size_t offset = im - m_manager.begin();
    result = m_value[ offset ];
  }

  return result ;
}

const void *
CSet::p_insert( const Manager & m , const void * v )
{
  std::vector< Manager >::iterator im = lower_bound( m_manager , * m.first );

  const size_t offset = im - m_manager.begin();

  std::vector<const void *>::iterator iv = m_value.begin();
  std::advance( iv , offset );

  if ( im == m_manager.end() || * m.first != * im->first ) {
    im = m_manager.insert( im , m );
    iv = m_value  .insert( iv , v );
  }

  return *iv ;
}

bool CSet::p_remove( const std::type_info & t , const void * v )
{
  const std::vector< Manager >::iterator im = lower_bound( m_manager , t );

  const size_t offset = im - m_manager.begin();

  std::vector<const void *>::iterator iv = m_value.begin();
  std::advance( iv , offset );

  const bool result = im != m_manager.end() && t == * im->first && v == * iv ;

  if ( result ) {
    m_manager.erase( im );
    m_value  .erase( iv );
  }

  return result ;
}

//----------------------------------------------------------------------

CSet::~CSet()
{
  try {
    const size_t n = m_manager.size();
    for ( size_t i = 0 ; i < n ; ++i ) {
      try {
        if ( m_manager[i].second ) {
          (*m_manager[i].second)( const_cast<void*>( m_value[i] ) );
        }
      } catch(...) {}
    }
  } catch(...) {}
}

CSet::CSet() : m_manager(), m_value() {}

} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------


#ifdef UNIT_TEST

namespace stk {
namespace unit_test {

class A {
public:
  virtual const char * name() const = 0 ;
  virtual ~A();
};

class B {
public:
  virtual const char * name() const = 0 ;
  virtual ~B();
};

void DoNotDelete( A * a )
{ std::cout << "DoNotDelete(" << a->name() << ")" << std::endl ; }

void DoDelete( A * a )
{
  std::cout << "DoDelete(" << a->name() << ")" << std::endl ;
  delete a ;
}

void DoNotDelete( B * b )
{ std::cout << "DoNotDelete(" << b->name() << ")" << std::endl ; }

void DoDelete( B * b )
{
  std::cout << "DoDelete(" << b->name() << ")" << std::endl ;
  delete b ;
}

class U : public A {
public:
  const char * name() const ;
  ~U() {}
};

class V : public B {
public:
  const char * name() const ;
  ~V() {}
};

class W : public B {
public:
  const char * name() const ;
  ~W() {}
};

class X : public A , public B {
public:
  const char * name() const ;
  ~X() {}
};

class Y : public A , public B {
public:
  const char * name() const ;
  ~Y() {}
};

class Z {
public:
  const char * name() const ;
  ~Z() {}
};

//----------------------------------------------------------------------

int cset()
{
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

    sa = cs.insert<A>(u,true);
    std::cout << "cs.insert<A>(u,true)->name() = " << sa->name() << std::endl ;

    sb = cs.insert<B>(v,true);
    std::cout << "cs.insert<B>(v,true)->name() = " << sb->name() << std::endl ;

    // Should not replace:
    sb = cs.insert<B>(w,true);
    std::cout << "cs.insert<B>(w,true)->name() = " << sb->name() << std::endl ;

    flag = cs.remove<A>( u );
    std::cout << "s.remove<A>(u) = " << flag << std::endl ;

    flag = cs.remove<B>( v );
    std::cout << "s.remove<B>(v) = " << flag << std::endl ;

    sa = cs.insert<A>(x);
    sb = cs.insert<B>(x);
    std::cout << "s.insert<A>(x)->name() = " << sa->name() << std::endl ;
    std::cout << "s.insert<B>(x)->name() = " << sb->name() << std::endl ;

    sa = cs.insert<A>(y);
    sb = cs.insert<B>(y);
    std::cout << "s.insert<A>(y)->name() = " << sa->name() << std::endl ;
    std::cout << "s.insert<B>(y)->name() = " << sb->name() << std::endl ;
  }

  delete x ; x = NULL ;
  delete y ; y = NULL ;
  delete w ; w = NULL ;
  delete v ; v = NULL ;
  delete u ; u = NULL ;

  return 0 ;
}

//----------------------------------------------------------------------

A::~A() {}
B::~B() {}

const char * U::name() const
{
  static const char n[] = "U" ;
  return n ;
}

const char * V::name() const
{
  static const char n[] = "V" ;
  return n ;
}

const char * W::name() const
{
  static const char n[] = "W" ;
  return n ;
}

const char * X::name() const
{
  static const char n[] = "X" ;
  return n ;
}

const char * Y::name() const
{
  static const char n[] = "Y" ;
  return n ;
}

const char * Z::name() const
{
  static const char n[] = "Z" ;
  return n ;
}

}
}

int main()
{
  return phdmesh::unit_test::cset();
}

#endif



