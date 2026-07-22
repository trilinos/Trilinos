// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Shards_ArrayVector_hpp
#define Shards_ArrayVector_hpp

//----------------------------------------------------------------------

#include <Shards_Array.hpp>

//----------------------------------------------------------------------

namespace shards {

/** \addtogroup  shards_package_array
 *  \{
 */

template< typename Scalar , ArrayOrder Order ,
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class ArrayVector ;

//----------------------------------------------------------------------

template< typename Scalar , ArrayOrder Order >
class ArrayVector<Scalar,Order,void,void,void,void, void,void,void,void>
  : public Array<Scalar,Order,void,void,void,void, void,void,void,void>
{
private:
  std::vector<Scalar> m_storage ;

  typedef
    Array<Scalar,Order,void,void,void,void, void,void,void,void>
      BaseType ;

  ArrayVector( const ArrayVector & );
  ArrayVector & operator = ( const ArrayVector & );

  Scalar * get_ptr() { return m_storage.empty() ? NULL : & m_storage[0] ; }

public:

  typedef array_traits::int_t size_type ;

  ArrayVector()
    : BaseType() , m_storage() {}

  ~ArrayVector() {}

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 ,
               const size_type n7 , const size_type n8 )
    {
      m_storage.resize( n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 );
      BaseType::template assign<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
        ( get_ptr() ,n1,n2,n3,n4,n5,n6,n7,n8);
    }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 >
  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 ,
               const size_type n7 )
    {
      m_storage.resize( n1 * n2 * n3 * n4 * n5 * n6 * n7 );
      BaseType::template assign<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7>
        ( get_ptr() ,n1,n2,n3,n4,n5,n6,n7);
    }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 >
  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 )
    {
      m_storage.resize( n1 * n2 * n3 * n4 * n5 * n6 );
      BaseType::template assign<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6>
        ( get_ptr() ,n1,n2,n3,n4,n5,n6);
    }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 >
  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 )
    {
      m_storage.resize( n1 * n2 * n3 * n4 * n5 );
      BaseType::template assign<Tag1,Tag2,Tag3,Tag4,Tag5>
       ( get_ptr(),n1,n2,n3,n4,n5);
    }

  template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 )
    {
      m_storage.resize( n1 * n2 * n3 * n4 );
      BaseType::template assign<Tag1,Tag2,Tag3,Tag4>
        ( get_ptr(),n1,n2,n3,n4);
    }

  template< class Tag1 , class Tag2 , class Tag3 >
  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 )
    {
      m_storage.resize( n1 * n2 * n3 );
      BaseType::template assign<Tag1,Tag2,Tag3>
        ( get_ptr(),n1,n2,n3);
    }

  template< class Tag1 , class Tag2 >
  void resize( const size_type n1 , const size_type n2 )
    {
      m_storage.resize( n1 * n2 );
      BaseType::template assign<Tag1,Tag2>( get_ptr(),n1,n2);
    }

  template< class Tag1 >
  void resize( const size_type n1 )
    {
      m_storage.resize( n1 );
      BaseType::template assign<Tag1>( get_ptr(),n1);
    }
};

//----------------------------------------------------------------------

template< typename Scalar , ArrayOrder Order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class ArrayVector
  : public Array<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
{
private:
  std::vector<Scalar> m_storage ;

  typedef
    array_traits::Helper<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
      help_type ;

  typedef
    Array<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
      BaseType ;

  ArrayVector( const ArrayVector & );
  ArrayVector & operator = ( const ArrayVector & );

  Scalar * get_ptr() { return m_storage.empty() ? NULL : & m_storage[0] ; }

public:

  typedef array_traits::int_t size_type ;

  ArrayVector()
    : BaseType() , m_storage() {}

  ~ArrayVector() {}

  void resize( const size_type * const dims )
    {
      help_type::assign( BaseType::m_stride , dims );
      const typename BaseType::size_type n = BaseType::size();
      m_storage.resize( n );
      BaseType::m_ptr = get_ptr();
    }

  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 ,
               const size_type n7 , const size_type n8 )
    {
      array_traits::CheckRank<8,BaseType::Rank>::ok();
      m_storage.resize( n1 * n2 * n3 * n4 * n5 * n6 * n7 * n8 );
      BaseType::assign( get_ptr() ,n1,n2,n3,n4,n5,n6,n7,n8);
    }

  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 ,
               const size_type n7 )
    {
      array_traits::CheckRank<7,BaseType::Rank>::ok();
      m_storage.resize( n1 * n2 * n3 * n4 * n5 * n6 * n7 );
      BaseType::assign( get_ptr(),n1,n2,n3,n4,n5,n6,n7);
    }

  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 )
    {
      array_traits::CheckRank<6,BaseType::Rank>::ok();
      m_storage.resize( n1 * n2 * n3 * n4 * n5 * n6 );
      BaseType::assign( get_ptr(),n1,n2,n3,n4,n5,n6);
    }

  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 )
    {
      array_traits::CheckRank<5,BaseType::Rank>::ok();
      m_storage.resize( n1 * n2 * n3 * n4 * n5 );
      BaseType::assign( get_ptr(),n1,n2,n3,n4,n5);
    }

  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 )
    {
      array_traits::CheckRank<4,BaseType::Rank>::ok();
      m_storage.resize( n1 * n2 * n3 * n4 );
      BaseType::assign( get_ptr(),n1,n2,n3,n4);
    }

  void resize( const size_type n1 , const size_type n2 ,
               const size_type n3 )
    {
      array_traits::CheckRank<3,BaseType::Rank>::ok();
      m_storage.resize( n1 * n2 * n3 );
      BaseType::assign( get_ptr(),n1,n2,n3);
    }

  void resize( const size_type n1 , const size_type n2 )
    {
      array_traits::CheckRank<2,BaseType::Rank>::ok();
      m_storage.resize( n1 * n2 );
      BaseType::assign( get_ptr(),n1,n2);
    }

  void resize( const size_type n1 )
    {
      array_traits::CheckRank<1,BaseType::Rank>::ok();
      m_storage.resize( n1 );
      BaseType::assign( get_ptr(),n1);
    }


  ArrayVector( const size_type * const dims )
    : BaseType(), m_storage() { resize( dims ); }

  ArrayVector( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 ,
               const size_type n7 , const size_type n8 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4,n5,n6,n7,n8); }

  ArrayVector( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 ,
               const size_type n7 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4,n5,n6,n7); }

  ArrayVector( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 , const size_type n6 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4,n5,n6); }

  ArrayVector( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 ,
               const size_type n5 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4,n5); }

  ArrayVector( const size_type n1 , const size_type n2 ,
               const size_type n3 , const size_type n4 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4); }

  ArrayVector( const size_type n1 , const size_type n2 ,
               const size_type n3 )
    : BaseType(), m_storage() { resize(n1,n2,n3); }

  ArrayVector( const size_type n1 , const size_type n2 )
    : BaseType(), m_storage() { resize(n1,n2); }

  ArrayVector( const size_type n1 )
    : BaseType(), m_storage() { resize(n1); }

};

/** \} */

}

#endif /* Shards_ArrayVector_hpp */

