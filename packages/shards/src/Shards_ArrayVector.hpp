/*------------------------------------------------------------------------*/
/*               shards : Shared Discretization Tools                     */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*                                                                        */
/* Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)             */
/*------------------------------------------------------------------------*/

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
    : m_storage()
  {}

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

