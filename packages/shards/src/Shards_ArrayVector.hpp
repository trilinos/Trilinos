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
    ArrayHelp<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
      help_type ;

  typedef
    Array<Scalar,Order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>
      BaseType ;

  void assign()
    {
      const typename BaseType::size_type n = BaseType::size();
      m_storage.resize( n );
      BaseType::m_ptr = n ? & m_storage[0] : NULL ;
    }

  ArrayVector();
  ArrayVector( const ArrayVector & );
  ArrayVector & operator = ( const ArrayVector & );

public:

  ~ArrayVector() {}

  void resize( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 ,
               const unsigned n5 , const unsigned n6 ,
               const unsigned n7 , const unsigned n8 )
    { help_type::assign(BaseType::m_stride,n1,n2,n3,n4,n5,n6,n7,n8); assign(); }

  void resize( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 ,
               const unsigned n5 , const unsigned n6 ,
               const unsigned n7 )
    { help_type::assign(BaseType::m_stride,n1,n2,n3,n4,n5,n6,n7); assign(); }

  void resize( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 ,
               const unsigned n5 , const unsigned n6 )
    { help_type::assign(BaseType::m_stride,n1,n2,n3,n4,n5,n6); assign(); }

  void resize( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 ,
               const unsigned n5 )
    { help_type::assign(BaseType::m_stride,n1,n2,n3,n4,n5); assign(); }

  void resize( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 )
    { help_type::assign(BaseType::m_stride,n1,n2,n3,n4); assign(); }

  void resize( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 )
    { help_type::assign(BaseType::m_stride,n1,n2,n3); assign(); }

  void resize( const unsigned n1 , const unsigned n2 )
    { help_type::assign(BaseType::m_stride,n1,n2); assign(); }

  void resize( const unsigned n1 )
    { help_type::assign(BaseType::m_stride,n1); assign(); }


  ArrayVector( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 ,
               const unsigned n5 , const unsigned n6 ,
               const unsigned n7 , const unsigned n8 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4,n5,n6,n7,n8); }

  ArrayVector( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 ,
               const unsigned n5 , const unsigned n6 ,
               const unsigned n7 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4,n5,n6,n7); }

  ArrayVector( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 ,
               const unsigned n5 , const unsigned n6 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4,n5,n6); }

  ArrayVector( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 ,
               const unsigned n5 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4,n5); }

  ArrayVector( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 , const unsigned n4 )
    : BaseType(), m_storage() { resize(n1,n2,n3,n4); }

  ArrayVector( const unsigned n1 , const unsigned n2 ,
               const unsigned n3 )
    : BaseType(), m_storage() { resize(n1,n2,n3); }

  ArrayVector( const unsigned n1 , const unsigned n2 )
    : BaseType(), m_storage() { resize(n1,n2); }

  ArrayVector( const unsigned n1 )
    : BaseType(), m_storage() { resize(n1); }

};

/** \} */

}

#endif /* Shards_ArrayVector_hpp */

