/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_BASEMAPPEDARRAY_HPP
#define KOKKOS_BASEMAPPEDARRAY_HPP

#include <cstddef>
#include <list>

namespace Kokkos {

class BaseDeviceMemory ;
class BaseMappedArrayRepository ;

class BaseMappedArray {
public:

  typedef size_t size_type ;

  /*------------------------------------------------------------------*/

  inline
  size_type rank() const { return m_rank ; }

  inline
  const size_type * dimension() const { return m_dimension ; }

  template < typename iType >
  inline
  size_type dimension( const iType & ordinal ) const
    {
      KOKKOS_MDARRAY_CHECK( require_ordinal_in_bounds(m_rank,ordinal) );
      return m_dimension[ordinal] ;
    }

  /*------------------------------------------------------------------*/

  template< typename ValueType >
  inline
  ValueType * pointer_on_device()
    { return (ValueType *) m_ptr_on_device ; }

  /*------------------------------------------------------------------*/

  void require_multi_index_in_bounds( ptrdiff_t i0 , ptrdiff_t i1 ,
                                      ptrdiff_t i2 , ptrdiff_t i3 ,
                                      ptrdiff_t i4 , ptrdiff_t i5 ,
                                      ptrdiff_t i6 , ptrdiff_t iP );

  void require_multi_index_in_bounds( ptrdiff_t i0 , ptrdiff_t i1 ,
                                      ptrdiff_t i2 , ptrdiff_t i3 ,
                                      ptrdiff_t i4 , ptrdiff_t i5 ,
                                      ptrdiff_t iP );

  /*------------------------------------------------------------------*/

  BaseMappedArray()
    : m_map_on_host( NULL )
    , m_next_on_host( this )
    , m_ptr_on_device( NULL )
    , m_rank( 0 ) 
    {
      m_dimension[0] = 0 ;
      m_dimension[1] = 0 ;
      m_dimension[2] = 0 ;
      m_dimension[3] = 0 ;
      m_dimension[4] = 0 ;
      m_dimension[5] = 0 ;
      m_dimension[6] = 0 ;
      m_dimension[7] = 0 ;
    }

  /** \brief  Construct a view of the array */
  BaseMappedArray( const BaseMappedArray & rhs )
    : m_map_on_host( NULL )
    , m_next_on_host(  rhs.m_next_on_host )
    , m_ptr_on_device( rhs.m_ptr_on_device )
    , m_rank(          rhs.m_rank )
    {
      m_dimension[0] = rhs.m_dimension[0] ;
      m_dimension[1] = rhs.m_dimension[1] ;
      m_dimension[2] = rhs.m_dimension[2] ;
      m_dimension[3] = rhs.m_dimension[3] ;
      m_dimension[4] = rhs.m_dimension[4] ;
      m_dimension[5] = rhs.m_dimension[5] ;
      m_dimension[6] = rhs.m_dimension[6] ;
      m_dimension[7] = rhs.m_dimension[7] ;

      const_cast<BaseMappedArray&>(rhs).m_next_on_host = this ;
    }

  /** \brief  Assign a view of the array */
  BaseMappedArray & operator = ( const BaseMappedArray & rhs )
    {
      prev_on_host()->m_next_on_host = m_next_on_host ;

      m_next_on_host  = rhs.m_next_on_host ;
      m_ptr_on_device = rhs.m_ptr_on_device ;
      m_rank          = rhs.m_rank ;
      m_dimension[0]  = rhs.m_dimension[0] ;
      m_dimension[1]  = rhs.m_dimension[1] ;
      m_dimension[2]  = rhs.m_dimension[2] ;
      m_dimension[3]  = rhs.m_dimension[3] ;
      m_dimension[4]  = rhs.m_dimension[4] ;
      m_dimension[5]  = rhs.m_dimension[5] ;
      m_dimension[6]  = rhs.m_dimension[6] ;
      m_dimension[7]  = rhs.m_dimension[7] ;

      const_cast<BaseMappedArray&>(rhs).m_next_on_host = this ;

      return *this ;
    }

  /**  Destroy this view of the array */
  inline
  ~BaseMappedArray()
    {
      prev_on_host()->m_next_on_host = m_next_on_host ;

      m_map_on_host   = NULL ;
      m_next_on_host  = this ;
      m_ptr_on_device = NULL ;
      m_rank          = 0 ;
      m_dimension[0]  = 0 ;
      m_dimension[1]  = 0 ;
      m_dimension[2]  = 0 ;
      m_dimension[3]  = 0 ;
      m_dimension[4]  = 0 ;
      m_dimension[5]  = 0 ;
      m_dimension[6]  = 0 ;
      m_dimension[7]  = 0 ;
    }

  /*------------------------------------------------------------------*/

  BaseMappedArrayRepository * repository() const ;

private:

  friend class BaseMappedArrayRepository ;

  void clear();

  void assign( BaseMappedArrayRepository * const map ,
               void * const ptr_on_device ,
               const size_type rank , const size_type dimension[] );

  BaseMappedArray * owned_view();

  void require_not_owned_view() const ;
  void require_next_exists() const ;

  BaseMappedArray * next_on_host() const ;
  BaseMappedArray * prev_on_host() const ;


  BaseMappedArrayRepository * m_map_on_host ;
  BaseMappedArray           * m_next_on_host ;
  void                      * m_ptr_on_device ;
  size_type                   m_rank ;
  size_type                   m_dimension[8] ;
};

/*------------------------------------------------------------------------*/

class BaseMappedArrayRepository {
public:
  typedef BaseMappedArray::size_type size_type ;

  /*------------------------------------------------------------------------*/
  /** \brief  Query the parallel work count for this mapping */
  size_type parallel_work_count() const
    { return m_parallel_work_count ; }

  /*------------------------------------------------------------------------*/
  /** \brief  Deallocate a mapped array owned by this host map object.
   *
   *  A warning will be issued if there exists a view to the array
   *  other than the input view.
   */
  void deallocate( BaseMappedArray & );

  /*------------------------------------------------------------------------*/
  /** \brief  Allocate a mapped array owned by this host map object.
   *
   *  Precondition: The array cannot already be allocated.
   */
  void allocate( BaseMappedArray & array , size_type sizeof_value ,
                         size_type rank , const size_type dimension[] );

  /*------------------------------------------------------------------------*/

  BaseMappedArrayRepository( BaseDeviceMemory & ,
                             size_type parallel_work_count );

  ~BaseMappedArrayRepository();

private:

  /** \brief  Require that the array is not owned by any host map */
  static void require_array_is_not_owned( BaseMappedArray & );

  /** \brief  Require that the array is owned by this host map
   *          and return the owned member from 'm_allocated_arrays'.
   */
  std::list< BaseMappedArray >::iterator
    require_array_is_member( BaseMappedArray & );

  BaseMappedArrayRepository();
  BaseMappedArrayRepository( const BaseMappedArrayRepository & );
  BaseMappedArrayRepository & operator = ( const BaseMappedArrayRepository & );

  BaseDeviceMemory         &  m_base_device ;
  std::list<BaseMappedArray>  m_allocated_arrays ;
  size_type                   m_parallel_work_count ;
};

/*--------------------------------------------------------------------------*/

class BaseDeviceMemory {
public:
  typedef size_t size_type ;

  virtual void * allocate( size_type value_size ,
                           size_type chunk_count ,
                           size_type work_count ) = 0 ;

  virtual void deallocate( void * ) = 0 ;

  virtual ~BaseDeviceMemory() {}
};

/*--------------------------------------------------------------------------*/

} // namespace Kokkos


/*------------------------------------------------------------------------*/
/*

class CudaMap {
public:

  template< typename Scalar ,
            typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iType0 & i0 ,
                      const iType1 & i1 ,
                      const iType2 & i2 ,
                      const iType3 & i3 ,
                      const iType4 & i4 ,
                      const iType5 & i5 ,
                      const iType6 & i6 )
                      const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK(
        array.require_multi_index_in_bounds( i0,i1,i2,i3,i4,i5,i6,iP) );

      return ((value_type*)array.m_ptr)[ ( iP + array.m_dimension[7] *
                                         ( i0 + array.m_dimension[0] * 
                                         ( i1 + array.m_dimension[1] *
                                         ( i2 + array.m_dimension[2] *
                                         ( i3 + array.m_dimension[3] *
                                         ( i4 + array.m_dimension[4] *
                                         ( i5 + array.m_dimension[5] *
                                         ( i6 ))))))));
    }

  template< typename Array , typename iType0 , ... >
  void allocate( BaseMappedArray & , size_type nB , n0, n1, n2, n3, n4, n5 , n6 );

  void deallocate( MappedArrayBase & );
};  

} // namespace Kokkos
*/

#endif /* KOKKOS_MAPPEDARRAY_HPP */

