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

#ifndef KOKKOS_ARRAYVIEWRAWDATA_HPP
#define KOKKOS_ARRAYVIEWRAWDATA_HPP

#include <cstddef>

#define KOKKOS_MDARRAY_CHECK( EXPR ) /* */

namespace Kokkos {

template< typename ValueType , class DeviceMap >
class MDArrayView ;

/** \brief  Implementation contiguous multidimensional array view.
 *
 * \param  ValueType  The scalar type of the array values.
 * \param  LocalDevice  The local computational device in
 *                      which the array values reside.
 *
 *  Two distinct view semantics are supported.
 *
 *  Multiple views to the same allocated memory are tracked
 *  in a ring such that new views are added to the ring,
 *  destroyed views are removed from the ring, and
 *  when the last view is destroyed the LocalDevice
 *  is called to destroy the allocated memory.
 *
 *  Views to non-allocated memory are NOT tracked and
 *  destruction of such views has no effect on the
 *  viewed memory.
 */
template< typename ValueType , class LocalDevice >
class MDArrayViewRawData {
private:

  template< typename , class > friend class MDArrayView ;

  typedef LocalDevice                           local_device_type ;
  typedef typename local_device_type::size_type size_type ;
  typedef ValueType                             value_type ;

  /*------------------------------------------------------------------*/

  /** \brief  Pointer to next view in the ring, if exists */
  MDArrayViewRawData * m_next_on_host ;

  /** \brief  Pointer to memory on the device */
  value_type         * m_ptr_on_device ;

  /** \brief  Rank of the multidimensional array view */
  size_type            m_rank ;

  /** \brief  Dimensions of the multidimensional array view */
  size_type            m_dimension[8] ;

  /*------------------------------------------------------------------*/
  /** \brief  Rank and multi-index bounds checking */

  void require_in_bounds( size_type i0 , size_type i1 ,
                          size_type i2 , size_type i3 ,
                          size_type i4 , size_type i5 ,
                          size_type i6 , size_type i7 );

  void require_in_bounds( size_type i0 , size_type i1 ,
                          size_type i2 , size_type i3 ,
                          size_type i4 , size_type i5 ,
                          size_type i6 );

  void require_in_bounds( size_type i0 , size_type i1 ,
                          size_type i2 , size_type i3 ,
                          size_type i4 , size_type i5 );

  void require_in_bounds( size_type i0 , size_type i1 ,
                          size_type i2 , size_type i3 ,
                          size_type i4 );

  void require_in_bounds( size_type i0 , size_type i1 ,
                          size_type i2 , size_type i3 );

  void require_in_bounds( size_type i0 , size_type i1 ,
                          size_type i2 );

  void require_in_bounds( size_type i0 , size_type i1 );

  void require_in_bounds( size_type i0 );

  /*------------------------------------------------------------------*/

  static value_type * allocate_on_device( size_type count , const std::string & label )
    { return (value_type *) local_device_type::allocate_memory( sizeof(ValueType) , count , label ); }

  static void deallocate_on_device( value_type * ptr_on_device )
    {  local_device_type::deallocate_memory( ptr_on_device ); }

  /*------------------------------------------------------------------*/
  /** \brief  Another view to the same memory. */
  inline
  void insert_view( const MDArrayViewRawData & rhs )
    {
      m_next_on_host  = NULL ;
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

      // If viewing allocated memory (exists and in a ring)
      // then insert this view into the ring of views.
      if ( NULL != rhs.m_ptr_on_device && NULL != rhs.m_next_on_host ) {
        MDArrayViewRawData & another_view = * rhs.m_next_on_host ;
        m_next_on_host = another_view.m_next_on_host ;
        another_view.m_next_on_host = this ;
      }
    }

  /** \brief  Remove this view of the ring of equivalent views.
   *          If the last view the request the map to destroy the array.
   */
  inline
  void clear_view()
    {
      if ( m_next_on_host == this ) {
        // This is the last view to allocated memory
        // so destroy the memory.
        deallocate_on_device( m_ptr_on_device );
      }
      else if ( m_next_on_host != NULL ) {
        // Remove this view from the ring of views
        // to the given allocated memory.
        MDArrayViewRawData * prev = m_next_on_host ;
        for ( ; this != prev->m_next_on_host ; prev = prev->m_next_on_host );
        prev->m_next_on_host = m_next_on_host ;
      }

      m_next_on_host  = NULL ;
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
  /** \brief  Construct a NULL view */
  inline
  MDArrayViewRawData()
    : m_next_on_host( NULL )
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

  /** \brief  Construct another view of the array */
  inline
  MDArrayViewRawData( const MDArrayViewRawData & rhs )
    { insert_view( rhs ); }

  /** \brief  Assign a view of the array */
  inline
  MDArrayViewRawData & operator = ( const MDArrayViewRawData & rhs )
    { clear_view(); insert_view( rhs ); return *this ; }

  /**  \brief  Destroy this view of the array */
  inline
  ~MDArrayViewRawData() { clear_view(); }

  /*------------------------------------------------------------------*/

  MDArrayViewRawData( size_type n0 , size_type n1 ,
                      size_type n2 , size_type n3 ,
                      size_type n4 , size_type n5 ,
                      size_type n6 , size_type n7 ,
                      const std::string & label )
    : m_next_on_host( this )
    , m_ptr_on_device( allocate_on_device( n0*n1*n2*n3*n4*n5*n6*n7 , label ) )
    , m_rank( 8 )
    {
      m_dimension[0] = n0 ;
      m_dimension[1] = n1 ;
      m_dimension[2] = n2 ;
      m_dimension[3] = n3 ;
      m_dimension[4] = n4 ;
      m_dimension[5] = n5 ;
      m_dimension[6] = n6 ;
      m_dimension[7] = n7 ;
    }

  MDArrayViewRawData( size_type n0 , size_type n1 ,
                      size_type n2 , size_type n3 ,
                      size_type n4 , size_type n5 ,
                      size_type n6 ,
                      const std::string & label )
    : m_next_on_host( this )
    , m_ptr_on_device( allocate_on_device( n0*n1*n2*n3*n4*n5*n6 , label ) )
    , m_rank( 7 )
    {
      m_dimension[0] = n0 ;
      m_dimension[1] = n1 ;
      m_dimension[2] = n2 ;
      m_dimension[3] = n3 ;
      m_dimension[4] = n4 ;
      m_dimension[5] = n5 ;
      m_dimension[6] = n6 ;
      m_dimension[7] = 0 ;
    }

  MDArrayViewRawData( size_type n0 , size_type n1 ,
                      size_type n2 , size_type n3 ,
                      size_type n4 , size_type n5 ,
                      const std::string & label )
    : m_next_on_host( this )
    , m_ptr_on_device( allocate_on_device( n0*n1*n2*n3*n4*n5 , label ) )
    , m_rank( 6 )
    {
      m_dimension[0] = n0 ;
      m_dimension[1] = n1 ;
      m_dimension[2] = n2 ;
      m_dimension[3] = n3 ;
      m_dimension[4] = n4 ;
      m_dimension[5] = n5 ;
      m_dimension[6] = 0 ;
      m_dimension[7] = 0 ;
    }

  MDArrayViewRawData( size_type n0 , size_type n1 ,
                      size_type n2 , size_type n3 ,
                      size_type n4 ,
                      const std::string & label )
    : m_next_on_host( this )
    , m_ptr_on_device( allocate_on_device( n0*n1*n2*n3*n4 , label ) )
    , m_rank( 5 )
    {
      m_dimension[0] = n0 ;
      m_dimension[1] = n1 ;
      m_dimension[2] = n2 ;
      m_dimension[3] = n3 ;
      m_dimension[4] = n4 ;
      m_dimension[5] = 0 ;
      m_dimension[6] = 0 ;
      m_dimension[7] = 0 ;
    }

  MDArrayViewRawData( size_type n0 , size_type n1 ,
                      size_type n2 , size_type n3 ,
                      const std::string & label )
    : m_next_on_host( this )
    , m_ptr_on_device( allocate_on_device( n0*n1*n2*n3 , label ) )
    , m_rank( 4 )
    {
      m_dimension[0] = n0 ;
      m_dimension[1] = n1 ;
      m_dimension[2] = n2 ;
      m_dimension[3] = n3 ;
      m_dimension[4] = 0 ;
      m_dimension[5] = 0 ;
      m_dimension[6] = 0 ;
      m_dimension[7] = 0 ;
    }

  MDArrayViewRawData( size_type n0 , size_type n1 ,
                      size_type n2 ,
                      const std::string & label )
    : m_next_on_host( this )
    , m_ptr_on_device( allocate_on_device( n0*n1*n2 , label ) )
    , m_rank( 3 )
    {
      m_dimension[0] = n0 ;
      m_dimension[1] = n1 ;
      m_dimension[2] = n2 ;
      m_dimension[3] = 0 ;
      m_dimension[4] = 0 ;
      m_dimension[5] = 0 ;
      m_dimension[6] = 0 ;
      m_dimension[7] = 0 ;
    }

  MDArrayViewRawData( size_type n0 , size_type n1 ,
                      const std::string & label )
    : m_next_on_host( this )
    , m_ptr_on_device( allocate_on_device( n0*n1 , label ) )
    , m_rank( 2 )
    {
      m_dimension[0] = n0 ;
      m_dimension[1] = n1 ;
      m_dimension[2] = 0 ;
      m_dimension[3] = 0 ;
      m_dimension[4] = 0 ;
      m_dimension[5] = 0 ;
      m_dimension[6] = 0 ;
      m_dimension[7] = 0 ;
    }

  MDArrayViewRawData( size_type n0 ,
                      const std::string & label )
    : m_next_on_host( this )
    , m_ptr_on_device( allocate_on_device( n0 , label ) )
    , m_rank( 1 )
    {
      m_dimension[0] = n0 ;
      m_dimension[1] = 0 ;
      m_dimension[2] = 0 ;
      m_dimension[3] = 0 ;
      m_dimension[4] = 0 ;
      m_dimension[5] = 0 ;
      m_dimension[6] = 0 ;
      m_dimension[7] = 0 ;
    }
};

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

} // namespace Kokkos

#endif /* KOKKOS_ARRAYVIEWRAWDATA_HPP */

