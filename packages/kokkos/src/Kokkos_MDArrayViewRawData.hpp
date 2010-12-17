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

#define KOKKOS_MDARRAY_CHECK( EXPR ) /* */

namespace Kokkos {

/** \brief  Implementation of data for a contiguous multidimensional array.
 *
 */
template< typename ValueType , class LocalDevice >
class MDArrayViewRawData {
public:
  typedef LocalDevice                           local_device_type ;
  typedef typename local_device_type::size_type size_type ;
  typedef ValueType                             value_type ;

  /*------------------------------------------------------------------*/

  MDArrayViewRawData * m_next_on_host ;
  value_type         * m_ptr_on_device ;
  size_type            m_rank ;
  size_type            m_dimension[8] ;

  /*------------------------------------------------------------------*/

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
  /** \brief  Insert this view into a ring of equivalent views */
  inline
  void insert_view( const MDArrayViewRawData & rhs )
    {
      m_rank         = rhs.m_rank ;
      m_dimension[0] = rhs.m_dimension[0] ;
      m_dimension[1] = rhs.m_dimension[1] ;
      m_dimension[2] = rhs.m_dimension[2] ;
      m_dimension[3] = rhs.m_dimension[3] ;
      m_dimension[4] = rhs.m_dimension[4] ;
      m_dimension[5] = rhs.m_dimension[5] ;
      m_dimension[6] = rhs.m_dimension[6] ;
      m_dimension[7] = rhs.m_dimension[7] ;

      if ( NULL != ( m_ptr_on_device = rhs.m_ptr_on_device ) ) {
        // Insert this new view into the ring of equivalent views.
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
      if ( m_ptr_on_device != NULL ) {
        if ( m_next_on_host != this ) {
          // Take this view out of the ring of views.
          MDArrayViewRawData * prev = m_next_on_host ;
          for ( ; this != prev->m_next_on_host ; prev = prev->m_next_on_host );
          prev->m_next_on_host = m_next_on_host ;
        }
        else {
          // Last view of created array, destroy it
          local_device_type::singleton().deallocate( *this );
        }

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
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  MDArrayViewRawData()
    : m_next_on_host( this )
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
    : m_next_on_host( this ) { insert_view( rhs ); }

  /** \brief  Assign a view of the array */
  inline
  MDArrayViewRawData & operator = ( const MDArrayViewRawData & rhs )
    { clear_view(); insert_view( rhs ); return *this ; }

  /**  \brief  Destroy this view of the array */
  inline
  ~MDArrayViewRawData() { clear_view(); }
};

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

} // namespace Kokkos

#endif /* KOKKOS_MAPPEDARRAY_HPP */

