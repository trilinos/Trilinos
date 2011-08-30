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

#ifndef KOKKOS_MDARRAYMIRROR_HPP
#define KOKKOS_MDARRAYMIRROR_HPP

#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_ViewMirror.hpp>

namespace Kokkos {
namespace Impl {

template< typename ValueType , class DeviceDst , class DeviceSrc >
class ViewTraits< MDArrayView< ValueType , DeviceDst > ,
                  MDArrayView< ValueType , DeviceSrc > >
{
public:

  typedef MDArrayView< ValueType , DeviceDst > MDArrayViewDst ;
  typedef MDArrayView< ValueType , DeviceSrc > MDArrayViewSrc ;

  /** \brief  Does their data reside in the same memory space ? */
  enum { same_memory_space =
           Impl::SameType< typename MDArrayViewDst::memory_space ,
                           typename MDArrayViewSrc::memory_space >::value };

  /** \brief  Does their data have the same mapping ? */
  enum { same_mdarray_map =
           Impl::SameType< typename MDArrayViewDst::mdarray_map ,
                           typename MDArrayViewSrc::mdarray_map >::value };

  /** \brief  The two view types can view the same array */
  enum { compatible = same_memory_space && same_mdarray_map };

  /** \brief  Are the dimensions equal? */
  inline
  static bool equal_dimension( const MDArrayViewDst & dst ,
                               const MDArrayViewSrc & src )
    {
      typedef typename MDArrayViewDst::size_type size_type ;

      return dst.rank()       == (size_type) src.rank() &&
             dst.dimension(0) == (size_type) src.dimension(0) &&
             dst.dimension(1) == (size_type) src.dimension(1) &&
             dst.dimension(2) == (size_type) src.dimension(2) &&
             dst.dimension(3) == (size_type) src.dimension(3) &&
             dst.dimension(4) == (size_type) src.dimension(4) &&
             dst.dimension(5) == (size_type) src.dimension(5) &&
             dst.dimension(6) == (size_type) src.dimension(6) &&
             dst.dimension(7) == (size_type) src.dimension(7);
    }
};

//----------------------------------------------------------------------------

template< class ValueType , class DeviceDst , class DeviceSrc >
class ViewMirror< MDArrayView< ValueType , DeviceDst > ,
                  MDArrayView< ValueType , DeviceSrc > ,
                  false /* Keep deep copies */ >
{
public:
  typedef MDArrayView< ValueType , DeviceDst > MDArrayViewDst ;
  typedef MDArrayView< ValueType , DeviceSrc > MDArrayViewSrc ;

  static MDArrayViewDst create( const MDArrayViewSrc & s )
  {
    return create_labeled_mdarray<MDArrayViewDst>(
             std::string() ,
             s.dimension(0) , s.dimension(1) ,
             s.dimension(2) , s.dimension(3) ,
             s.dimension(4) , s.dimension(5) ,
             s.dimension(6) , s.dimension(7) );
  }

  static void update( const MDArrayViewDst & d , const MDArrayViewSrc & s )
  {
    const bool equal_dim =
      ViewTraits< MDArrayViewDst , MDArrayViewSrc >::equal_dimension( d , s );
    if ( ! equal_dim ) { view_mirror_incompatible_throw(); }
    Impl::MDArrayDeepCopy< ValueType , DeviceDst , DeviceSrc >::run( d , s );
  }
};

//----------------------------------------------------------------------------

template< class ValueType , class DeviceDst , class DeviceSrc >
class ViewMirror< MDArrayView< ValueType , DeviceDst > ,
                  MDArrayView< ValueType , DeviceSrc > ,
                  true  /* Avoid deep copies */ >
{
public:
  typedef MDArrayView< ValueType , DeviceDst > MDArrayViewDst ;
  typedef MDArrayView< ValueType , DeviceSrc > MDArrayViewSrc ;

  inline
  static MDArrayViewDst create( const MDArrayViewSrc & s )
  { return MDArrayViewDst( s ); }

  static void update( const MDArrayViewDst & d , const MDArrayViewSrc & s )
  { if ( d != s ) { view_mirror_incompatible_throw(); } }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_MDARRAYMIRROR_HPP */

