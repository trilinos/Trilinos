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

#ifndef KOKKOS_MULTIVECTORMIRROR_HPP
#define KOKKOS_MULTIVECTORMIRROR_HPP

#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_ViewMirror.hpp>

namespace Kokkos {
namespace Impl {

template< typename ValueType , class DeviceDst , class DeviceSrc >
class ViewTraits< MultiVectorView< ValueType , DeviceDst > ,
                  MultiVectorView< ValueType , DeviceSrc > >
{
public:

  typedef MultiVectorView< ValueType , DeviceDst > MultiVectorViewDst ;
  typedef MultiVectorView< ValueType , DeviceSrc > MultiVectorViewSrc ;

  /** \brief  Does their data reside in the same memory space ? */
  enum { same_memory_space =
           Impl::SameType< typename MultiVectorViewDst::memory_space ,
                           typename MultiVectorViewSrc::memory_space >::value };

  /** \brief  The two view types can view the same array */
  enum { compatible = same_memory_space };

  /** \brief  Are the dimensions equal? */
  inline
  static bool equal_dimension( const MultiVectorViewDst & dst ,
                               const MultiVectorViewSrc & src )
    {
      typedef typename MultiVectorViewDst::size_type size_type ;

      return dst.length() == (size_type) src.length() &&
             dst.count()  == (size_type) src.count();
    }
};

//----------------------------------------------------------------------------

template< class ValueType , class DeviceDst , class DeviceSrc >
class ViewMirror< MultiVectorView< ValueType , DeviceDst > ,
                  MultiVectorView< ValueType , DeviceSrc > ,
                  false /* Keep deep copies */ >
{
public:
  typedef MultiVectorView< ValueType , DeviceDst > MultiVectorViewDst ;
  typedef MultiVectorView< ValueType , DeviceSrc > MultiVectorViewSrc ;

  static MultiVectorViewDst create( const MultiVectorViewSrc & s )
  {
    return create_labeled_multivector<MultiVectorViewDst>(
             std::string() , s.length(), s.count() );
  }

  static void update( const MultiVectorViewDst & d , const MultiVectorViewSrc & s )
  {
    const bool equal_dim =
      ViewTraits< MultiVectorViewDst , MultiVectorViewSrc >::equal_dimension( d , s );
    if ( ! equal_dim ) { view_mirror_incompatible_throw(); }
    Impl::MultiVectorDeepCopy< ValueType , DeviceDst , DeviceSrc >::run( d , s );
  }
};

//----------------------------------------------------------------------------

template< class ValueType , class DeviceDst , class DeviceSrc >
class ViewMirror< MultiVectorView< ValueType , DeviceDst > ,
                  MultiVectorView< ValueType , DeviceSrc > ,
                  true  /* Avoid deep copies */ >
{
public:
  typedef MultiVectorView< ValueType , DeviceDst > MultiVectorViewDst ;
  typedef MultiVectorView< ValueType , DeviceSrc > MultiVectorViewSrc ;

  inline
  static MultiVectorViewDst create( const MultiVectorViewSrc & s )
  { return MultiVectorViewDst( s ); }

  static void update( const MultiVectorViewDst & d , const MultiVectorViewSrc & s )
  { if ( d != s ) { view_mirror_incompatible_throw(); } }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_MULTIVECTORMIRROR_HPP */

