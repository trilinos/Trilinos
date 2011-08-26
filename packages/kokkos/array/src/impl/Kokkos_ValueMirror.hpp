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

#ifndef KOKKOS_VALUEMIRROR_HPP
#define KOKKOS_VALUEMIRROR_HPP

#include <Kokkos_ValueView.hpp>
#include <Kokkos_ViewMirror.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename ValueType , class DeviceDst , class DeviceSrc >
class ViewTraits< ValueView< ValueType , DeviceDst > ,
                  ValueView< ValueType , DeviceSrc > >
{
public:

  typedef ValueView< ValueType , DeviceDst > ValueViewDst ;
  typedef ValueView< ValueType , DeviceSrc > ValueViewSrc ;

  /** \brief  Does their data reside in the same memory space ? */
  enum { same_memory_space =
           Impl::SameType< typename ValueViewDst::memory_space ,
                           typename ValueViewSrc::memory_space >::value };

  /** \brief  The two view types can view the same array */
  enum { compatible = same_memory_space };
};

//----------------------------------------------------------------------------

template< class ValueType , class DeviceDst , class DeviceSrc >
class ViewMirror< ValueView< ValueType , DeviceDst > ,
                  ValueView< ValueType , DeviceSrc > ,
                  false /* Keep deep copies */ >
{
public:
  static ValueView< ValueType , DeviceDst > ValueViewDst ;
  static ValueView< ValueType , DeviceSrc > ValueViewSrc ;

  static ValueViewDst create( const ValueViewSrc & s )
  { return create_labeled_value<ValueViewDst>( std::string() ); }

  static void update( const MDArrayViewDst & d , const MDArrayViewSrc & s )
  { Impl::ValueDeepCopy< ValueType , DeviceDst , DeviceSrc >::run( d , s ); }
};

//----------------------------------------------------------------------------

template< class ValueType , class DeviceDst , class DeviceSrc >
class ViewMirror< ValueView< ValueType , DeviceDst > ,
                  ValueView< ValueType , DeviceSrc > ,
                  true  /* Avoid deep copies */ >
{
public:
  static ValueView< ValueType , DeviceDst > ValueViewDst ;
  static ValueView< ValueType , DeviceSrc > ValueViewSrc ;

  inline
  static ValueViewDst create( const ValueViewSrc & s )
  { return ValueViewDst( s ); }

  static void update( const ValueViewDst & d , const ValueViewSrc & s )
  { if ( d != s ) { view_mirror_incompatible_throw(); } }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_VALUEMIRROR_HPP */

