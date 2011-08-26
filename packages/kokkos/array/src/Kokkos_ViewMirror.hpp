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

#ifndef KOKKOS_VIEWMIRROR_HPP
#define KOKKOS_VIEWMIRROR_HPP

#if ! defined( KOKKOS_MIRROR_VIEW_OPTIMIZE )
#define KOKKOS_MIRROR_VIEW_OPTIMIZE 0
#endif /* ! defined( KOKKOS_MIRROR_VIEW_OPTIMIZE ) */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class ViewDest , class ViewSrc > class ViewTraits ;

template< class ViewDst , class ViewSrc ,
          bool AvoidDeepCopy = ViewTraits< ViewDst , ViewSrc >::compatible &&
                               KOKKOS_MIRROR_VIEW_OPTIMIZE >
class ViewMirror ;

void view_mirror_incompatible_throw();

} // namespace Impl

template< class ViewType >
typename ViewType::HostView mirror_create( const ViewType & s )
{ return Impl::ViewMirror< typename ViewType::HostView , ViewType >::create( s ); }

template< class ViewTypeDst , class ViewTypeSrc >
void mirror_update( const ViewTypeDst & d , const ViewTypeSrc & s )
{ return Impl::ViewMirror< ViewTypeDst , ViewTypeSrc >::update( d , s ); }

} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* KOKKOS_VIEWMIRROR_HPP */

