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

#ifndef KOKKOS_CREATEMIRROR_HPP
#define KOKKOS_CREATEMIRROR_HPP

#if ! defined( KOKKOS_MIRROR_VIEW_OPTIMIZE )
#define KOKKOS_MIRROR_VIEW_OPTIMIZE 0
#endif /* ! defined( KOKKOS_MIRROR_VIEW_OPTIMIZE ) */

#include <impl/Kokkos_StaticAssert.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class ViewType , bool Optimize > class CreateMirror ;

} // namespace Impl

template< class ViewType >
typename ViewType::HostView create_mirror( const ViewType & s )
{
  typedef typename ViewType::device_type      device ;
  typedef typename device::memory_space       memory ;
  typedef typename ViewType::HostView         host_view ;
  typedef typename host_view::device_type     host_device ;
  typedef typename host_device::memory_space  host_memory ;

  enum { compatible = Impl::SameType< memory , host_memory >::value };
  enum { optimize   = compatible && KOKKOS_MIRROR_VIEW_OPTIMIZE };

  return Impl::CreateMirror< ViewType , optimize >::create( s );
}

} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* KOKKOS_CREATEMIRROR_HPP */

