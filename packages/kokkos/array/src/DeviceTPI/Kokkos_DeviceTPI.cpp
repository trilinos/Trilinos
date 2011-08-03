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

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <TPI.h>
#include <Kokkos_DeviceTPI.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace {

class DeviceTPI_Impl {
public:

  ~DeviceTPI_Impl();

  static DeviceTPI_Impl & singleton();
};

DeviceTPI_Impl & DeviceTPI_Impl::singleton()
{
  static DeviceTPI_Impl self ;
  return self ;
}

DeviceTPI_Impl::~DeviceTPI_Impl()
{}

}

/*--------------------------------------------------------------------------*/

void DeviceTPI::initialize( size_type nthreads )
{
  TPI_Init( nthreads );
}

void DeviceTPI::finalize()
{
  TPI_Finalize();
}

/*--------------------------------------------------------------------------*/

unsigned int DeviceTPI::m_launching_kernel = false ;

void DeviceTPI::set_dispatch_functor()
{
  if ( m_launching_kernel ) {
    std::string msg ;
    msg.append( "Kokkos::DeviceTPI::set_dispatch_functor FAILED: " );
    msg.append( "kernel dispatch is already in progress, " );
    msg.append( "a recursive call or forgotten 'clear_dispatch_kernel" );
    throw std::runtime_error( msg );
  }
  m_launching_kernel = true ;
}

void DeviceTPI::clear_dispatch_functor()
{
  if ( ! m_launching_kernel ) {
    std::string msg ;
    msg.append( "Kokkos::DeviceTPI::clear_dispatch_functor FAILED: " );
    msg.append( "no kernel dispatch in progress." );
    throw std::runtime_error( msg );
  }
  m_launching_kernel = false ;
}


} // namespace Kokkos

