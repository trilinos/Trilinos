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

#include <string>
#include <stdexcept>
#include <cstdlib>
#include <ostream>
#include <Kokkos_CudaMap.hpp>

namespace Kokkos {

/*--------------------------------------------------------------------------*/

CudaDevice::CudaDevice()
{}

CudaDevice::~CudaDevice()
{}

void CudaDevice::deallocate( void * pointer )
{
  cudaFree( pointer );
}

void * CudaDevice::allocate( CudaDevice::size_type  sizeof_value ,
                             CudaDevice::size_type  chunk_count ,
                             CudaDevice::size_type  work_count )
{
  void * pointer = NULL ;

  // TODO: Error check return value
  cudaMalloc( & pointer , chunk_count * work_count * sizeof_value );

  // TODO: Initialize memory to zero.

  return pointer ;
}

/*--------------------------------------------------------------------------*/

CudaMap::CudaMap( CudaMap::size_type parallel_work_count )
  : m_device()
  , m_arrays( m_device , parallel_work_count )
{}

CudaMap::~CudaMap()
{}

} // namespace Kokkos


