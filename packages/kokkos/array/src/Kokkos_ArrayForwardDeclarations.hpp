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

#ifndef KOKKOS_ARRAYFORWARDDECLARATIONS_HPP
#define KOKKOS_ARRAYFORWARDDECLARATIONS_HPP

namespace Kokkos {

//----------------------------------------------------------------------------
// Devices

class DeviceHost ;
class DeviceTPI ;
class DeviceCuda ;

//----------------------------------------------------------------------------
// Simple value

template< typename ValueType , class DeviceType >
class ValueView ;

template< typename ValueType , class DeviceType >
void deep_copy( const ValueView< ValueType , DeviceType > & dst ,
                const ValueType & src );

template< typename ValueType , class DeviceType >
void deep_copy( ValueType & dst ,
                const ValueView< ValueType , DeviceType > & src );

template< typename ValueType , class DeviceDst , class DeviceSrc >
void deep_copy( const ValueView< ValueType , DeviceDst > & dst ,
                const ValueView< ValueType , DeviceSrc > & src );

//----------------------------------------------------------------------------
// Multivector

template< typename ValueType , class DeviceType >
class MultiVectorView ;

template< typename ValueType , class DeviceDst , class DeviceSrc >
void deep_copy( const MultiVectorView< ValueType , DeviceDst > & dst ,
                const MultiVectorView< ValueType , DeviceSrc > & src );

//----------------------------------------------------------------------------
// Multidimensional array and mapping options

enum { MDArrayMaxRank = 8 };

template< typename ValueType ,
          class DeviceType ,
          class MapOption = typename DeviceType::default_mdarray_map >
class MDArrayView ;

template< typename ValueType , class DeviceDst , class MapDst ,
                               class DeviceSrc , class MapSrc >
void deep_copy( const MDArrayView< ValueType , DeviceDst , MapDst > & dst ,
                const MDArrayView< ValueType , DeviceSrc , MapSrc > & src );

class MDArrayIndexMapLeft ;
class MDArrayIndexMapRight ;

//----------------------------------------------------------------------------
// View to device memory

template< typename ValueType , class DeviceType >
class MemoryView ;

//----------------------------------------------------------------------------

} // namespace Kokkos

#endif /* #ifndef KOKKOS_ARRAYFORWARDDECLARATIONS_HPP */


