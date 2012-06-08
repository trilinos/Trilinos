/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_IMPL_INDEXMAP_HPP
#define KOKKOS_IMPL_INDEXMAP_HPP

#include <impl/KokkosArray_ArrayTraits.hpp>

namespace KokkosArray {
namespace Impl {

enum { IndexMapMaxRank = 8 };

/** \brief  Mapping of a multi-index to an offset.
 *
 *  If the Rank is zero then then all dimensions are runtime knowledge.
 *  If the Rank is non-zero then all but the leading dimension is
 *  compile-time knowledge from the template arguments.
 *
 *  The leading dimension is always set at runtime.
 */
template< class MemorySpace , unsigned Rank ,
          unsigned N1 = 0, unsigned N2 = 0, unsigned N3 = 0,
          unsigned N4 = 0, unsigned N5 = 0, unsigned N6 = 0, unsigned N7 = 0 >
class IndexMapRight ;

template< class MemorySpace , unsigned Rank ,
          unsigned N1 = 0, unsigned N2 = 0, unsigned N3 = 0,
          unsigned N4 = 0, unsigned N5 = 0, unsigned N6 = 0, unsigned N7 = 0 >
class IndexMapLeft ;

} // namespace Impl
} // namespace KokkosArray

#endif /* KOKKOS_IMPL_INDEXMAP_HPP */


