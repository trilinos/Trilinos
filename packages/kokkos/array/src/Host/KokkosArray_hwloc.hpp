/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_HWLOC_HPP
#define KOKKOSARRAY_HWLOC_HPP

#include <vector>
#include <ostream>

namespace KokkosArray {
namespace Impl {

struct hwloc {

  static const unsigned max_depth = 4 ;

  /** \brief  Hierarchy of hardware thread capacity.
   *  
   *  Total thread capacity is product of non-zero values.
   */
  static void get_thread_capacity( unsigned capacity[] );

  static void print_thread_capacity( std::ostream & );

  static unsigned get_thread_capacity_depth();

  enum BindingPolicy { SPREAD , PACK };

  /** \brief  Query coordinate to bind thread according to policy */
  static void map_thread( const BindingPolicy policy ,
                          const unsigned      rank ,
                          const unsigned      count ,
                                unsigned      coordinate[] );


  static void map_thread( const BindingPolicy policy ,
                          const unsigned      capacity_depth ,
                          const unsigned      capacity[] ,
                          const unsigned      rank ,
                          const unsigned      count ,
                                unsigned      coordinate[] );


  /** \brief  Bind the current thread to the hierarchical coordinate */
  static bool bind_this_thread( const unsigned coordinate[] );
};

} // namespace Impl
} // namespace KokkosArray

#endif /* #define KOKKOSARRAY_HWLOC_HPP */

