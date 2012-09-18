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

#ifndef KOKKOSARRAY_IMPLWALLTIME_HPP
#define KOKKOSARRAY_IMPLWALLTIME_HPP

#include <stddef.h>

#ifdef _MSC_VER
#include <gettimeofday.c>
#else
#include <sys/time.h>
#endif

namespace KokkosArray {
namespace Impl {

/** \brief  Time since construction */

class Timer {
private:
  struct timeval m_old ;
  Timer( const Timer & );
  Timer & operator = ( const Timer & );
public:

  inline
  void reset() { gettimeofday( & m_old , ((struct timezone *) NULL ) ); }

  inline
  ~Timer() {}

  inline
  Timer() { reset(); }

  inline
  double seconds() const
  {
    struct timeval m_new ;

    ::gettimeofday( & m_new , ((struct timezone *) NULL ) );

    return ( (double) ( m_new.tv_sec  - m_old.tv_sec ) ) +
           ( (double) ( m_new.tv_usec - m_old.tv_usec ) * 1.0e-6 );
  }
};

} // namespace Impl
} // namespace KokkosArray

#endif /* #ifndef KOKKOSARRAY_IMPLWALLTIME_HPP */

