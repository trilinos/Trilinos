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

#ifndef KOKKOS_MACRO_DEVICE
#error "KOKKOS_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <impl/Kokkos_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace {

template< typename T, class > class UnitTestMDArrayDeepCopy ;

template<typename T>
class UnitTestMDArrayDeepCopy< T, KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE device ;
  typedef Kokkos::DeviceHost  host ;

  typedef Kokkos::MDArrayView< T , host > host_dView ;
  typedef Kokkos::MDArrayView< T , device > dView ;

  UnitTestMDArrayDeepCopy() { run_test(); }

  void run_test()
  {
    enum { dN = 1000 };

    dView dx , dy ;
    host_dView host_dx , host_dy ;

    host_dx = Kokkos::create_labeled_mdarray< host_dView > ( "dx" , dN );
    host_dy = Kokkos::create_labeled_mdarray< host_dView > ( "dy" , dN );

    dx = Kokkos::create_labeled_mdarray< dView > ( "dx" , dN );
    dy = Kokkos::create_labeled_mdarray< dView > ( "dy" , dN );

    for ( size_t i = 0 ; i < dN ; ++i ) { host_dx(i) = i ; }

    Kokkos::deep_copy( dx , host_dx );
    Kokkos::deep_copy( dy , dx );
    Kokkos::deep_copy( host_dy , dy );


    for ( size_t i = 0 ; i < dN ; ++i ) {
      ASSERT_EQ(host_dx(i), host_dy(i));
    }
  }
};

}

/*--------------------------------------------------------------------------*/

