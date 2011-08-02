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

template< class > class UnitTestDeviceMemoryManagement ;

template<>
class UnitTestDeviceMemoryManagement< KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE device ;

  typedef Kokkos::MemoryView< double , device::memory_space > dView ;
  typedef Kokkos::MemoryView< int ,    device::memory_space > iView ;

  static std::string name()
  {
    std::string tmp ;
    tmp.append( "UnitTestDeviceMemoryManagement< Kokkos::" );
    tmp.append( KOKKOS_MACRO_TO_STRING( KOKKOS_MACRO_DEVICE ) );
    tmp.append( " >" );
    return tmp ;
  }

  void error( const char * msg ) const
  {
    std::string tmp = name();
    tmp.append( msg );
    throw std::runtime_error( tmp );
  }

  UnitTestDeviceMemoryManagement()
  {
    std::ostringstream s_init ;
    std::ostringstream s_noop ;
    std::ostringstream s_alloc ;
    std::ostringstream s_assign ;
    std::ostringstream s_clear1 ;
    std::ostringstream s_clear3 ;
  
    device::memory_space::print_memory_view( s_init );

    dView dx , dy , dz ;
    iView ix , iy , iz ;
  
    if ( dx || dx.test_support_view_count() ||
         dy || dy.test_support_view_count() ||
         dz || dz.test_support_view_count() ||
         ix || ix.test_support_view_count() ||
         iy || iy.test_support_view_count() ||
         iz || iz.test_support_view_count() ) {
      error("FAILED Initialization");
    }
  
    device::memory_space::clear_memory_view( dx );
    device::memory_space::clear_memory_view( ix );
    device::memory_space::assign_memory_view( dy , dz );
    device::memory_space::assign_memory_view( iy , iz );
    device::memory_space::print_memory_view( s_noop );
  
    if ( dx || dx.test_support_view_count() ||
         dy || dy.test_support_view_count() ||
         dz || dz.test_support_view_count() ||
         ix || ix.test_support_view_count() ||
         iy || iy.test_support_view_count() ||
         iz || iz.test_support_view_count() ||
         s_init.str() != s_noop.str() ) {
      error("FAILED No-op");
    }
  
    device::memory_space::allocate_memory_view( dx , 10 , "dx" );
    device::memory_space::allocate_memory_view( ix , 20 , "ix" );
    device::memory_space::print_memory_view( s_alloc );
  
    if ( ! dx || 1 != dx.test_support_view_count() ||
         ! ix || 1 != ix.test_support_view_count() ||
         s_alloc.str().length() <= s_init.str().length() ) {
      error("FAILED Allocation");
    }
  
    device::memory_space::assign_memory_view( dy , dx );
    device::memory_space::assign_memory_view( iy , ix );
    device::memory_space::assign_memory_view( iz , iy );
    device::memory_space::print_memory_view( s_assign );
  
    if ( ! dx || 2 != dx.test_support_view_count() ||
         ! dy || 2 != dy.test_support_view_count() ||
         ! ix || 3 != ix.test_support_view_count() ||
         ! iy || 3 != iy.test_support_view_count() ||
         ! iz || 3 != iz.test_support_view_count() ||
         dx != dy ||
         ix != iy ||
         ix != iz ||
         s_assign.str() != s_alloc.str() ) {
      error("FAILED Assign view");
    }
  
    device::memory_space::clear_memory_view( dx );
    device::memory_space::clear_memory_view( ix );
    device::memory_space::print_memory_view( s_clear1 );
    
    if ( dx || dx.test_support_view_count() ||
         ix || ix.test_support_view_count() ||
         ! dy || 1 != dy.test_support_view_count() ||
         ! iy || 2 != iy.test_support_view_count() ||
         ! iz || 2 != iz.test_support_view_count() ||
         iy != iz ||
         s_clear1.str() != s_alloc.str() ) {
      error("FAILED Clear view #1");
    }
  
    device::memory_space::clear_memory_view( dy );
    device::memory_space::clear_memory_view( iy );
  
    if ( dx || dx.test_support_view_count() ||
         dy || dy.test_support_view_count() ||
         dz || dz.test_support_view_count() ||
         ix || ix.test_support_view_count() ||
         iy || iy.test_support_view_count() ||
         ! iz || 1 != iz.test_support_view_count() ) {
      error("FAILED Clear view #2");
    }
  
    device::memory_space::clear_memory_view( iz );
    device::memory_space::print_memory_view( s_clear3 );
  
    if ( dx || dx.test_support_view_count() ||
         dy || dy.test_support_view_count() ||
         dz || dz.test_support_view_count() ||
         ix || ix.test_support_view_count() ||
         iy || iy.test_support_view_count() ||
         iz || iz.test_support_view_count() ||
         s_init.str() != s_clear3.str() ) {
      error("FAILED Clear view #3");
    }
  
    // optional:

    // std::cout << name() << std::endl << s_alloc.str();
  }
};

}

/*--------------------------------------------------------------------------*/

