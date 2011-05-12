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
#include <macros/Kokkos_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace {

template< class > class UnitTestDeviceMemoryManagement ;

template<>
class UnitTestDeviceMemoryManagement< Kokkos :: KOKKOS_MACRO_DEVICE >
{
public:
  typedef Kokkos:: KOKKOS_MACRO_DEVICE device ;

  typedef Kokkos::MemoryView< double , device > dView ;
  typedef Kokkos::MemoryView< int ,    device > iView ;

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
    throw std::runtime_error( msg );
  }

  UnitTestDeviceMemoryManagement()
  {
    std::ostringstream s_empty ;
    std::ostringstream s_alloc ;
    std::ostringstream s_assign ;
    std::ostringstream s_clear1 ;
  
    dView dx , dy , dz ;
    iView ix , iy , iz ;
  
    if ( dx.ptr_on_device() || dx.test_support_view_count() ||
         dy.ptr_on_device() || dy.test_support_view_count() ||
         dz.ptr_on_device() || dz.test_support_view_count() ||
         ix.ptr_on_device() || ix.test_support_view_count() ||
         iy.ptr_on_device() || iy.test_support_view_count() ||
         iz.ptr_on_device() || iz.test_support_view_count() ) {
      error("FAILED Initialization");
    }
  
    device::clear_memory_view( dx );
    device::clear_memory_view( ix );
    device::assign_memory_view( dy , dz );
    device::assign_memory_view( iy , iz );
    device::print_memory_view( s_empty );
  
    if ( dx.ptr_on_device() || dx.test_support_view_count() ||
         dy.ptr_on_device() || dy.test_support_view_count() ||
         dz.ptr_on_device() || dz.test_support_view_count() ||
         ix.ptr_on_device() || ix.test_support_view_count() ||
         iy.ptr_on_device() || iy.test_support_view_count() ||
         iz.ptr_on_device() || iz.test_support_view_count() ||
         s_empty.str().length() ) {
      error("FAILED No-op");
    }
  
    device::allocate_memory_view( dx , 10 , "dx" );
    device::allocate_memory_view( ix , 20 , "ix" );
    device::print_memory_view( s_alloc );
  
    if ( ! dx.ptr_on_device() || 1 != dx.test_support_view_count() ||
         ! ix.ptr_on_device() || 1 != ix.test_support_view_count() ||
         ! s_alloc.str().length() ) {
      error("FAILED Allocation");
    }
  
    device::assign_memory_view( dy , dx );
    device::assign_memory_view( iy , ix );
    device::assign_memory_view( iz , iy );
    device::print_memory_view( s_assign );
  
    if ( ! dx.ptr_on_device() || 2 != dx.test_support_view_count() ||
         ! dy.ptr_on_device() || 2 != dy.test_support_view_count() ||
         ! ix.ptr_on_device() || 3 != ix.test_support_view_count() ||
         ! iy.ptr_on_device() || 3 != iy.test_support_view_count() ||
         ! iz.ptr_on_device() || 3 != iz.test_support_view_count() ||
         dx.ptr_on_device() != dy.ptr_on_device() ||
         ix.ptr_on_device() != iy.ptr_on_device() ||
         ix.ptr_on_device() != iz.ptr_on_device() ||
         s_assign.str() != s_alloc.str() ) {
      error("FAILED Assign view");
    }
  
    device::clear_memory_view( dx );
    device::clear_memory_view( ix );
    device::print_memory_view( s_clear1 );
    
    if ( dx.ptr_on_device() || dx.test_support_view_count() ||
         ix.ptr_on_device() || ix.test_support_view_count() ||
         ! dy.ptr_on_device() || 1 != dy.test_support_view_count() ||
         ! iy.ptr_on_device() || 2 != iy.test_support_view_count() ||
         ! iz.ptr_on_device() || 2 != iz.test_support_view_count() ||
         iy.ptr_on_device() != iz.ptr_on_device() ||
         s_clear1.str() != s_alloc.str() ) {
      error("FAILED Clear view #1");
    }
  
    device::clear_memory_view( dy );
    device::clear_memory_view( iy );
  
    if ( dx.ptr_on_device() || dx.test_support_view_count() ||
         dy.ptr_on_device() || dy.test_support_view_count() ||
         dz.ptr_on_device() || dz.test_support_view_count() ||
         ix.ptr_on_device() || ix.test_support_view_count() ||
         iy.ptr_on_device() || iy.test_support_view_count() ||
         ! iz.ptr_on_device() || 1 != iz.test_support_view_count() ) {
      error("FAILED Clear view #2");
    }
  
    device::clear_memory_view( iz );
    device::print_memory_view( s_empty );
  
    if ( dx.ptr_on_device() || dx.test_support_view_count() ||
         dy.ptr_on_device() || dy.test_support_view_count() ||
         dz.ptr_on_device() || dz.test_support_view_count() ||
         ix.ptr_on_device() || ix.test_support_view_count() ||
         iy.ptr_on_device() || iy.test_support_view_count() ||
         iz.ptr_on_device() || iz.test_support_view_count() ||
         s_empty.str().length() ) {
      error("FAILED Clear view #3");
    }
  
    // optional:

    // std::cout << name() << std::endl << s_alloc.str();
  }
};

}

/*--------------------------------------------------------------------------*/

