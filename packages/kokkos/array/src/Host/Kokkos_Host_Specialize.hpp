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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Partial specializations for the device

#include <Kokkos_Host_macros.hpp>

#if defined( KOKKOS_MEMORYVIEW_HPP ) && ! defined( KOKKOS_HOST_MEMORYVIEW )
#define KOKKOS_HOST_MEMORYVIEW
#include <impl/Kokkos_MemoryView_macros.hpp>
#endif

//----------------------------------------------------------------------------

#if defined( KOKKOS_VALUE_HPP ) && ! defined( KOKKOS_HOST_VALUE )
#define KOKKOS_HOST_VALUE
#include <impl/Kokkos_Value_macros.hpp>
#include <Host/Kokkos_Host_Value.hpp>
#endif

//----------------------------------------------------------------------------

#if defined( KOKKOS_MULTIVECTOR_HPP ) && ! defined( KOKKOS_HOST_MULTIVECTOR )
#define KOKKOS_HOST_MULTIVECTOR
#include <impl/Kokkos_MultiVector_macros.hpp>
#include <Host/Kokkos_Host_MultiVector.hpp>
#endif

//----------------------------------------------------------------------------

#if defined( KOKKOS_CRSMAP_HPP ) && ! defined( KOKKOS_HOST_CRSMAP )
#define KOKKOS_HOST_CRSMAP
#include <impl/Kokkos_CrsMap_macros.hpp>
#include <Host/Kokkos_Host_CrsMap.hpp>
#endif

//----------------------------------------------------------------------------

#if defined( KOKKOS_MDARRAY_HPP ) && ! defined( KOKKOS_HOST_MDARRAY )
#define KOKKOS_HOST_MDARRAY
#include <Host/Kokkos_Host_MDArray.hpp>
#include <impl/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <impl/Kokkos_MDArrayIndexMapRight_macros.hpp>

// For MDArray< ValueType , Host >
#include <impl/Kokkos_MDArray_macros.hpp>

// For MDArray< ValueType , HostMapped< MDArrayMap > >
#undef  KOKKOS_MACRO_DEVICE
#define KOKKOS_MACRO_DEVICE  Impl::HostMapped< MDArrayMap >
#define KOKKOS_MACRO_MDARRAY_TEMPLATE_ARGUMENT  class MDArrayMap
#include <impl/Kokkos_MDArray_macros.hpp>
#undef KOKKOS_MACRO_MDARRAY_TEMPLATE_ARGUMENT
#endif

//----------------------------------------------------------------------------

#include <Kokkos_Clear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


