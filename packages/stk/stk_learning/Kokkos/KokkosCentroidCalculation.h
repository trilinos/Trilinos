/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef _KokkosCentroidCalculation_h_
#define _KokkosCentroidCalculation_h_

#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/ioUtils.hpp>
#include "mtk_kokkos.h"
#include <stk_util/stk_config.h>


namespace {

#if KOKKOS_ENABLE_CUDA
typedef double my_double;
#else
typedef long double my_double;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
typedef Kokkos::OpenMP   ExecSpace ;
#elif KOKKOS_ENABLE_CUDA
typedef Kokkos::Cuda     ExecSpace ;
#else
typedef Kokkos::Serial   ExecSpace ;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
typedef Kokkos::OpenMP       MemSpace;
#elif KOKKOS_ENABLE_CUDA
typedef Kokkos::CudaSpace    MemSpace;
#else
typedef Kokkos::HostSpace    MemSpace;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
typedef Kokkos::OpenMP       UVMMemSpace;
#elif KOKKOS_ENABLE_CUDA
typedef Kokkos::CudaUVMSpace UVMMemSpace;
#else
typedef Kokkos::HostSpace    UVMMemSpace;
#endif


typedef Kokkos::RangePolicy<ExecSpace> range_policy ;

typedef Kokkos::View<my_double*, MemSpace>   DeviceViewVectorType;
typedef Kokkos::View<my_double*, Kokkos::HostSpace>   HostViewVectorType;

typedef Kokkos::TeamPolicy<ExecSpace>               team_policy ;
typedef Kokkos::TeamPolicy<ExecSpace>::member_type  member_type ;

#ifdef KOKKOS_ENABLE_CUDA
typedef Kokkos::LayoutLeft   Layout ;
#else
typedef Kokkos::LayoutRight   Layout ;
#endif

typedef Kokkos::View<double**, Layout, MemSpace>   DeviceViewMatrixType;
typedef Kokkos::View<const double**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess> >   ConstDeviceViewMatrixType;
//typedef Kokkos::View<double**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::Atomic> >   DeviceViewAtomicMatrixType;
typedef Kokkos::View<double**, Layout, MemSpace>   DeviceViewAtomicMatrixType;

typedef Kokkos::View<double**, Layout, MemSpace>   DeviceArray2DType;
typedef Kokkos::View<const double**, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess> >   ConstDeviceArray2DType;

typedef Kokkos::View<double***, Layout, MemSpace>   DeviceArray3DType;
typedef Kokkos::View<const double***, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess> >   ConstDeviceArray3DType;

typedef Kokkos::View<int*, Layout, MemSpace> DeviceViewIntType;
typedef Kokkos::View<const int*, Layout, MemSpace, Kokkos::MemoryTraits<Kokkos::RandomAccess>> ConstDeviceViewIntType;


} // namespace

#endif
