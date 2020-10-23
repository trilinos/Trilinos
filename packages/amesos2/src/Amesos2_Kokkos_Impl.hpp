// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_Kokkos_Impl.hpp
  \author
  \date

  \brief  ETI for Solvers using Kokkos adapter
*/

#ifndef AMESOS2_KOKKOS_IMPL_HPP
#define AMESOS2_KOKKOS_IMPL_HPP

#include <type_traits>
#include "Amesos2_KokkosMultiVecAdapter_decl.hpp"
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#define AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(S,LO,EXEC_SPACE)                             \
  template class Amesos2::AMESOS2_KOKKOS_IMPL_SOLVER_NAME<KokkosSparse::CrsMatrix<S, LO,         \
    typename EXEC_SPACE::device_type>,                                                           \
    Kokkos::View<S**, Kokkos::LayoutLeft, typename EXEC_SPACE::device_type> >;

#ifdef KOKKOS_ENABLE_CUDA_UVM
#define AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(S,LO)       \
  template class Amesos2::AMESOS2_KOKKOS_IMPL_SOLVER_NAME<KokkosSparse::CrsMatrix<S, LO,          \
    Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>>,                                              \
    Kokkos::View<S**, Kokkos::LayoutLeft, Kokkos::Device<Kokkos::Cuda,Kokkos::CudaSpace>> >;
#else
#define AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(S,LO)
#endif

#if defined(KOKKOS_ENABLE_SERIAL)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Kokkos::Serial)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Kokkos::Serial)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Kokkos::Serial)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Kokkos::Serial)
#endif
#endif

#if defined(KOKKOS_ENABLE_THREADS)
#define EXEC_SPACE Kokkos::Threads
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Kokkos::Threads)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Kokkos::Threads)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Kokkos::Threads)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Kokkos::Threads)
#endif
#endif // KOKKOS_ENABLE_THREADS

#if defined(KOKKOS_ENABLE_OPENMP)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Kokkos::OpenMP)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Kokkos::OpenMP)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Kokkos::OpenMP)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Kokkos::OpenMP)
#endif
#endif // KOKKOS_ENABLE_OPENMP

#if defined(KOKKOS_ENABLE_CUDA)
#ifdef HAVE_TPETRA_INST_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(float, int, Kokkos::Cuda)
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(float, int)
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(double, int, Kokkos::Cuda)
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(double, int)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<float>, int, Kokkos::Cuda)
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(Kokkos::complex<float>, int)
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER(Kokkos::complex<double>, int, Kokkos::Cuda)
    AMESOS2_KOKKOS_LOCAL_INSTANT_KOKKOS_ADAPTER_UVM_OFF(Kokkos::complex<double>, int)
#endif
#endif // KOKKOS_ENABLE_CUDA

#endif // AMESOS2_KOKKOS_IMPL_HPP
