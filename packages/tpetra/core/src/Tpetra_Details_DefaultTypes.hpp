// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_DEFAULTTYPES_HPP
#define TPETRA_DETAILS_DEFAULTTYPES_HPP

#include "TpetraCore_config.h"
#include "KokkosClassic_DefaultNode_config.h"
#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"

//! Namespace for Tpetra classes and methods
namespace Tpetra {

/// \brief Namespace for Tpetra implementation details.
/// \warning Do NOT rely on the contents of this namespace.
namespace Details {

/// \brief Declarations of values of Tpetra classes' default template parameters.
///
/// \warning Don't use this directly.  Get defaults from Tpetra classes.
///   For example: <tt>Tpetra::MultiVector<>::scalar_type</tt>,
///   <tt>Tpetra::Map<>::local_ordinal_type</tt>.
namespace DefaultTypes {
  //! Default value of Scalar template parameter.
#if defined(HAVE_TPETRA_INST_DOUBLE)
  using scalar_type = double;
#elif defined(HAVE_TPETRA_INST_FLOAT)
  using scalar_type = float;
#else
#  error "Tpetra: No scalar types in the set {float, double} have been enabled."
#endif

  //! Default value of LocalOrdinal template parameter.
  using local_ordinal_type = int;

  /// \typedef global_ordinal_type
  /// \brief Default value of GlobalOrdinal template parameter.
#if defined(HAVE_TPETRA_INST_INT_LONG_LONG)
    using global_ordinal_type = long long;
#elif defined(HAVE_TPETRA_INST_INT_INT)
    using global_ordinal_type = int;
#elif defined(HAVE_TPETRA_INST_INT_LONG)
    using global_ordinal_type = long;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED_LONG)
    using global_ordinal_type = unsigned long;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED)
    using global_ordinal_type = unsigned;
#else
    #error "Tpetra: No global ordinal types in the set {int, long long, long, unsigned long, unsigned} have been enabled."
#endif

  /// \typedef execution_space
  /// \brief Default Tpetra execution space and Node type.
#if defined(HAVE_TPETRA_DEFAULTNODE_HIPWRAPPERNODE)
  using execution_space = ::Kokkos::Experimental::HIP;
  using node_type = ::Kokkos::Compat::KokkosHIPWrapperNode;
#elif defined(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE)
  using execution_space = ::Kokkos::Cuda;
  using node_type = ::Kokkos::Compat::KokkosCudaWrapperNode;
#elif defined(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE)
  using execution_space = ::Kokkos::OpenMP;
  using node_type = ::Kokkos::Compat::KokkosOpenMPWrapperNode;
#elif defined(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE)
  using execution_space = ::Kokkos::Threads;
  using node_type = ::Kokkos::Compat::KokkosThreadsWrapperNode;
#elif defined(HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE)
  using execution_space = ::Kokkos::Serial;
  using node_type = ::Kokkos::Compat::KokkosSerialWrapperNode;
#else
#    error "No default Tpetra Node type specified.  Please set the CMake option Tpetra_DefaultNode to a valid Node type."
#endif

  /// \brief Memory space used for MPI communication buffers.
  ///
  /// See #1088 for why this is not just ExecutionSpace::memory_space

  template<typename ExecutionSpace>
  struct CommBufferMemorySpace
  {
    using type = typename ExecutionSpace::memory_space;
  };

#ifdef KOKKOS_ENABLE_CUDA
  template<>
  struct CommBufferMemorySpace<Kokkos::Cuda>
  {
    using type = Kokkos::CudaSpace;
  };
#endif

#ifdef KOKKOS_ENABLE_HIP
  template<>
  struct CommBufferMemorySpace<Kokkos::Experimental::HIPHostPinnedSpace>
  {
    using type = Kokkos::Experimental::HIPHostPinnedSpace;
  };
#endif

  template<typename Device>
  using comm_buffer_memory_space = typename CommBufferMemorySpace<typename Device::execution_space>::type;

} // namespace DefaultTypes

} // namespace Details

} // namespace Tpetra

#endif // TPETRA_DETAILS_DEFAULTTYPES_HPP
