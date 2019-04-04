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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
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
  typedef double scalar_type;
#elif defined(HAVE_TPETRA_INST_FLOAT)
  typedef float scalar_type;
#else
#  error "Tpetra: No scalar types in the set {float, double} have been enabled."
#endif

  //! Default value of LocalOrdinal template parameter.
  typedef int local_ordinal_type;

  /// \typedef global_ordinal_type
  /// \brief Default value of GlobalOrdinal template parameter.
#if defined(HAVE_TPETRA_INST_INT_INT)
  typedef int global_ordinal_type;
#elif defined(HAVE_TPETRA_INST_INT_LONG_LONG)
  typedef long long global_ordinal_type;
#elif defined(HAVE_TPETRA_INST_INT_LONG)
  typedef long global_ordinal_type;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED_LONG)
  typedef unsigned long global_ordinal_type;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED)
  typedef unsigned global_ordinal_type;
#else
#  error "Tpetra: No global ordinal types in the set {int, long long, long, unsigned long, unsigned} have been enabled."
#endif

  /// \typedef execution_space
  /// \brief Default Tpetra execution space.
#if defined(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE)
  using execution_space = ::Kokkos::Cuda;
#elif defined(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE)
  using execution_space = ::Kokkos::OpenMP;
#elif defined(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE)
  using execution_space = ::Kokkos::Threads;
#elif defined(HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE)
  using execution_space = ::Kokkos::Serial;
#else
#    error "No default Tpetra Node type specified.  Please set the CMake option Tpetra_DefaultNode to a valid Node type."
#endif

  //! Default value of Node template parameter.
  using node_type = ::Kokkos::Compat::KokkosDeviceWrapperNode<execution_space>;
} // namespace DefaultTypes

} // namespace Details

} // namespace Tpetra

#endif // TPETRA_DETAILS_DEFAULTTYPES_HPP
