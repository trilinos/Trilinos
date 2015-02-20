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

#ifndef KOKKOS_DEFAULT_KERNELS_HPP
#define KOKKOS_DEFAULT_KERNELS_HPP

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_AltSparseOps.hpp"
#include "Kokkos_DefaultSparseOps.hpp"
#include "Kokkos_DefaultRelaxation.hpp"

namespace Kokkos {
  namespace Compat {
    // Forward declaration (to avoid circular subpackage dependencies).
    template<class DeviceType>
    class KokkosDeviceWrapperNode;
  } // namespace Compat
  // Forward declaration (to avoid circular subpackage dependencies).
  class Cuda;
} // namespace Kokkos

namespace KokkosClassic {

  /// \brief Traits class providing default kernel types for CRS,
  ///   block CRS and relaxation kernels.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries in the sparse matrix; same as
  ///   the Scalar template parameter of Tpetra objects.
  /// \tparam Ordinal The type of indices in the sparse matrix; same
  ///   as the LocalOrdinal template parameter of Tpetra objects.
  /// \tparam Node The Kokkos Node type; same as the Node template
  ///   parameter of Tpetra objects.
  ///
  /// This class provides typedefs to implementations of local sparse
  /// matrix kernels.  Its typedefs are used by Tpetra objects by
  /// default (hence the name).  If you wish your Tpetra objects to
  /// use different kernels, you may implement your own class with the
  /// necessary typedefs.
  ///
  /// This class includes at least the "SparseOps" typedef, but may
  /// otherwise include any of the following typedefs:
  /// - SparseOps: A valid fifth template parameter for
  ///   Tpetra::CrsGraph and Tpetra::CrsMatrix.  Implements local
  ///   sparse matrix-(multi)vector multiply and sparse triangular
  ///   solve, for a CrsMatrix.
  /// - Relaxations: Implements local relaxation kernels.
  ///
  /// It is not necessary to implement a DefaultKernels-like class
  /// with typedefs if you wish your Tpetra objects to use different
  /// kernels.  We provide DefaultKernels mainly as a convenience for
  /// Tpetra developers, so they can find all the default kernels in
  /// one place.  We also use specializations of DefaultKernels for
  /// various Scalar, Ordinal, and Node types to ensure that Tpetra
  /// uses high-performance TPLs whenever possible.
  template <class Scalar, class Ordinal, class Node>
  struct DefaultKernels {
    typedef DefaultHostSparseOps<void, Ordinal, Node> SparseOps;
    typedef DefaultRelaxation<Scalar, Ordinal, Node> Relaxations;
  };

#if defined(HAVE_TPETRACLASSIC_SERIAL)
  /// \brief Partial specialization for Node=SerialNode.
  ///
  /// AltSparseOps doesn't use KokkosClassic's parallel programming
  /// programming model, doesn't do deep copies for first touch, and
  /// doesn't rely so heavily on inlining.  Thus, it's a reasonable
  /// choice when not using threads.
  template <class Scalar, class Ordinal>
  struct DefaultKernels<Scalar, Ordinal, SerialNode> {
    typedef AltSparseOps<void, Ordinal, SerialNode,
                         details::AltSparseOpsDefaultAllocator<Ordinal, SerialNode> > SparseOps;
    typedef DefaultRelaxation<Scalar, Ordinal, SerialNode> Relaxations;
  };
#endif // defined(HAVE_TPETRACLASSIC_SERIAL)

#if defined(HAVE_TPETRACLASSIC_TBB)
  class TBBNode;
  //! Partial specialization for TBBNode, using first-touch allocation.
  template <class Scalar, class Ordinal>
  struct DefaultKernels<Scalar, Ordinal, TBBNode> {
    typedef DefaultHostSparseOps<void, Ordinal, TBBNode, details::FirstTouchCRSAllocator> SparseOps;
    typedef DefaultRelaxation<Scalar, Ordinal, TBBNode> Relaxations;
  };
#endif // HAVE_TPETRACLASSIC_TBB

#if defined(HAVE_TPETRACLASSIC_TPI)
  class TPINode;
  //! Partial specialization for TPINode, using first-touch allocation.
  template <class Scalar, class Ordinal>
  struct DefaultKernels<Scalar, Ordinal, TPINode> {
    typedef DefaultHostSparseOps<void, Ordinal, TPINode, details::FirstTouchCRSAllocator> SparseOps;
    typedef DefaultRelaxation<Scalar, Ordinal, TPINode> Relaxations;
  };
#endif // HAVE_TPETRACLASSIC_TPI

#if defined(HAVE_TPETRACLASSIC_OPENMP)
  class OpenMPNode;
  //! Partial specialization for OpenMPNode, using first-touch allocation.
  template <class Scalar, class Ordinal>
  struct DefaultKernels<Scalar, Ordinal, OpenMPNode> {
    typedef DefaultHostSparseOps<void, Ordinal, OpenMPNode, details::FirstTouchCRSAllocator> SparseOps;
    typedef DefaultRelaxation<Scalar, Ordinal, OpenMPNode> Relaxations;
  };
#endif // HAVE_TPETRACLASSIC_OPENMP

} // namespace KokkosClassic

#endif // KOKKOS_DEFAULT_KERNELS_HPP
