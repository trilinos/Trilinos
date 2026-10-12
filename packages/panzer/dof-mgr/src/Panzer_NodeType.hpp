// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_KOKKOS_NODE_TYPE_HPP
#define PANZER_KOKKOS_NODE_TYPE_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <type_traits>

namespace panzer {

  /**
   * \brief The Kokkos node type used to instantiate all Tpetra objects (Map, MultiVector, CrsGraph, etc.) throughout panzer.
   *
   * Tpetra classes are templated on a Kokkos node type in addition to
   * their scalar/ordinal types. This typedef follows Phalanx: it pairs
   * `PHX::ExecutionSpace` with `PHX::MemorySpace`, so that Tpetra
   * objects created by panzer run where Phalanx evaluators run and
   * allocate where Phalanx fields allocate. Use this typedef wherever
   * a Tpetra node type template argument is needed, e.g.
   * `Tpetra::Map<LocalOrdinal,GlobalOrdinal,TpetraNodeType>`, rather
   * than naming a node or a space directly.
   *
   * Tpetra must have been configured to instantiate this node type,
   * which is what the static_assert below checks. The two memory space
   * choices have to agree as well: `Phalanx_ENABLE_SHARED_SPACE` and
   * `Tpetra_ALLOCATE_IN_SHARED_SPACE` both select the same
   * `Kokkos::SharedSpace`, so enabling one without the other produces a
   * node type Tpetra did not instantiate.
   */
  using TpetraNodeType =
      Tpetra::KokkosCompat::KokkosDeviceWrapperNode<PHX::ExecutionSpace,
                                                    PHX::MemorySpace>;

  /** \brief Whether Tpetra explicitly instantiated its classes over \c Node.
   *
   * Specialized below for each node type Tpetra was configured to
   * instantiate. Tpetra instantiates exactly one execution space by
   * default, and HIP and SYCL are never the default, so a GPU build
   * that does not set the matching Tpetra_INST_<BACKEND> gets a host
   * node while Phalanx picks the GPU.
   *
   * Comparing node types rather than execution spaces is deliberate:
   * Tpetra's node typedefs already carry the memory space, including
   * its Tpetra_ALLOCATE_IN_SHARED_SPACE choice, so one comparison
   * covers the execution space and both shared space switches at once.
   * It also compares types rather than spellings, so it cannot be
   * fooled by one type named two ways.
   */
  template <typename Node>
  struct TpetraNodeIsInstantiated : std::false_type {};

#if defined(HAVE_TPETRA_INST_SERIAL) && defined(KOKKOS_ENABLE_SERIAL)
  template <>
  struct TpetraNodeIsInstantiated<Tpetra::KokkosCompat::KokkosSerialWrapperNode>
      : std::true_type {};
#endif
#if defined(HAVE_TPETRA_INST_OPENMP) && defined(KOKKOS_ENABLE_OPENMP)
  template <>
  struct TpetraNodeIsInstantiated<Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>
      : std::true_type {};
#endif
#if defined(HAVE_TPETRA_INST_PTHREAD) && defined(KOKKOS_ENABLE_THREADS)
  template <>
  struct TpetraNodeIsInstantiated<Tpetra::KokkosCompat::KokkosThreadsWrapperNode>
      : std::true_type {};
#endif
#if defined(HAVE_TPETRA_INST_CUDA) && defined(KOKKOS_ENABLE_CUDA)
  template <>
  struct TpetraNodeIsInstantiated<Tpetra::KokkosCompat::KokkosCudaWrapperNode>
      : std::true_type {};
#endif
#if defined(HAVE_TPETRA_INST_HIP) && defined(KOKKOS_ENABLE_HIP)
  template <>
  struct TpetraNodeIsInstantiated<Tpetra::KokkosCompat::KokkosHIPWrapperNode>
      : std::true_type {};
#endif
#if defined(HAVE_TPETRA_INST_SYCL) && defined(KOKKOS_ENABLE_SYCL)
  template <>
  struct TpetraNodeIsInstantiated<Tpetra::KokkosCompat::KokkosSYCLWrapperNode>
      : std::true_type {};
#endif

  static_assert(TpetraNodeIsInstantiated<TpetraNodeType>::value,
                "panzer: the Tpetra node type panzer derives from Phalanx' "
                "execution and memory spaces is not one that Tpetra "
                "instantiated, so every Tpetra object panzer creates would "
                "fail to link.  Either set the Tpetra_INST_<BACKEND> matching "
                "Phalanx_DEFAULT_EXECUTION_SPACE -- note Tpetra_INST_HIP and "
                "Tpetra_INST_SYCL are OFF by default, so a HIP or SYCL build "
                "needs one set explicitly -- or make "
                "Tpetra_ALLOCATE_IN_SHARED_SPACE agree with "
                "Phalanx_ENABLE_SHARED_SPACE, since the two select the same "
                "Kokkos::SharedSpace and enabling one alone changes the node's "
                "memory space to one Tpetra did not instantiate.");

}

#endif
