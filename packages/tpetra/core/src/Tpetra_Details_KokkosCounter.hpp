// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// clang-format off
#ifndef TPETRA_DETAILS_KOKKOS_COUNTER_HPP
#define TPETRA_DETAILS_KOKKOS_COUNTER_HPP

/// \file Tpetra_Details_KokkosCounter.hpp
/// \brief Declaration of various tools for counting Kokkos calls of various
/// types using the Kokkos Profiling Library

#include <string>
#include <Teuchos_FancyOStream.hpp>

namespace Tpetra {
namespace Details {

/// \brief Counter for Kokkos::deep_copy calls
namespace DeepCopyCounter {
  /// \brief Start the deep_copy counter
  void start();

  /// \brief Reset the deep_copy counter
  void reset();

  /// \brief Stop the deep_copy counter
  void stop();

  /// \brief Query the deep_copy counter for copies in the same space
  size_t get_count_same_space();

  /// \brief Query the deep_copy counter for copies between different spaces
  size_t get_count_different_space();

}

/// \brief Counter for Kokkos::fence calls
namespace FenceCounter {
  /// \brief Start the fence counter
  void start();

  /// \brief Reset the fence counter
  void reset();

  /// \brief Stop the fence counter
  void stop();

  /// \brief Query the fence counter for given device, for an exec_space_instance.fence()
  size_t get_count_instance(const std::string & device);

  /// \brief Query the fence counter for given device, for an Kokkos::fence()
  size_t get_count_global(const std::string & device);
}

// clang-format on

/// \brief Counter for Kokkos regions representing third-party library usage
namespace KokkosRegionCounter {
/// \brief Start the counter
void start();

/// \brief Reset the counter
void reset();

/// \brief Stop the counter
void stop();

/// \brief How many regions containing `substr` have been seen
size_t get_count_region_contains(const std::string &substr);

/// \brief Print all observed region labels, separated by newline
void dump_regions(std::ostream &os);
void dump_regions(Teuchos::FancyOStream &os);
} // namespace KokkosRegionCounter

// clang-format off



} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_KOKKOS_COUNTER_HPP
