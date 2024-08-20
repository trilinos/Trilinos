// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_KOKKOS_TEUCHOS_TIMER_INJECTION_HPP
#define TPETRA_DETAILS_KOKKOS_TEUCHOS_TIMER_INJECTION_HPP

/// \file Tpetra_Details_KokkosTeuchosTimerInjection.hpp
/// \brief Declaration functions that use Kokkos' profiling library to add deep copies between memory spaces, 
/// and Kokkos fences to the Teuchos::TimeMonitor system. This does have the side effect of making Kokkos::deep_copy()
///  calls on the host also call Kokkos::fence()



namespace Tpetra {
namespace Details {

  // The force option overrides the environment variable control via TPETRA_TIME_KOKKOS_DEEP_COPY
  // This is used for unit testing the capability
  void AddKokkosDeepCopyToTimeMonitor(bool force = false);

  // The force option overrides the environment variable control via TPETRA_TIME_KOKKOS_FENCE
  // This is used for unit testing the capability
  void AddKokkosFenceToTimeMonitor(bool force = false);

  // The force option overrides the environment variable control via TPETRA_TIME_KOKKOS_FUNCTIONS
  // This is used for unit testing the capability
  void AddKokkosFunctionsToTimeMonitor(bool force = false);


} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_KOKKOS_TEUCHOS_TIMER_INJECTION_HPP
