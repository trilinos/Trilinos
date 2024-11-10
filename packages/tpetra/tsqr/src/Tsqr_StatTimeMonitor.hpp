// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_StatTimeMonitor_hpp
#define __TSQR_StatTimeMonitor_hpp

#include "Teuchos_Time.hpp"
#include "Tsqr_TimeStats.hpp"

namespace TSQR {

  /// \class StatTimeMonitor
  /// \brief Like Teuchos::TimeMonitor, but collects running stats
  ///
  /// Like Teuchos::TimeMonitor, this class uses the RAII idiom to
  /// time a scope.  However, it also maintains running statistics,
  /// via a reference to a TimeStats object.
  ///
  /// \note Implementers: You may safely add new statistics to
  ///   TimeStats without needing to change this class.
  class StatTimeMonitor {
  public:
    /// \brief Constructor
    ///
    /// \param timer [in/out] Reference to the raw timer object
    ///
    /// \param stats [in/out] Running statistics, to be updated when
    ///   this StatTimeMonitor falls out of scope
    StatTimeMonitor (Teuchos::Time& timer, TimeStats& stats);

    /// \brief Destructor
    ///
    /// Updates running statistics via the TimeStats object.
    ~StatTimeMonitor ();

  private:
    // Copying syntactically forbidden
    StatTimeMonitor (const StatTimeMonitor&);
    StatTimeMonitor& operator= (const StatTimeMonitor&);

    Teuchos::Time& timer_;
    TimeStats& stats_;
  };

} // namespace TSQR

#endif // __TSQR_StatTimeMonitor_hpp
