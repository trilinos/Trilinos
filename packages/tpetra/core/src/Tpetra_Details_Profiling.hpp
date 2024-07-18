// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_PROFILING_HPP
#define TPETRA_DETAILS_PROFILING_HPP

/// \file Tpetra_Details_Profiling.hpp
/// \brief Declaration of Tpetra::Details::Profiling, a scope guard
///   for Kokkos Profiling.

#include "TpetraCore_config.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_RCP.hpp"


namespace Tpetra {
namespace Details {

/// \brief Profile the given scope.
///
/// This class implements a "scope guard" for profiling a region of
/// code, using Kokkos Profiling.
///
/// If you disable Kokkos Profiling by setting some configuration
/// option such that KOKKOS_ENABLE_PROFILING is not defined, then this
/// class is harmless, but will do NOTHING.
///
/// Tpetra currently uses Kokkos::Profiling::pushRegion(name)
/// to start the region to profile, and Kokkos::Profiling::popRegion()
/// to close the region to profile.  The "scope guard" idiom means
/// that the destructor automatically closes the region, no matter how
/// the scope exits, even if an exception is thrown.
///
/// Here is an example of how to use this class:
/// \code
/// void someFunctionYouWrote () {
///   ProfilingRegion wholeFuncRegion ("someFunctionYouWrote");
///   doStuff ();
///   {
///     ProfilingRegion innerRegion1 ("inner-region-1");
///     doOtherStuff ();
///   }
///   doSomeMoreStuff ();
///   {
///     ProfilingRegion innerRegion2 ("yet-another-inner-region");
///     doYetMoreStuff ();
///   }
///   finallyDoEvenStillMoreStuff ();
/// }
/// \endcode
///
/// It's critical that you not forget to name the ProfilingRegion
/// variable.  If you forget, you will create a temporary
/// ProfilingRegion object that will first open, then immediately
/// close the "region" before anything else happens.  This is not what
/// you want.
///
/// You may name regions whatever you like.  They don't have to have
/// the same name as the region variable.  Kokkos Profiling has a
/// length limit of 512 characters for region names.
///
/// For details about Kokkos Profiling, see the kokkos-tools wiki:
///
/// https://github.com/kokkos/kokkos-tools/wiki
class ProfilingRegion {
public:
  //! Default constructor does not construct a region.
  ProfilingRegion();
  //! Open region to profile; name the region \c name.
  ProfilingRegion (const char name[]);
  //! Open region to profile, if the group name \c group is enabled by the
  //! TPETRA_TIMING variable; name the region \c name.
  ProfilingRegion (const char name[], const char group[]);
  //! Close region to profile.
  ~ProfilingRegion ();

private:
  bool kokkos_region_active_;
  Teuchos::RCP<Teuchos::TimeMonitor> tm;

};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PROFILING_HPP
