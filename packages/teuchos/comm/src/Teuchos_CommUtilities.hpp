// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_COMM_UTILTIES_HPP
#define TEUCHOS_COMM_UTILTIES_HPP

#include "Teuchos_TimeMonitor.hpp"
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
#include "Kokkos_Profiling_ScopedRegion.hpp"
#endif

#ifdef HAVE_TEUCHOS_COMM_TIMERS

#define TEUCHOS_COMM_TIME_MONITOR( FUNCNAME ) \
  TEUCHOS_FUNC_TIME_MONITOR( FUNCNAME )

#else // HAVE_TEUCHOS_COMM_TIMERS

#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)

#define TEUCHOS_COMM_TIME_MONITOR( FUNCNAME ) \
  std::ostringstream TEUCHOS_COMM_TIME_MONITOR_oss; \
  TEUCHOS_COMM_TIME_MONITOR_oss << FUNCNAME; \
  ::Kokkos::Profiling::ScopedRegion TEUCHOS_COMM_TIME_MONITOR_scopedRegion( TEUCHOS_COMM_TIME_MONITOR_oss.str().c_str() )

#else

#define TEUCHOS_COMM_TIME_MONITOR( FUNCNAME )

#endif

#endif // HAVE_TEUCHOS_COMM_TIMERS

#endif // TEUCHOS_COMM_UTILTIES_HPP
