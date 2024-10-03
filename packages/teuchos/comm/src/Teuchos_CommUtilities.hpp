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

#ifdef HAVE_TEUCHOS_COMM_TIMERS

#define TEUCHOS_COMM_TIME_MONITOR( FUNCNAME ) \
  TEUCHOS_FUNC_TIME_MONITOR( FUNCNAME )

#else // HAVE_TEUCHOS_COMM_TIMERS

#define TEUCHOS_COMM_TIME_MONITOR( FUNCNAME )

#endif // HAVE_TEUCHOS_COMM_TIMERS

#endif // TEUCHOS_COMM_UTILTIES_HPP
