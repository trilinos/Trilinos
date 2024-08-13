// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKHELPER_TIMERS_HPP
#define IFPACK2_BLOCKHELPER_TIMERS_HPP


namespace Ifpack2 {

  namespace BlockHelperDetails {

#if defined(HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS)
#define IFPACK2_BLOCKHELPER_TIMER(label) TEUCHOS_FUNC_TIME_MONITOR(label);
#define IFPACK2_BLOCKHELPER_TIMER_FENCE(execution_space) execution_space().fence();
#define IFPACK2_BLOCKHELPER_TIMER_DEFAULT_FENCE() Kokkos::DefaultExecutionSpace().fence();
#else
#define IFPACK2_BLOCKHELPER_TIMER(label)
#define IFPACK2_BLOCKHELPER_TIMER_FENCE(execution_space)
#define IFPACK2_BLOCKHELPER_TIMER_DEFAULT_FENCE()
#endif

#define IFPACK2_BLOCKHELPER_TIMER_WITH_FENCE(label, execution_space) \
  IFPACK2_BLOCKHELPER_TIMER_FENCE(execution_space) \
  IFPACK2_BLOCKHELPER_TIMER(label)

  } // namespace BlockHelperDetails

} // namespace Ifpack2

#endif
