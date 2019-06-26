#ifndef __TACHO_TEST_HPP__
#define __TACHO_TEST_HPP__

#if defined(TACHO_USE_DEPRECATED_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::DeprecatedTaskScheduler<T>;
#endif
#if defined(TACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE)
template<typename T> using TaskSchedulerType = Kokkos::DeprecatedTaskSchedulerMultiple<T>;
#endif
#if defined(TACHO_USE_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::TaskScheduler<T>;
#endif
#if defined(TACHO_USE_TASKSCHEDULER_MULTIPLE)
template<typename T> using TaskSchedulerType = Kokkos::TaskSchedulerMultiple<T>;
#endif
#if defined(TACHO_USE_CHASELEV_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::ChaseLevTaskScheduler<T>;
#endif


#include "Tacho_TestCrsMatrixBase.hpp"
#include "Tacho_TestGraph.hpp"
#include "Tacho_TestSymbolic.hpp"
#include "Tacho_TestNumeric.hpp"
//#include "Tacho_TestTaskFunctor.hpp"

#include "Tacho_TestDenseMatrixView.hpp"
#include "Tacho_TestDenseByBlocks.hpp"

#include "Tacho_TestDenseLinearAlgebra.hpp"

#endif
