#ifndef __TACHO_TEST_HPP__
#define __TACHO_TEST_HPP__

#if defined(TACHO_USE_DEPRECATED_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::DeprecatedTaskScheduler<T>;
static const char * scheduler_name = "DeprecatedTaskScheduler";
#endif
#if defined(TACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE)
template<typename T> using TaskSchedulerType = Kokkos::DeprecatedTaskSchedulerMultiple<T>;
static const char * scheduler_name = "DeprecatedTaskSchedulerMultiple";
#endif
#if defined(TACHO_USE_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::TaskScheduler<T>;
static const char * scheduler_name = "TaskScheduler";
#endif
#if defined(TACHO_USE_TASKSCHEDULER_MULTIPLE)
template<typename T> using TaskSchedulerType = Kokkos::TaskSchedulerMultiple<T>;
static const char * scheduler_name = "TaskSchedulerMultiple";
#endif
#if defined(TACHO_USE_CHASELEV_TASKSCHEDULER)
template<typename T> using TaskSchedulerType = Kokkos::ChaseLevTaskScheduler<T>;
static const char * scheduler_name = "ChaseLevTaskScheduler";
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
