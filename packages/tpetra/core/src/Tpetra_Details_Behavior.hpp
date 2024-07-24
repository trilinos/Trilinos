// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_BEHAVIOR_HPP
#define TPETRA_DETAILS_BEHAVIOR_HPP

#include <stddef.h>

/// \file Tpetra_Details_Behavior.hpp
/// \brief Declaration of Tpetra::Details::Behavior, a class that
///   describes Tpetra's behavior.

namespace Tpetra {
namespace Details {

/// \brief Description of Tpetra's behavior.
///
/// "Behavior" means things like whether to do extra debug checks or
/// print debug output.  These depend both on build options and on
/// environment variables.  Build options generally control the
/// default behavior.
///
/// This class' methods have the following properties:
///
/// <ul>
/// <li>They may read the environment.</li>
/// <li>They read the environment at most once, on first call.</li>
/// <li>They are idempotent.</li>
/// <li>They are reentrant.  For example, they may use std::call_once
///     if appropriate.</li>
/// </ul>
///
/// We intended for it to be inexpensive to call this class' methods
/// repeatedly.  The idea is that you don't have to cache variables;
/// you should just call the functions freely.  In the common case,
/// the \c bool methods should just perform an 'if' test and just
/// return the \c bool value.  We spent some time thinking about how
/// to make the methods reentrant without a possibly expensive
/// mutex-like pthread_once / std::call_once cost on each call.
///
/// Tpetra does <i>not</i> promise to see changes to environment
/// variables made after using any Tpetra class or calling any Tpetra
/// function.  Best practice would be to set any environment variables
/// that you want to set, before starting the executable.
///
/// Our main goal with this class is to give both users and developers
/// more run-time control in determining Tpetra's behavior, by setting
/// environment variables.  This makes debugging much more efficient,
/// since before, enabling debugging code would have required
/// reconfiguring and recompiling.  Not all of Tpetra has bought into
/// this system yet; some debug code is still protected by macros like
/// <tt>HAVE_TPETRA_DEBUG</tt>.  However, our goal is that as much
/// Tpetra debugging code as possible can be enabled or disabled via
/// environment variable.  This will have the additional advantage of
/// avoiding errors due to only building and testing in debug or
/// release mode, but not both.
///
/// The behavior of Tpetra can be modified at runtime through two environment
/// variables:
///
/// TPETRA_DEBUG: flags Tpetra to turn on debug checking.
/// TPETRA_VERBOSE: flags Tpetra to turn on debug _output_.
/// TPETRA_TIMING: flags Tpetra to turn on timing code.
///
/// These are two different things.  For example, TPETRA_DEBUG may do extra MPI
/// communication in order to ensure correct error state propagation, but
/// TPETRA_DEBUG should never print copious debug output if no errors occurred.
/// The idea is that if users get a mysterious error or hang, they can rerun
/// with TPETRA_DEBUG set.  TPETRA_VERBOSE is for Tpetra developers to use for
/// debugging Tpetra. TPETRA_TIMING is for Tpetra developers to use for timing
/// Tpetra.
///
/// The environment variables are understood to be "on" or "off" and recognized
/// if specified in one of two ways. The first is to specify the variable
/// unconditionally ON or OFF. e.g., TPETRA_[VERBOSE,DEBUG,TIMING]=ON or
/// TPETRA_[VERBOSE,DEBUG,TIMING]=OFF. The default value of TPETRA_VERBOSE and
/// TPETRA_TIMING is always OFF. The default value for TPETRA_DEBUG is ON if
/// Tpetra is configured with Tpetra_ENABLE_DEBUG, otherwise it is OFF.
///
/// The second is to specify the variable on a per class/object basis, e.g.,
/// TPETRA_VERBOSE=CrsGraph,CrsMatrix,Distributor means that verbose output
/// will be enabled for CrsGraph, CrsMatrix, and Distributor classes.  For this
/// second method, the default values of both TPETRA_VERBOSE and TPETRA_DEBUG
/// is OFF.
class Behavior {
public:
  /// \brief Whether Tpetra is in debug mode.
  ///
  /// "Debug mode" means that Tpetra does extra error checks that may
  /// require more MPI communication or local computation.  It may
  /// also produce more detailed error messages, and more copious
  /// debug output.
  static bool debug ();

  /// \brief Whether the given Tpetra object is in debug mode.
  ///
  /// \param name [in] Name of the Tpetra object.  Typically, the object would
  ///        be a class name, e.g., "CrsGraph" or method, e.g.,
  ///        "CrsGraph::insertLocalIndices".
  static bool debug (const char name[]);

  /// \brief Whether Tpetra is in verbose mode.
  ///
  /// "Verbose mode" means that Tpetra prints copious debug output to
  /// std::cerr on every MPI process.  This is a LOT of output!  You
  /// really don't want to do this when running on many MPI processes.
  static bool verbose ();

  /// \brief Whether the given Tpetra object is in verbose mode.
  ///
  /// \param name [in] Name of the Tpetra object.  Typically, the object would
  ///        be a class name, e.g., "CrsGraph" or method, e.g.,
  ///        "CrsGraph::insertLocalIndices".
  static bool verbose (const char name[]);

  /// \brief Disable verbose mode, programatically
  static void disable_verbose_behavior ();

  /// \brief Enable verbose mode, programatically
  static void enable_verbose_behavior ();

  /// \brief Whether Tpetra is in timing mode.
  ///
  /// "Timing mode" means that Tpetra enables code that instruments internal timing.
  static bool timing ();

  /// \brief Whether the given Tpetra object is in timing mode.
  ///
  /// \param name [in] Name of the Tpetra object.  Typically, the object would
  ///        be a class name, e.g., "CrsGraph" or method, e.g.,
  ///        "CrsGraph::insertLocalIndices".
  static bool timing (const char name[]);

  /// \brief Disable timing, programatically
  static void disable_timing();

  /// \brief Enable timing, programatically
  static void enable_timing();

  /// \brief Whether to assume that MPI is CUDA aware.
  ///
  /// An MPI implementation is "CUDA aware" if it can accept CUDA
  /// device buffers (Kokkos::CudaSpace) as send and receive buffers.
  /// You may control this behavior at run time via the
  /// <tt>TPETRA_ASSUME_GPU_AWARE_MPI</tt> environment variable.
  ///
  /// For a discussion, see Trilinos GitHub issues #1571 and #1088.
  static bool assumeMpiIsGPUAware ();

  /// \brief Whether the CUDA_LAUNCH_BLOCKING environment variable has been set.
  static bool cudaLaunchBlocking ();

  /// \brief MPI process count above which
  ///   Tpetra::CrsMatrix::transferAndFillComplete will attempt to do
  ///   advanced neighbor discovery.
  ///
  /// This is platform dependent, and the user/developer should test
  /// each new platform for the correct value.  You may control this
  /// at run time via the <tt>MM_TAFC_OptimizationCoreCount</tt>
  /// environment variable.
  static int TAFC_OptimizationCoreCount ();

  /// \brief Number of entries below which arrays, lists, etc. will be
  ///   printed in debug mode.
  ///
  /// You may control this at run time via the
  /// <tt>TPETRA_VERBOSE_PRINT_COUNT_THRESHOLD</tt> environment
  /// variable.
  static size_t verbosePrintCountThreshold ();

  /// \brief Threshold for deciding if a local matrix is "imbalanced" in
  /// the number of entries per row. The threshold is compared against
  /// the difference between maximum row length and average row length.
  ///
  /// The threshold is measured in max number of entries in excess of the
  /// average (it is not a proportion between max and average).
  ///
  /// If the "imbalance" of a local matrix is greater than this threshold,
  /// a different algorithm may be used for some operations like 
  /// sparse matrix-vector multiply, packAndPrepare, and
  /// unpackAndCombine.  You may control this at run time via the
  /// <tt>TPETRA_ROW_IMBALANCE_THRESHOLD</tt> environment variable.
  static size_t rowImbalanceThreshold ();

  /// \brief Whether to use the cuSPARSE merge path algorithm to perform
  ///  sparse matrix-multivector products, one vector at a time. Depending on
  ///  the matrix and the number of vectors in the multivector, this may
  ///  be better than just applying the default SpMV algorithm to the entire
  ///  multivector at once.
  ///
  ///  Note: full support for merge path SPMV on multivectors
  ///  is coming soon.
  ///
  /// You may control this at run time via the
  /// <tt>TPETRA_MULTIVECTOR_USE_MERGE_PATH</tt> environment variable (default: false)
  static bool useMergePathMultiVector();

  /// \brief Unpack rows of a matrix using hierarchical unpacking
  static bool hierarchicalUnpack ();

  /// \brief Size of batch for hierarchical unpacking
  static size_t hierarchicalUnpackBatchSize ();

  /// \brief Size of team for hierarchical unpacking
  static size_t hierarchicalUnpackTeamSize ();

  /// \brief the threshold for transitioning from device to host
  ///
  /// If the number of elements in the multivector does not exceed this 
  /// threshold and the data is on host, then run the calculation on
  /// host.  Otherwise, run on device.
  /// By default this is 10000, but may be altered by the environment
  /// variable TPETRA_VECTOR_DEVICE_THRESHOLD
  static size_t multivectorKernelLocationThreshold ();

  /// \brief Use Teuchos::Timer in Tpetra::ProfilingRegion
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_USE_TEUCHOS_TIMERS</tt> environment variable.
  static bool profilingRegionUseTeuchosTimers();

  /// \brief Use Kokkos::Profiling in Tpetra::ProfilingRegion
  ///
  /// This is enabled by default if KOKKOS_ENABLE_PROFILING is defined.
  /// You may control this at run time via the <tt>TPETRA_USE_KOKKOS_PROFILING</tt>
  /// environment variable.
  static bool profilingRegionUseKokkosProfiling();

  /// \brief Fusing SpMV and update in residual instead of using 2 kernel launches.
  /// Fusing kernels implies that no TPLs (CUSPARSE, ROCSPARSE, ...) will be used for the residual.
  ///
  /// This is enabled by default.  You may control this at run time via the
  /// <tt>TPETRA_FUSED_RESIDUAL</tt> environment variable.
  static bool fusedResidual();

  /// \brief Skip copyAndPermute if possible
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_SKIP_COPY_AND_PERMUTE</tt> environment variable.
  static bool skipCopyAndPermuteIfPossible();

  /// \brief Overlap communication and computation.
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_OVERLAP</tt> environment variable.
  static bool overlapCommunicationAndComputation();

  /// \brief Add Teuchos timers for all host calls to Kokkos::deep_copy().
  /// This is especially useful for identifying host/device data transfers
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_TIME_KOKKOS_DEEP_COPY</tt> environment variable.
  static bool timeKokkosDeepCopy();
  
  /// \brief Adds verbose output to Kokkos deep_copy timers
  /// by appending source and destination.
  /// This is especially useful for identifying host/device data transfers
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_TIME_KOKKOS_DEEP_COPY_VERBOSE1</tt> environment variable.
  static bool timeKokkosDeepCopyVerbose1();

  
  /// \brief Adds verbose output to Kokkos deep_copy timers
  /// by appending source, destination, and size.
  /// This is especially useful for identifying host/device data transfers
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_TIME_KOKKOS_DEEP_COPY_VERBOSE2</tt> environment variable.
  static bool timeKokkosDeepCopyVerbose2();

  /// \brief Add Teuchos timers for all host calls to Kokkos::fence().
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_TIME_KOKKOS_FENCE</tt> environment variable.
  static bool timeKokkosFence();  

  /// \brief Add Teuchos timers for all host calls to Kokkos::parallel_for(), 
  /// Kokkos::parallel_reduce() and Kokkos::parallel_scan().
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_TIME_KOKKOS_FUNCTIONS</tt> environment variable.
  static bool timeKokkosFunctions();  

  /// \brief Warn if more than this many Kokkos spaces are accessed.
  ///
  /// This is disabled by default.  You may control this at run time via the
  /// <tt>TPETRA_SPACE_ID_WARN_LIMIT</tt> environment variable.
  static size_t spacesIdWarnLimit();

  /// \brief Search the environment for TPETRA_ variables and reject unrecognized ones
  static void reject_unrecognized_env_vars();
};



} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_BEHAVIOR_HPP
