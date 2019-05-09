#ifndef TPETRA_DETAILS_BEHAVIOR_HPP
#define TPETRA_DETAILS_BEHAVIOR_HPP

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
///
/// These are two different things.  For example, TPETRA_DEBUG may do extra MPI
/// communication in order to ensure correct error state propagation, but
/// TPETRA_DEBUG should never print copious debug output if no errors occurred.
/// The idea is that if users get a mysterious error or hang, they can rerun
/// with TPETRA_DEBUG set.  TPETRA_VERBOSE is for Tpetra developers to use for
/// debugging Tpetra.
///
/// The environment variables are understood to be "on" or "off" and recognized
/// if specified in one of two ways.  The first is to specify the variable
/// unconditionally ON or OFF.  e.g., TPETRA_VERBOSE=ON or TPETRA_VERBOSE=OFF.
/// The default value of TPETRA_VERBOSE is always OFF.  The default value for
/// TPETRA_DEBUG is ON if Tpetra is configured with Tpetra_ENABLE_DEBUG,
/// otherwise it is OFF
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

  /// \brief Whether to assume that MPI is CUDA aware.
  ///
  /// An MPI implementation is "CUDA aware" if it can accept CUDA
  /// device buffers (Kokkos::CudaSpace) as send and receive buffers.
  /// You may control this behavior at run time via the
  /// <tt>TPETRA_ASSUME_CUDA_AWARE_MPI</tt> environment variable.
  ///
  /// For a discussion, see Trilinos GitHub issues #1571 and #1088.
  static bool assumeMpiIsCudaAware ();

  /// \brief The core count above which Tpetra::CrsMatrix::transferAndFillComplere
  /// will attempt to do advanced neighbor discovery. This is platform 
  /// dependent, and the user/developer should test each new platform 
  /// for the correct value. The 3000 value was found on Serrano summer 2018 

  static int TAFC_OptimizationCoreCount ();

};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_BEHAVIOR_HPP
