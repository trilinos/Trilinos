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
class Behavior {
public:
  /// \brief Whether Tpetra is in debug mode.
  ///
  /// "Debug mode" means that Tpetra does extra error checks that may
  /// require more MPI communication or local computation.  It may
  /// also produce more detailed error messages, and more copious
  /// debug output.
  static bool debug ();

  /// \brief Whether Tpetra is in verbose mode.
  ///
  /// "Verbose mode" means that Tpetra prints copious debug output to
  /// std::cerr on every MPI process.  This is a LOT of output!  You
  /// really don't want to do this when running on many MPI processes.
  static bool verbose ();
};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_BEHAVIOR_HPP
