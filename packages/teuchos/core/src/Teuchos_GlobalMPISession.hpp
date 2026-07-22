// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_GLOBAL_MPI_SESSION_HPP
#define TEUCHOS_GLOBAL_MPI_SESSION_HPP

/*! \file Teuchos_GlobalMPISession.hpp
    \brief A MPI utilities class, providing methods for initializing,
        finalizing, and querying the global MPI session
*/

#include "TeuchosCore_ConfigDefs.hpp"

#include "Teuchos_ArrayView.hpp"


namespace Teuchos {

/// \class GlobalMPISession
/// \brief Initialize, finalize, and query the global MPI session.
///
/// This class insulates basic <tt>main()</tt> program type of code from
/// having to know if MPI is enabled or not.  The typical use case is to
/// replace an explicit call to MPI_Init() in your main() routine with
/// creation of a GlobalMPISession instance.  The instance's destructor (which
/// in this case will be called at the end of main()) then calls
/// MPI_Finalize().  So, instead of writing:
//
/// \code
/// int main () {
///   (void) MPI_Init (&argc, &argv);
///   // Your code goes here ...
///   (void) MPI_Finalize ();
///   return 0;
/// }
/// \endcode
/// you would write:
/// \code
/// #include <Teuchos_GlobalMPISession.hpp>
///
/// int main () {
///   Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
///   // Your code goes here ...
///   return 0;
/// }
///
/// \endcode
///
/// This saves you from needing to remember to call MPI_Init() or
/// MPI_Finalize().  Also, having the GlobalMPISession object's constructor
/// call MPI_Finalize() allows destructors from other objects to call MPI
/// functions.  That wold never be possible if you were to directly call
/// MPI_Finalize() at the end of main().
///
/// This class even works if you have not built Teuchos with MPI support.  In
/// that case, it behaves as if MPI_COMM_WORLD had one process, which is
/// always the calling process.  Thus, you can use this class to insulate your
/// code from needing to know about MPI.  You don't even have to include
/// mpi.h, as long as your code doesn't directly use MPI routines or types.
/// Teuchos implements wrappers for MPI communicators (see the Teuchos::Comm
/// class and its subclasses in the TeuchosComm subpackage) which allow you to
/// use a very very small subset of MPI functionality without needing to
/// include mpi.h or depend on MPI in any way.
///
/// This class also contains the most minimal of other static member functions
/// that are needed for only the most simplistic of tasks needed by other
/// TeuchosCore software.  For example, you can do a barrier or sum an int
/// across processes.  These are needed by the most basic operations involving
/// output or determining success or failure across processes for unit tests.
///
/// GlobalMPISession's static functions cleverly checks whether MPI has been
/// initialized already before calling any MPI functions.  Therefore, you can
/// use it in your libraries without requiring that a GlobalMPISession object
/// was created in main().
class TEUCHOSCORE_LIB_DLL_EXPORT GlobalMPISession
{
public:

  //! @name Public constructor and destructor
  //@{

  /** \brief Calls <tt>MPI_Init()</tt> if MPI is enabled.
   *
   * \param argc [in] Address of the argument passed into
   * <tt>main(argc,argv)</tt>.  Same as the first argument of MPI_Init().
   *
   * \param argv [in] Address of the argument passed into
   * <tt>main(argc,argv)</tt>.  Same as the second argument of MPI_Init().
   *
   * \param out [in] If <tt> out != NULL</tt>, then a small message on will be
   * printed to this stream on <i>each</i> process in <tt>MPI_COMM_WORLD</tt>.
   * The default is <tt>&std::cout</tt>.
   *
   * If the command-line arguments include the option
   * <tt>--teuchos-suppress-startup-banner</tt>, the this option will be
   * removed from <tt>argv[]</tt> before being passed to
   * <tt>MPI_Init(...)</tt>, and the startup output message to <tt>*out</tt>
   * will be suppressed.
   *
   * If Teuchos was <i>not</i> built with MPI support, the constructor
   * just prints a startup banner (unless the banner was suppressed --
   * see previous paragraph).
   *
   * \warning The default third parameter may result in a lot of lines printed
   * to std::cout, if <tt>MPI_COMM_WORLD</tt> is large!  Users should
   * generally pass in <tt>NULL</tt> for the third argument.  On the other
   * hand, it can be useful to see that startup banner, just to know that MPI
   * is working.
   *
   * \warning If MPI_Initialized() returns true before calling this
   * constructor, then a error message is printed to <tt>*out</tt>
   * std::terminate() is called.  This is the only sane behavior for this
   * constuctor.  The very nature of the GlboalMPISession object is to be
   * constructed at the tope of main() outside of a try {} block.  If MPI
   * can't be initialized, then the only thing to do is to abort the program.
   * It would not be reasonble to simply not not call MPI_Initialize() becuase
   * this would ignore the input arguments that the user (you) should be
   * expecting to be read.
   *
   * \warning Any other MPI functions called direclty from within this
   * constructor will result in an error message to be printed and for abort
   * to be called.
   */
  GlobalMPISession( int* argc, char*** argv, std::ostream *out = &std::cout );

  //! Call <tt>MPI_Finalize()</tt> if MPI is enabled.
  ~GlobalMPISession();

  //@}

  //! @name Static functions
  //@{

  /// \brief abort the program
  ///
  /// Calls MPI_Abort for HAVE_MPI
  /// Otherwise calls std::abort
  static void abort();

  //! @name Static functions
  //@{

  /// \brief Return whether MPI was initialized.
  ///
  /// This is always true if the constructor returned.  If the
  /// constructor was not called, it may or may not be true, depending
  /// on whether the user called MPI_Init() themselves.  If the
  /// constructor was called but threw an exception, then some MPI
  /// function returned an error code.
  static bool mpiIsInitialized();

  /// \brief Return whether MPI was already finalized.
  ///
  /// This is always true if the destructor was called.  If the
  /// destructor was not called, it may or may not be true, depending
  /// on whether the user called MPI_Init() themselves.
  static bool mpiIsFinalized();

  /** \brief The rank of the calling process in <tt>MPI_COMM_WORLD</tt>.
   *
   * \return <tt>0</tt> if MPI has not yet been initialized, else the
   *   rank of the calling process in <tt>MPI_COMM_WORLD</tt>.
   *
   * You may call this method even if the constructor was never
   * called.  Thus, it is safe to use no matter how MPI_Init() was
   * called.  However, MPI_Init() must have been called somehow in
   * order for this method to return a sensible result.
   */
  static int getRank();

  /** \brief The number of processes in <tt>MPI_COMM_WORLD</tt>.
   *
   * \return <tt>1</tt> if MPI has not yet been initialized, else the
   *   number of processes in <tt>MPI_COMM_WORLD</tt>.
   *
   * You may call this method even if the constructor was never
   * called.  Thus, it is safe to use no matter how MPI_Init() was
   * called.  However, MPI_Init() must have been called somehow in
   * order for this method to return a sensible result.
   */
  static int getNProc();

  /// \brief Call MPI_Barrier() on <tt>MPI_COMM_WORLD</tt>.
  ///
  /// This method must be called collectively on all processes in
  /// <tt>MPI_COMM_WORLD</tt>.
  ///
  /// \note Users should invoke barrier through the Teuchos::Comm
  ///   interface.  We only expose this method for Teuchos-internal
  ///   functionality.
  static void barrier();

  /** \brief Sum a set of integers across processes.
   *
   * This performs an MPI_Allreduce() of localVal over
   * <tt>MPI_COMM_WORLD</tt>, and returns the result (which is the
   * same on all processes).
   *
   * This method must be called collectively on all processes in
   * <tt>MPI_COMM_WORLD</tt>.
   *
   * \param localVal [in] Value on local process to sum across processes.
   * \return The global sum (on all processes).
   *
   * \note Users should invoke reductions through the Teuchos::Comm
   *   interface.  We only expose this method for Teuchos-internal
   *   functionality.
   */
  static int sum(int localVal);

  /** \brief Global all-to-all of a set of integers across processes.
   *
   * This performs an MPI_Allgather() of localVal over
   * <tt>MPI_COMM_WORLD</tt>, and writes the results (which are the
   * same on all processes) to allVals.
   *
   * This method must be called collectively on all processes in
   * <tt>MPI_COMM_WORLD</tt>.
   *
   * \param localVal [in] Value on local process to pass to all processes.
   *
   * \param allVals [out] Array (length getNProc()) that gives the
   *   value of localVal for each process.  On output, allVals[k] is
   *   localVal for process k.
   */
  static void allGather(int localVal, const ArrayView<int> &allVals);

#ifdef HAVE_TEUCHOSCORE_KOKKOS
  /// \brief Fetch a deep copy of the input arguments to \c main()
  ///   (that is, \c argv), as given to GlobalMPISession's
  ///   constructor.
  ///
  /// \return If GlobalMPISession's constructor hasn't been called
  ///   yet, return an empty vector.  Else, return the input
  ///   arguments.
  static std::vector<std::string> getArgv ();
#endif // HAVE_TEUCHOSCORE_KOKKOS
  //@}

private:

  static bool haveMPIState_;
  static bool mpiIsFinalized_;
  static int rank_;
  static int nProc_;
#ifdef HAVE_TEUCHOSCORE_KOKKOS
  /// \brief Deep copy of the input arguments.
  ///
  /// This is useful if we want to call Kokkos::initialize later with
  /// the correct command-line arguments.  We keep a deep copy because
  /// applications may choose to modify the command-line arguments
  /// after calling this object's constructor.  That could mess up
  /// indexing if we just keep a pointer to the original.
  static std::vector<std::string> argvCopy_;
#endif // HAVE_TEUCHOSCORE_KOKKOS

  static void initialize( std::ostream *out );

  static void justInTimeInitialize();

};

} // namespace Teuchos

#endif // TEUCHOS_GLOBAL_MPI_SESSION_HPP
