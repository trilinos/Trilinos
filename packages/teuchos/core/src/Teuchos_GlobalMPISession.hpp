// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_GLOBAL_MPI_SESSION_HPP
#define TEUCHOS_GLOBAL_MPI_SESSION_HPP

/*! \file Teuchos_GlobalMPISession.hpp
    \brief A MPI utilities class, providing methods for initializing,
        finalizing, and querying the global MPI session
*/
#include "Teuchos_ArrayView.hpp"


namespace Teuchos {

/// \class GlobalMPISession
/// \brief Initialize, finalize, and query the global MPI session.
///
/// This class insulates basic <tt>main()</tt> program type of code
/// from having to know if MPI is enabled or not.  The typical use
/// case is to replace an explicit call to MPI_Init() in your main()
/// routine with creation of a GlobalMPISession instance.  The
/// instance's destructor (which in this case will be called at the
/// end of main()) then calls MPI_Finalize().  So, instead of writing:
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
///   Teuchos::GlobalMPISession (&argc, &argv, NULL);
///   // Your code goes here ...
///   return 0;
/// }
/// \endcode
/// This saves you from needing to remember to call MPI_Init() or
/// MPI_Finalize().  GlobalMPISession cleverly checks whether MPI has
/// been initialized already before calling MPI_Init(), so you can use
/// it in your libraries without needing to know whether users have
/// called MPI_Init() yet.
///
/// This class even works if you have not built Trilinos with MPI support.  In
/// that case, it behaves as if MPI_COMM_WORLD had one process, which is
/// always the calling process.  Thus, you can use this class to insulate your
/// code from needing to know about MPI.  You don't even have to include
/// mpi.h, as long as your code doesn't directly use MPI routines or types.
/// Teuchos implements wrappers for MPI communicators (see the Comm class and
/// its subclasses in the TeuchosComm subpackage) which allow you to use a
/// very very small subset of MPI functionality without needing to include
/// mpi.h or depend on MPI in any way.
///
/// This class also contains the the most minimal of other members that are
/// needed for only the most simplistic of tasks needed by other TeuchosCore
/// software.  For example, you can do a barrier or sum an int across
/// processes.  These are needed by the most basic operations involving output
/// or determiing success or failure across processes for unit tests.
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
   * \param out [in] If <tt>out!=NULL</tt>, then a small message on each
   * processor will be printed to this stream.  The default is
   * <tt>&std::cout</tt>.
   *
   * If the option <tt>--teuchos-suppress-startup-banner</tt> is found, the
   * this option will be removed from <tt>argv[]</tt> before being passed to
   * <tt>MPI_Init(...)</tt> and the startup output message to <tt>*out</tt>
   * will be suppressed.
   *
   * If Teuchos was <i>not</i> built with MPI support, the constructor
   * just prints a startup banner (unless the banner was suppressed --
   * see previous paragraph).
   *
   * \warning This constructor can only be called once per executable.
   * Otherwise, an error is printed to <tt>*out</tt> and an std::exception
   * will be thrown!
   */
  GlobalMPISession( int* argc, char*** argv, std::ostream *out = &std::cout );

  //! Call <tt>MPI_Finalize()</tt> if MPI is enabled.
  ~GlobalMPISession();

  //@}

  //! @name Static functions
  //@{

  //! Return whether MPI was initialized.
  static bool mpiIsInitialized();

  //! Return whether MPI was already finalized.
  static bool mpiIsFinalized();

  /** \brief Returns the process rank relative to <tt>MPI_COMM_WORLD</tt>
   *
   * Returns <tt>0</tt> if MPI is not enabled.
   *
   * Note, this function can be called even if the above constructor was never
   * called so it is safe to use no matter how <tt>MPI_Init()</tt> got called
   * (but it must have been called somewhere).
   */
  static int getRank();

  /** \brief Returns the number of processors relative to
   * <tt>MPI_COMM_WORLD</tt>
   *
   * Returns <tt>1</tt> if MPI is not enabled.
   *
   * Note, this function can be called even if the above constructor was never
   * called so it is safe to use no matter how <tt>MPI_Init()</tt> got called
   * (but it must have been called somewhere).
   */
  static int getNProc();

  /** \brief Perform a barrier so that all processes line up.
   */
  static void barrier();

  /** \brief Sum a set of integers across processes.
   *
   * \param localVal [in] Value on local process to sum across processes.
   */
  static int sum(int localVal);

  /** \brief Global all-to-all of a set of integers across processes.
   *
   * \param localVal [in] Value on local process to pass to all processes.
   *
   * \param allVals [out] Array (length getNProc()) that gives the value of
   * localVal for each process.  On output, allVals[k] is localVal for process
   * k.
   */
  static void allGather(int localVal, const ArrayView<int> &allVals);

  //@}

private:

  static bool haveMPIState_;
  static bool mpiIsFinalized_;
  static int rank_;
  static int nProc_;

  static void initialize( std::ostream *out );

  static void justInTimeInitialize();

};

} // namespace Teuchos

#endif // TEUCHOS_GLOBAL_MPI_SESSION_HPP
