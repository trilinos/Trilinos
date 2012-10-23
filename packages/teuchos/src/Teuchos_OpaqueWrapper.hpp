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

#ifndef TEUCHOS_OPAQUE_WRAPPER_HPP
#define TEUCHOS_OPAQUE_WRAPPER_HPP


#include "Teuchos_RCP.hpp"


//#define TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE

#ifdef TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE
#  include "Teuchos_VerboseObject.hpp"
#endif // TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE


namespace Teuchos {


/** \class OpaqueWrapper
 * \brief Base class for wrapped opaque objects.
 * \tparam Opaque Type of the opaque object (a.k.a. handle).
 *
 * \subsection Teuchos_OpaqueWrapper_Summary Summary
 *
 * If you want to <i>create</i> an RCP to an opaque object, use the
 * opaqueWrapper() nonmember template function.  The <i>type</i> of an
 * RCP to an opaque object T is <tt>RCP<OpaqueWrapper<T> ></tt>.
 *
 * \section Teuchos_OpaqueWrapper_Prereq Prerequisites
 *
 * In order to understand this documentation, you must first have
 * learned how to use RCP (Teuchos' reference-counted pointer class)
 * to manage dynamically allocated memory and other resources.  It
 * also helps to be familiar with MPI (the Message Passing Interface
 * for distributed-memory parallel programming), but this is not
 * required.
 *
 * \section Teuchos_OpaqueWrapper_Handles What are opaque objects (a.k.a. handles)?
 *
 * Many different software libraries use the <i>opaque object</i> or
 * opaque handle idiom to hide the internals of a data structure from
 * users.  This standard technique allows users to treat an instance
 * of a data structure as a handle.  Users may pass the handle around
 * as if it were a simple value type (like int), and must call
 * nonmember functions in order to create, operate on, use, or destroy
 * instances of the data structure.  The MPI (Message Passing
 * Interface) standard is full of examples of the opaque object idiom,
 * including MPI_Comm (for communicators), MPI_Datatype (for standard
 * and custom data types), and MPI_Op (for standard and custom
 * reduction operations).
 *
 * In general, opaque handles (corresponding to the Opaque template
 * parameter) must be <i>assignable</i>.  This means that copy
 * construction and assignment (operator=) must be syntactically
 * correct for instances of Opaque.  This is certainly true of MPI's
 * opaque handle types.
 *
 * Opaque handles are a useful technique, but they interfere with
 * correct use of RCP, as we will explain below.  This base class
 * allows opaque objects to be wrapped by a real object, whose address
 * you can take.  This is needed in order to wrap an opaque object in
 * a RCP, for example.
 *
 * \section Teuchos_OpaqueWrapper_Special Why do opaque handles need special treatment?
 *
 * This class was motivated in particular by MPI's common use of the
 * opaque object idiom.  For MPI, passing MPI_Comm, MPI_Datatype, and
 * MPI_Op objects around by handles hides implementation details from
 * the user.  Handles also make it easier to access MPI functionality
 * from Fortran, so that C, C++, and Fortran can all share the same
 * handle mechanism.  In fact, some MPI implementations (such as
 * MPICH, at least historically if not currently) simply implement
 * these handles all as integers.  (As the MPI standard's advice to
 * implementers suggests, such an implementation would likely maintain
 * a table for each MPI process that maps the integer value to a
 * pointer to the corresponding object.)  For example, MPI_Comm might
 * be a typedef to int, and MPI_COMM_WORLD might be a literal integer
 * value:
 *
 \code
 typedef int MPI_Comm;
 #define MPI_COMM_WORLD 42
 \endcode

 * In this case, the expression <tt>rcp(&MPI_COMM_WORLD)</tt> would
 * not even compile, since one cannot take the address of an integer
 * literal.  To make this this compile, one might try the following:

 \code

 // THIS FUNCTION IS WRONG.  IT MAY SEGFAULT.
 Teuchos::RCP<MPI_Comm> getMpiCommPtr()
 {
   MPI_Comm comm = MPI_COMM_WORLD;
   // WRONG!!!  comm is a stack variable!
   return Teuchos::rcp (&comm, false);
 }

 \endcode

 * Using the returned communicator would result in undefined behavior,
 * which in practice might be a segfault, memory corruption, or MPI
 * getting severely confused.  This is because the stack variable
 * <tt>comm</tt>, which may be just an integer, disappears at the end
 * of the function.  Its address would no longer point to valid memory
 * after the function returns.
 *
 * The following code is syntactically correct, but may leak memory:

 \code

 // THIS CODE LEAKS MEMORY FOR GENERAL MPI_Comm OBJECTS.
 Teuchos::RCP<MPI_Comm> getMpiCommPtr (MPI_Comm comm)
 {
   MPI_Comm *pComm = new MPI_Comm (comm);
   // Works for comm==MPI_COMM_WORLD or MPI_COMM_SELF;
   // leaks memory for user-created MPI_Comm objects.
   return Teuchos::rcp (pComm);
 }

 \endcode

 * The above implementation of getMpiCommPtr() is correct only for the
 * standard MPI_Comm objects provided by MPI, like MPI_COMM_WORLD and
 * MPI_COMM_SELF.  It is <i>not</i> correct, and in fact may leak
 * memory, for custom MPI_Comm objects that the user creates by
 * calling functions like MPI_Comm_split().  This is because
 * user-created MPI_Comm objects must be freed by MPI_Comm_free().
 * Other kinds of opaque objects, like MPI_Datatype and MPI_Op, have
 * their own free functions.  Thus, even if opaque handles have the
 * type integer, they really behave like pointers or references.  Some
 * of them can and should be freed at the end of their useful lives;
 * others must not.  (Compare std::ostream; std::cout should never be
 * closed by typical user code, but an output file should be closed.)
 *
 * \section Teuchos_OpaqueWrapper_How How to use OpaqueWrapper
 *
 * We fix this problem by providing the OpaqueWrapper template base
 * class and the opaqueWrapper() nonmember template function.  Use
 * this function to wrap an opaque handle (like an MPI_Comm) in an
 * RCP.  This ensures that the RCP does the right thing in case the
 * handle must be freed.  For example, to wrap MPI_COMM_WORLD in a
 * RCP, just do this:

 \code
 RCP<OpaqueWrapper<MPI_Comm> > comm = opaqueWrapper (MPI_COMM_WORLD);
 \endcode

 * If you instead want to create a custom MPI_Comm using a function
 * like MPI_Comm_split(), then you may wrap it in an RCP as follows:

 \code
 MPI_Comm rawComm;
 // We omit all arguments but the last of MPI_Comm_split, for clarity.
 int errCode = MPI_Comm_split (..., &rawComm);
 if (errCode != MPI_SUCCESS) {
   // ... Handle the error ...
 }
 // See note below on using MPI_Comm_free here.
 RCP<OpaqueWrapper<MPI_Comm> > comm = opaqueWrapper (rawComm, MPI_Comm_free);
 \endcode

 * The optional second argument to opaqueWrapper() is a "free"
 * function, of type OpaqueFree which is a template parameter.  If the
 * free function is provided, then when the RCP's reference count goes
 * to zero, that function is called to "free" the handle.  If
 * opaqueFree is a free function, then the following must be
 * syntactically valid:
 \code
 Opaque opaque;
 opaqueFree (&opaque);
 \endcode
 * The function's return value, if any, is ignored.  Users are
 * responsible for knowing whether to provide a free function to
 * opaqueWrapper().  In this case, because we created an MPI_Comm
 * dynamically using a communicator "constructor" function, the
 * MPI_Comm must be "freed" after use.  RCP will automatically call
 * the "free" function once the reference count of \c comm reaches
 * zero.
 *
 * Users are not allowed to construct an OpaqueWrapper object
 * explicitly.  You must use the opaqueWrapper() nonmember function to
 * do so.
 *
 * \section Teuchos_OpaqueWrapper_Note Note on <tt>MPI_*_free</tt> functions, memory leaks, and MPI_Finalize
 *
 * This is a technical note that does not apply to most users.  It
 * will only apply to you if you are very concerned about small memory
 * leaks at the end of your program.  The explanations require
 * familiarity with RCP semantics and some familiarity with MPI's
 * attributes (a.k.a. key/value pair) facility.
 *
 * Functions of the form <tt>MPI_*_free</tt> (such as MPI_Comm_free)
 * should not be called after MPI_Finalize has been called.  
 * This issue has come up before; see e.g., 
 * <a href="https://software.sandia.gov/bugzilla/show_bug.cgi?id=5724">Bug 5724</a>.
 * There are two ways to deal with this:
 * 1. Make sure that all references to wrapped Opaque instances go
 *    away before MPI_Finalize is called.
 * 2. Write (or use) a "free" function that checks first whether
 *    MPI_Finalize has been called (by calling MPI_Finalized), before
 *    invoking any <tt>MPI_*_free</tt> functions.
 *
 * The details::safeCommFree function is a "free" function for
 * MPI_Comm that does the latter; it is declared in
 * Teuchos_DefaultMpiComm.hpp.  The "check MPI_Finalized first"
 * approach may permit a small memory leak, if your program allows
 * references to wrapped Opaque instances to persist past
 * MPI_Finalize, and if those Opaque instances need to be freed.  This
 * should not matter if the MPI implementation correctly frees all
 * resources before your program exits.
 *
 * If small memory leaks are not acceptable, even at program completion, 
 * then you should start with the approach discussed in Section 8.7.1 of the 
 * <a href="http://www.mpi-forum.org/docs/mpi-3.0/mpi30-report.pdf">MPI 3.0 Standard</a>,
 * namely, attaching attributes to MPI_COMM_SELF, with functions
 * that are called automatically by MPI_Finalize before it tears down
 * MPI_COMM_SELF.  Here are two ways to use this technique.
 *
 * First, you could forgo providing a "free" function to
 * opaqueWrapper(), and simply attach the raw opaque handle as an
 * attribute to MPI_COMM_SELF.  This makes the MPI_Comm persist
 * until MPI_Finalize is called.
 *
 * Second, if you don't want the MPI_Comm to persist that long, you
 * could develop a custom reference counting mechanism that
 * interoperates with MPI's attributes.  You could do this with just
 * a combination of RCP and the MPI attributes facility.  For example:
 * 1. Give a "free" function to opaqueWrapper() that first checks
 *    whether the Opaque handle is invalid (e.g.,
 *    <tt>MPI_*_NULL</tt>), frees it if it is not invalid, and then
 *    invalidates the Opaque handle after freeing it (e.g., by setting
 *    it to <tt>MPI_*_NULL</tt>).  Many (all?) <tt>MPI_*_free</tt>
 *    functions invalidate their input handle anyway.
 * 2. Make the <tt>RCP<OpaqueWrapper<Opaque> ></tt> itself an
 *    attribute of MPI_COMM_SELF.  (You'll have to pass in as a "new"
 *    <tt>RCP*</tt>.)
 * 3. Set the attribute's "free" function to call the appropriate
 *    <tt>MPI_*_free</tt> function if the underlying Opaque handle has
 *    not been invalidated.  This bypasses the OpaqueWrapper's normal
 *    "free" function mechanism and invalidates all references.  This
 *    is what you want, since no user code should use MPI objects
 *    after MPI_Finalize has been called.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template <class Opaque>
class OpaqueWrapper {
public:
  /// \brief Constructor that accepts and wraps a raw handle.
  ///
  /// Users typically never have to invoke the constructor explicitly.
  /// The opaqueWrapper() nonmember template function does this for them.
  OpaqueWrapper( Opaque opaque )
    : opaque_(opaque)
    {}
  /// \brief Implicit type conversion from wrapper to raw handle.
  ///
  /// Users typically never have to convert directly from an
  /// OpaqueWrapper to the raw handle that it wraps.  For example, if
  /// you have an <tt>RCP<OpaqueHandle<T> ></tt>, just deferencing the
  /// RCP will return the raw handle via this implicit type conversion
  /// operator:
  /// \code
  /// // We omit the right-hand side of this assignment, for simplicity.
  /// RCP<OpaqueWrapper<T> > wrapped = ...;
  /// // RCP's operator* returns OpaqueWrapper<T>&.
  /// // In turn, the operator below automatically converts to T.
  /// T raw = *wrapped;
  /// \endcode
  operator Opaque () const
    { return opaque_; }
  /// \brief Explicit type conversion from wrapper to raw handle.
  ///
  /// Users typically never have to convert directly from an
  /// OpaqueWrapper to the raw handle that it wraps.  However,
  /// in case they do, we provide this operator.
  Opaque operator()() const
    { return opaque_; }
protected:
  /// \brief The actual handle.
  ///
  /// This is protected and not private so that OpaqueWrapperWithFree
  /// can access it.  In general, one should avoid using protected
  /// data, but it would be silly to add member functions just for
  /// this simple use case.
  Opaque  opaque_;
private:
  OpaqueWrapper(); // Not defined
  OpaqueWrapper(const OpaqueWrapper&); // Not defined
  OpaqueWrapper& operator=(const OpaqueWrapper&); // Not defined
};

/** \class OpaqueWrapperWithFree
 * \brief Subclass for wrapped opaque objects with a free function.
 * \tparam Opaque Type of the opaque object (a.k.a. handle).
 * \tparam OpaqueFree Type of the function for freeing the handle.
 *
 * \note Teuchos users do not normally need to use or refer to this
 *   class.  The opaqueWrapper() nonmember template function suffices
 *   for nearly all use cases.  Please see the documentation of the
 *   base class OpaqueWrapper.
 *
 * This subclass allows a client to easily wrap any opaque object that
 * needs a function to free it.  This function (or function object)
 * must be callable as:

 \code
  opaqueFree(&opaque);
 \endcode

 * Again, this is typical for the opaque objects implemented in MPI
 * for instance.  For example, in order to delete an MPI_Comm object
 * created by the user (not MPI_COMM_WORLD), you must use the function
 * MPI_Comm_free().  See the documentation of OpaqueWrapper for
 * examples of how to supply a function for freeing an opaque handle.
 *
 * \relates OpaqueWrapper
 */
template <class Opaque, class OpaqueFree>
class OpaqueWrapperWithFree : public OpaqueWrapper<Opaque> {
public:
  //! Constructor: takes the opaque handle, and its free function.
  OpaqueWrapperWithFree( Opaque opaque, OpaqueFree opaqueFree )
    : OpaqueWrapper<Opaque>(opaque), opaqueFree_(opaqueFree)
  {}
  //! Destructor: invokes the free function.
  ~OpaqueWrapperWithFree()
  {
    // FIXME (mfh 10 Sep 2012) This only works if the free function is
    // a raw function pointer, not if it is a general "function object"
    // (i.e., something callable via operator()).
    if(opaqueFree_) {
#ifdef TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE
      Teuchos::RCP<Teuchos::FancyOStream>
	out = Teuchos::VerboseObjectBase::getDefaultOStream();
      Teuchos::OSTab tab(out);
      *out << "\nOpaqueWrapperWithFree::~OpaqueWrapperWithFree(): Freeing opaque object"
	   << " of type " << TypeNameTraits<Opaque>::name() << "!\n";
#endif // TEUCHOS_OPAQUE_WRAPPER_ANNOUNCE_FREE
      opaqueFree_(&this->opaque_);
    }
  }
private:
  //! Function (or function object) for freeing the handle.
  OpaqueFree opaqueFree_;
  OpaqueWrapperWithFree(); // Not defined
  OpaqueWrapperWithFree(const OpaqueWrapperWithFree&); // Not defined
  OpaqueWrapperWithFree& operator=(const OpaqueWrapperWithFree&); // Not defined
};


/** \brief Create a new <tt>OpaqueWrapper</tt> object without a free function.
 *
 * See the documentation of OpaqueWrapper for a detailed explanation
 * of why and how to use this function.
 *
 * \relates OpaqueWrapper
 */
template <class Opaque>
inline
RCP<OpaqueWrapper<Opaque> >
opaqueWrapper( Opaque opaque)
{
  return rcp(new OpaqueWrapper<Opaque>(opaque));
}


/** \brief Create a new <tt>OpaqueWrapper</tt> object with a free function.
 *
 * See the documentation of OpaqueWrapper for a detailed explanation
 * of why and how to use this function.
 *
 * \relates OpaqueWrapper
 */
template <class Opaque, class OpaqueFree>
inline
RCP<OpaqueWrapper<Opaque> >
opaqueWrapper( Opaque opaque, OpaqueFree opaqueFree)
{
  return rcp(new OpaqueWrapperWithFree<Opaque,OpaqueFree>(opaque,opaqueFree));
}


} // namespace Teuchos


#endif  // TEUCHOS_OPAQUE_WRAPPER_HPP
