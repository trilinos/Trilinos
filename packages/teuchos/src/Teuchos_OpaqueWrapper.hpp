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
 * \subsection Summary
 *
 * If you want to <i>create</i> an RCP to an opaque object, use the
 * opaqueWrapper() nonmember template function.  The <i>type</i> of an
 * RCP to an opaque object T is <tt>RCP<OpaqueWrapper<T> ></tt>.
 *
 * \subsection What are opaque objects (a.k.a. handles)?
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
 * Opaque handles are a useful technique, but they interfere with
 * correct use of RCP, as we will explain below.  This base class
 * allows opaque objects to be wrapped by a real object, whose address
 * you can take.  This is needed in order to wrap an opaque object in
 * a RCP, for example.
 *
 * \subsection Why do opaque handles need special treatment?
 *
 * This class was motivated in particular by MPI's common use of the
 * opaque object idiom.  For MPI, passing MPI_Comm, MPI_Datatype, and
 * MPI_Op objects around by handles hides implementation details from
 * the user.  Handles also make it easier to access MPI functionality
 * from Fortran, so that C, C++, and Fortran can all share the same
 * handle mechanism.  In fact, some MPI implementations (such as
 * MPICH, at least historically if not currently) simply implement
 * these handles all as integers.  (Presumably, such an implementation
 * would need to maintain a table on each process that maps the
 * integer value to a pointer to the corresponding object.)  For
 * example, MPI_Comm might be a typedef to int, and MPI_COMM_WORLD
 * might be a literal integer value:
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
 * We work around this by providing the OpaqueWrapper template base
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
 int errCode = MPI_Comm_split (..., &rawComm);
 if (errCode != MPI_SUCCESS) {
   // Handle the error
 }
 RCP<OpaqueWrapper<MPI_Comm> > comm = opaqueWrapper (rawComm, MPI_Comm_free);
 \endcode

 * The second argument to opaqueWrapper() is a "free" function.  It
 * takes a pointer to the opaque handle, and its return value (if any)
 * is ignored.  Users are responsible for knowing whether to provide a
 * free function to opaqueWrapper().
 *
 * Users are not allowed to construct an OpaqueWrapper object
 * explicitly.  You must use the opaqueWrapper() nonmember function to
 * do so.
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
  /// RCP<OpaqueWrapper<T> > wrapped = ...;
  /// T raw = *wrapped;
  /// \endcode
  operator Opaque () const
    { return opaque_; }
  /// \brief Explicit type conversion from wrapper to raw handle.
  ///
  /// Users typically never have to convert directly from an
  /// OpaqueWrapper to the raw handle that it wraps.
  Opaque operator()() const
    { return opaque_; }
protected:
  /// The actual handle.
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
  OpaqueWrapperWithFree( Opaque opaque, OpaqueFree opaqueFree )
    : OpaqueWrapper<Opaque>(opaque), opaqueFree_(opaqueFree)
    {}
  ~OpaqueWrapperWithFree()
    {
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
  OpaqueFree opaqueFree_;
  OpaqueWrapperWithFree(); // Not defined
  OpaqueWrapperWithFree(const OpaqueWrapperWithFree&); // Not defined
  OpaqueWrapperWithFree& operator=(const OpaqueWrapperWithFree&); // Not defined
};


/** \brief Helper function created a new <tt>OpaqueWrapper</tt> object without
 * a free function.
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


/** \brief Helper function created a new <tt>OpaqueWrapper</tt> object with a
 * free function.
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
