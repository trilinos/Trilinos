// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
 * \ingroup teuchos_mem_mng_grp
 *
 * \section Teuchos_OpaqueWrapper_Summary Summary
 *
 * If you want to create an RCP to an opaque handle (an instance of a
 * type like MPI_Comm), use the opaqueWrapper() nonmember template
 * function.  If that opaque handle needs to be freed after all
 * references to it go away, then supply a "free" function to
 * opaqueWrapper().  The type returned by opaqueWrapper() is
 * <tt>RCP<OpaqueWrapper<T> ></tt>, an RCP (reference-counted "smart"
 * pointer) to a wrapper of the opaque handle type T.  Users are not
 * allowed to construct an OpaqueWrapper object explicitly.  You must
 * use the opaqueWrapper() nonmember function to do so.
 *
 * \section Teuchos_OpaqueWrapper_Prereq Prerequisites
 *
 * In order to understand this documentation, you must first have
 * learned how to use RCP (Teuchos' reference-counted "smart" pointer
 * class) to manage dynamically allocated memory and other resources.
 * It also helps to be familiar with MPI (the Message Passing
 * Interface for distributed-memory parallel programming), but this is
 * not required.
 *
 * \section Teuchos_OpaqueWrapper_Handles What are opaque handles?
 *
 * Many different software libraries use the <i>opaque handle</i>
 * (a.k.a. opaque object) idiom to hide the internals of a data
 * structure from users.  This standard technique allows users to
 * treat an instance of a data structure as a handle.  Users may pass
 * the handle around as if it were a simple value type (like int), and
 * must call nonmember functions in order to create, operate on, use,
 * or destroy instances of the data structure.  The MPI (Message
 * Passing Interface) standard is full of examples of the opaque
 * handle idiom, including MPI_Comm (for communicators), MPI_Datatype
 * (for standard and custom data types), and MPI_Op (for standard and
 * custom reduction operations).
 *
 * In general, opaque handles (corresponding to the Opaque template
 * parameter) must be <i>assignable</i>.  This means that copy
 * construction and assignment (operator=) must be syntactically
 * correct for instances of Opaque.  This is certainly true of MPI's
 * opaque handle types.
 *
 * Opaque handles are a useful technique, but they interfere with
 * correct use of reference-counted "smart" pointer types such as
 * Teuchos' RCP or std::shared_ptr.  We will explain below why this is
 * the case.  The OpaqueWrapper base class allows opaque handles to be
 * wrapped by a real object, whose address you can take.  This is
 * needed in order to wrap an opaque object in a RCP, for example.
 *
 * \section Teuchos_OpaqueWrapper_Special Why do opaque handles need special treatment?
 *
 * The OpaqueWrapper class was motivated by MPI's common use of the
 * opaque handle idiom.  For MPI, passing MPI_Comm, MPI_Datatype, and
 * MPI_Op objects around by handles hides implementation details from
 * the user.  Handles also make it easier to access MPI functionality
 * from Fortran, so that C, C++, and Fortran can all share the same
 * handle mechanism.  In fact, some MPI implementations (such as
 * MPICH, at least historically if not currently) simply implement
 * these handles all as integers.  (As the MPI standard's advice to
 * implementers suggests, such an implementation would likely maintain
 * a table for each MPI process that maps the integer value to a
 * pointer to the corresponding object.)  For example, MPI_Comm might
 * be a typedef to int, and MPI_COMM_WORLD might be a C preprocessor
 * macro for a literal integer value:
 *
 \code
 typedef int MPI_Comm;
 #define MPI_COMM_WORLD 42
 \endcode

 * In this case, the expression <tt>rcp(&MPI_COMM_WORLD)</tt> would
 * not even compile, since one cannot take the address of an integer
 * literal such as 42.  (Remember that preprocessor macros get
 * replaced with their values before the C++ compiler does its work.)
 * To make this expression compile, one might try the following:

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
 * like MPI_Comm_split(), then you may wrap it in an RCP as follows
 * (please see discussion later about MPI_Comm_free()):

 \code
 MPI_Comm rawComm;
 // We omit all arguments but the last of MPI_Comm_split, for clarity.
 int errCode = MPI_Comm_split (..., &rawComm);
 if (errCode != MPI_SUCCESS) {
   // ... Handle the error ...
 }
 RCP<OpaqueWrapper<MPI_Comm> > comm = opaqueWrapper (rawComm, MPI_Comm_free);
 \endcode

 * The optional second argument to opaqueWrapper() is a "free"
 * function.  It has type OpaqueFree which is a template parameter.
 * If the free function is provided, then when the RCP's reference
 * count goes to zero, that function is called to "free" the handle.
 * If opaqueFree is a free function, then the following must be
 * syntactically valid, where <tt>opaque</tt> has type <tt>Opaque</tt>:

 \code
 opaqueFree (&opaque);
 \endcode

 * The function's return value, if any, is ignored.  Furthermore, the
 * OpaqueFree type must be copy constructible.  (A function pointer is
 * trivally copy constructible.)
 *
 * Users are responsible for knowing whether to provide a free
 * function to opaqueWrapper().  In this case, because we created an
 * MPI_Comm dynamically using a communicator "constructor" function,
 * the MPI_Comm must be "freed" after use.  RCP will automatically
 * call the "free" function once the reference count of \c comm
 * reaches zero.
 *
 * \note The above example is only correct if the reference count of
 * \c comm will go to zero before MPI_Finalize is called.  This is
 * because it's not valid to call MPI_Comm_free after MPI_Finalize has
 * been called.  The details::safeCommFree function checks whether
 * MPI_Finalize has been called (via MPI_Finalized) before calling
 * MPI_Comm_free; you may use this function as the free function if
 * you are concerned about this.
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

 * It must also be copy constructible.  (A function pointer is
 * trivially copy constructible.)  Please refer to the documentation
 * of OpaqueWrapper for examples of how to supply a function for
 * freeing an opaque handle.
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
