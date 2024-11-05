// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_DEFAULT_COMM_HPP
#define TEUCHOS_DEFAULT_COMM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif

namespace Teuchos {

#ifdef HAVE_MPI
namespace Details {

template<class OrdinalType>
int
mpiFreeDefaultComm (MPI_Comm, int, void*, void*);

template<class OrdinalType>
int
mpiFreeDefaultSerialComm (MPI_Comm, int, void*, void*);

} // namespace Details
#endif // HAVE_MPI

/// \class DefaultComm
/// \brief Return a default global communicator appropriate for the build.
///
/// Use this class to get a Comm instance representing the default
/// global communicator.  If Teuchos was built with MPI (i.e., if the
/// HAVE_MPI macro is defined), then the default communicator wraps
/// MPI_COMM_WORLD.  Otherwise, it is a "serial" communicator
/// (containing one process, whose rank is zero).
///
/// \tparam OrdinalType The ordinal type for the Comm communicator
///   wrapper template class.  Comm uses \c OrdinalType to represent
///   things like array lengths and indices.
///
/// \note (mfh 19 Jul 2011, 22 Dec 2011) \c OrdinalType is called
///   <tt>OrdinalType</tt> and not \c Ordinal, because of a bug in
///   Intel's C++ compiler (version 11.1).  This compiler confuses the
///   \c Ordinal template parameter of DefaultComm with the
///   Teuchos::Ordinal typedef.  The \c Ordinal template parameter
///   should actually shadow the typedef in Teuchos, and it does with
///   GCC 4.5.1, but does not with Intel's compiler.  This may be the
///   case with other compilers as well, but I haven't tested them
///   yet.  If this class' template parameter were named \c Ordinal,
///   then the following line of code would result in a compile error:
///   \code
///   RCP<const Comm<int> pComm = DefaultComm<int>::getDefaultSerialComm (null);
///   \endcode
template<typename OrdinalType>
class DefaultComm {
public:
  /// \brief Return the default global communicator.
  ///
  /// \warning When running with MPI, do not call this function until
  ///   after MPI_Init has been called.  You can use GlobalMPISesssion
  ///   to initialize MPI without explicitly depending on the MPI
  ///   interface or the mpi.h header file.  (If Trilinos was not
  ///   built with MPI, GlobalMPISession will do the right thing, so
  ///   you can use it unconditionally.)
  ///
  /// \warning Do not use or refer to the returned object after
  ///   MPI_Finalize has been called.  In a non-MPI build, do not use
  ///   or refer to the returned object after main() has returned.
  static Teuchos::RCP<const Comm<OrdinalType> > getComm ();

  /// \brief Return a serial Comm if the input Comm is null.
  ///
  /// If the input communicator \c comm is null, return the default
  /// serial communicator.  Otherwise, just return the input.
  ///
  /// \warning The same warnings as for getComm() apply here.
  static Teuchos::RCP<const Comm<OrdinalType> >
  getDefaultSerialComm (const Teuchos::RCP<const Comm<OrdinalType> >& comm);

private:
  /// \brief The default global communicator.
  ///
  /// If Teuchos was built with MPI, this is a wrapper for
  /// MPI_COMM_WORLD.  Otherwise, this is a "serial" communicator
  /// (containing one process, whose rank is zero).
  static const Comm<OrdinalType>* comm_;

  //! A "serial" communicator (containing one process, whose rank is zero).
  static const Comm<OrdinalType>* defaultSerialComm_;

#ifdef HAVE_MPI
  //! MPI_Finalize hook that calls freeDefaultComm().
  template<class OT>
  friend int
  Details::mpiFreeDefaultComm (MPI_Comm, int, void*, void*);

  //! MPI_Finalize hook that calls freeDefaultSerialComm().
  template<class OT>
  friend int
  Details::mpiFreeDefaultSerialComm (MPI_Comm, int, void*, void*);
#endif // HAVE_MPI

  //! Free the default Comm object.
  static void freeDefaultComm () {
    if (comm_ != NULL) {
      delete comm_;
      comm_ = NULL;
    }
  }

  //! Free the default serial Comm object.
  static void freeDefaultSerialComm () {
    if (defaultSerialComm_ != NULL) {
      delete defaultSerialComm_;
      defaultSerialComm_ = NULL;
    }
  }
};

#ifdef HAVE_MPI
namespace Details {

template<class OrdinalType>
int
mpiFreeDefaultComm (MPI_Comm, int, void*, void*)
{
  try {
    ::Teuchos::DefaultComm<OrdinalType>::freeDefaultComm ();
  } catch (...) {
    // Destructors must not throw exceptions, so we must accept the
    // possible memory leak and move on.
    std::cerr << "Teuchos::DefaultComm: Failed to free default Comm!  We can't "
      "throw an exception here because this is a singleton destructor that "
      "should only be called at MPI_Finalize or (if not building with MPI) at "
      "exit from main()." << std::endl;
    // FIXME (mfh 16 Nov 2014) There might be some way to create a
    // custom return code with MPI error reporting.  For now, we just
    // pick some error code not equal to MPI_SUCCESS.  It could
    // perhaps overlap with some existing error code.
    return (MPI_SUCCESS == 0) ? -1 : 0;
  }
  return MPI_SUCCESS;
}

template<class OrdinalType>
int
mpiFreeDefaultSerialComm (MPI_Comm, int, void*, void*)
{
  try {
    ::Teuchos::DefaultComm<OrdinalType>::freeDefaultSerialComm ();
  } catch (...) {
    // Destructors must not throw exceptions, so we must accept the
    // possible memory leak and move on.
    std::cerr << "Teuchos::DefaultComm: Failed to free default serial Comm!  "
      "We can't throw an exception here because this is a singleton destructor "
      "that should only be called at MPI_Finalize or (if not building with MPI)"
      " at exit from main()." << std::endl;
    // FIXME (mfh 16 Nov 2014) There might be some way to create a
    // custom return code with MPI error reporting.  For now, we just
    // pick some error code not equal to MPI_SUCCESS.  It could
    // perhaps overlap with some existing error code.
    return (MPI_SUCCESS == 0) ? -1 : 0;
  }
  return MPI_SUCCESS;
}

} // namespace Details
#endif // HAVE_MPI


template<typename OrdinalType>
Teuchos::RCP<const Teuchos::Comm<OrdinalType> >
DefaultComm<OrdinalType>::getComm ()
{
  if (comm_ == NULL) {
#ifdef HAVE_MPI
#  if MPI_VERSION >= 2

    comm_ = new MpiComm<OrdinalType> (MPI_COMM_WORLD);

    // We want comm_ to be deallocated when MPI_Finalize is called.
    // The standard idiom for this (look in the MPI standard) is to
    // register an attribute ((key,value) pair) with MPI_COMM_SELF,
    // with a custom "destructor" to be called at MPI_Finalize.

    // 'key' is an output argument of MPI_Comm_create_keyval.
    int key = MPI_KEYVAL_INVALID;
    int err =
      MPI_Comm_create_keyval (MPI_COMM_NULL_COPY_FN,
                              Details::mpiFreeDefaultComm<OrdinalType>,
                              &key,
                              NULL); // no extra state
    if (err != MPI_SUCCESS) {
      if (comm_ != NULL) { // clean up if MPI call fails
        delete comm_;
        comm_ = NULL;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Teuchos::DefaultComm::getComm: MPI_Comm_create_keyval failed!");
    }
    int val = key; // doesn't matter

    // Attach the attribute to MPI_COMM_SELF.
    err = MPI_Comm_set_attr (MPI_COMM_SELF, key, &val);
    if (err != MPI_SUCCESS) {
      // MPI (versions up to and including 3.0) doesn't promise
      // correct behavior after any function returns something other
      // than MPI_SUCCESS.  Thus, it's not required to try to free the
      // new key via MPI_Comm_free_keyval.  Furthermore, if something
      // went wrong with MPI_Comm_set_attr, it's likely that the
      // attribute mechanism is broken.  Thus, it would be unwise to
      // call MPI_Comm_free_keyval.  However, we can still clean up
      // other data.
      if (comm_ != NULL) { // clean up if MPI call fails
        delete comm_;
        comm_ = NULL;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Teuchos::DefaultComm::getComm: MPI_Comm_set_attr failed!");
    }

    // It looks weird to "free" the key right away.  However, this
    // does not actually cause the "destructor" to be called.  It only
    // gets called at MPI_FINALIZE.  See MPI 3.0 standard, Section
    // 6.7.2, MPI_COMM_FREE_KEYVAL:
    //
    // "Note that it is not erroneous to free an attribute key that is
    // in use, because the actual free does not transpire until after
    // all references (in other communicators on the process) to the
    // key have been freed.  These references need to be explicitly
    // freed by the program, either via calls to MPI_COMM_DELETE_ATTR
    // that free one attribute instance, or by calls to MPI_COMM_FREE
    // that free all attribute instances associated with the freed
    // communicator."
    //
    // We rely here on the latter mechanism.  MPI_FINALIZE calls
    // MPI_COMM_FREE on MPI_COMM_SELF, so we do not need to call it
    // explicitly.
    //
    // It's not clear what to do if the MPI_* calls above succeeded,
    // but this call fails (i.e., returns != MPI_SUCCESS).  We could
    // throw; this would make sense to do, because MPI (versions up to
    // and including 3.0) doesn't promise correct behavior after any
    // MPI function returns something other than MPI_SUCCESS.  We
    // could also be optimistic and just ignore the return value,
    // hoping that if the above calls succeeded, then the communicator
    // will get freed at MPI_FINALIZE, even though the unfreed key may
    // leak memory (see Bug 6338).  I've chosen the latter.
    (void) MPI_Comm_free_keyval (&key);

#  else // MPI_VERSION < 2
#    error "Sorry, you need an MPI implementation that supports at least MPI 2.0 in order to build this code.  MPI 2.0 came out in 1997.  I wrote this comment in 2017.  If you really _really_ want MPI 1.x support, please file a GitHub issue for this feature request at github.com/trilinos/trilinos/issues with an expression of its priority and we will get to it as soon as we can."
#  endif // MPI_VERSION >= 2

#else // NOT HAVE_MPI
    comm_ = new SerialComm<OrdinalType> ();
    // We want comm_ to be deallocated when main exits, so register
    // its deallocation function as an atexit handler.
    //
    // The POSIX standard allows atexit to fail, in particular if it
    // lacks space for registering more functions.  "[T]he application
    // should call sysconf() to obtain the value of {ATEXIT_MAX}, the
    // [maximum] number of functions that can be registered. There is
    // no way for an application to tell how many functions have
    // already been registered with atexit()."
    //
    // We don't do this here.  Instead, we just check atexit's return
    // code.  If it fails, we throw.
    int err = atexit (freeDefaultComm);
    if (err != 0) {
      if (comm_ != NULL) { // clean up if atexit fails
        delete comm_;
        comm_ = NULL;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Teuchos::DefaultComm::getComm: atexit failed!");
    }
#endif // HAVE_MPI
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (comm_ == NULL, std::logic_error, "Teuchos::DefaultComm::getComm: "
     "comm_ == NULL before return.  This should never happen.  "
     "Please report this bug to the Teuchos developers.");

  // Return a nonowning RCP, because we need to ensure that
  // destruction happens at MPI_Finalize (or at exit of main(), if not
  // building with MPI).
  return rcp (comm_, false);
}

template<typename OrdinalType>
Teuchos::RCP<const Teuchos::Comm<OrdinalType> >
DefaultComm<OrdinalType>::
getDefaultSerialComm (const Teuchos::RCP<const Comm<OrdinalType> >& comm)
{
  if (! comm.is_null ()) {
    return comm;
  } else {
    if (defaultSerialComm_ == NULL) {
#ifdef HAVE_MPI
#  if MPI_VERSION >= 2
      //defaultSerialComm_ = new MpiComm<OrdinalType> (MPI_COMM_SELF);
      defaultSerialComm_ = new SerialComm<OrdinalType> ();

      // Register an MPI_Finalize hook to free defaultSerialComm_.
      // (See getComm implementation above in this file for details.)

      int key = MPI_KEYVAL_INVALID;
      int err =
        MPI_Comm_create_keyval (MPI_COMM_NULL_COPY_FN,
                                Details::mpiFreeDefaultSerialComm<OrdinalType>,
                                &key,
                                NULL); // no extra state
      if (err != MPI_SUCCESS) {
        if (defaultSerialComm_ != NULL) { // clean up if MPI call fails
          delete defaultSerialComm_;
          defaultSerialComm_ = NULL;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "Teuchos::DefaultComm::getDefaultSerialComm"
          ": MPI_Comm_create_keyval failed!");
      }
      int val = key; // doesn't matter

      // Attach the attribute to MPI_COMM_SELF.
      err = MPI_Comm_set_attr (MPI_COMM_SELF, key, &val);
      if (err != MPI_SUCCESS) {
        // See comments in getComm implementation above to see why we
        // don't call MPI_Comm_free_keyval here.
        if (defaultSerialComm_ != NULL) { // clean up if MPI call fails
          delete defaultSerialComm_;
          defaultSerialComm_ = NULL;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::runtime_error, "Teuchos::DefaultComm::getDefaultSerialComm"
          ": MPI_Comm_set_attr failed!");
      }

      // See comments in getComm implementation above to see why we
      // _do_ call MPI_Comm_free_keyval here, and why we don't check
      // the return code.
      (void) MPI_Comm_free_keyval (&key);

#  else // MPI_VERSION < 2
#    error "Sorry, you need an MPI implementation that supports at least MPI 2.0 in order to build this code.  MPI 2.0 came out in 1997.  I wrote this comment in 2017.  If you really _really_ want MPI 1.x support, please file a GitHub issue for this feature request at github.com/trilinos/trilinos/issues with an expression of its priority and we will get to it as soon as we can."
#  endif // MPI_VERSION >= 2

#else // NOT HAVE_MPI
      defaultSerialComm_ = new SerialComm<OrdinalType> ();
      // We want defaultSerialComm_ to be deallocated when main exits,
      // so register its deallocation function as an atexit handler.
      //
      // The POSIX standard allows atexit to fail, in particular if it
      // lacks space for registering more functions.  "[T]he
      // application should call sysconf() to obtain the value of
      // {ATEXIT_MAX}, the [maximum] number of functions that can be
      // registered. There is no way for an application to tell how
      // many functions have already been registered with atexit()."
      //
      // We don't do this here.  Instead, we just check atexit's
      // return code.  If it fails, we throw.
      int err = atexit (freeDefaultComm);
      if (err != 0) {
        if (defaultSerialComm_ != NULL) { // clean up if atexit fails
          delete defaultSerialComm_;
          defaultSerialComm_ = NULL;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "Teuchos::DefaultComm::getDefaultSerialComm: atexit failed!");
      }
#endif // HAVE_MPI
    }

    TEUCHOS_TEST_FOR_EXCEPTION
      (defaultSerialComm_ == NULL, std::logic_error, "Teuchos::DefaultComm::"
       "getDefaultSerialComm: defaultSerialComm_ == NULL before return.  This sh"
       "ould never happen.  Please report this bug to the Teuchos developers.");

    // Return a nonowning RCP, because we need to ensure that
    // destruction happens at MPI_Finalize (or at exit of main(), if not
    // building with MPI).
    return rcp (defaultSerialComm_, false);
  }
}

template<typename OrdinalType>
const Teuchos::Comm<OrdinalType>*
DefaultComm<OrdinalType>::comm_ = NULL;

template<typename OrdinalType>
const Teuchos::Comm<OrdinalType>*
DefaultComm<OrdinalType>::defaultSerialComm_ = NULL;

} // namespace Teuchos

#endif // TEUCHOS_DEFAULT_COMM_HPP
