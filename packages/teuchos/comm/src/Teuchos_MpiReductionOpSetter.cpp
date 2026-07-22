// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_MpiReductionOpSetter.hpp"

#ifdef HAVE_MPI
#  ifdef MPIAPI
#    define CALL_API MPIAPI
#  else
#    define CALL_API
#  endif

//
// mfh 23 Nov 2014: My commits over the past day or two attempt to
// address Bug 6263.  In particular, the code as I found it had the
// following issues:
//
//   1. Static RCP instances (that persist past return of main())
//   2. Static MPI_Op datum (that persists past MPI_Finalize())
//   3. Code won't work with MPI_THREAD_{SERIALIZED,MULTIPLE},
//      because it assumes that only one MPI_Op for reductions
//      is needed at any one time
//
// I'm neglecting Issue #3 for now and focusing on the first two
// issues.  #1 goes away if one doesn't use RCPs and handles
// deallocation manually (we could also use std::shared_ptr, but that
// would require C++11).  #2 goes away with the standard idiom of an
// MPI_Finalize() hook (attach a (key,value) pair to MPI_COMM_SELF).
//

extern "C" {

// The MPI_Op that implements the reduction or scan operation will
// call this function.  We only need to create the MPI_Op once
// (lazily, on demand).  This function in turn will invoke
// theReductOp_ (see below), which gets set to the current reduction
// operation.  Thus, we only never need to create one MPI_Op, but we
// swap out the function.  This is meant to save overhead in creating
// and freeing MPI_Op for each reduction or scan.
void CALL_API
Teuchos_MPI_reduction_op (void* invec, void* inoutvec,
                          int* len, MPI_Datatype* datatype);
} // extern "C"

namespace { // anonymous

//
// theMpiOp_: The MPI_Op singleton that implements the Teuchos
// reduction or scan operation.  We only need to create the MPI_Op
// once (lazily, on demand).  When we create the MPI_Op, we stash its
// "destructor" in MPI_COMM_SELF so that it gets freed at
// MPI_Finalize.  (This is a standard MPI idiom.)
//
// This variable is global, persistent (until MPI_Finalize is called),
// and initialized lazily.
MPI_Op theMpiOp_ = MPI_OP_NULL;

// The current reduction or scan "function."  (It's actually a class
// instance.)
//
// This static variable is _NOT_ persistent.  It does not need
// deallocation.
const Teuchos::Details::MpiReductionOpBase* theReductOp_ = NULL;

// Free the given MPI_Op, and return the error code returned by MPI_Op_free.
int
freeMpiOp (MPI_Op* op)
{
  // If this function is called as an MPI_Finalize hook, MPI should
  // still be initialized at this point, and it should be OK to call
  // MPI functions.  Thus, we don't need to check if MPI is
  // initialized.
  int err = MPI_SUCCESS;
  if (op != NULL) {
    err = MPI_Op_free (op);
    if (err == MPI_SUCCESS) {
      // No externally visible side effects unless the above function succeeded.
      *op = MPI_OP_NULL;
    }
  }
  return err;
}

// Free the MPI_Op singleton (theMpiOp_), and return the error code
// returned by freeMpiOp().  As a side effect, if freeing succeeds,
// set theMpiOp_ to MPI_OP_NULL.
//
// This is the singleton's "destructor" that we attach to
// MPI_COMM_SELF as an MPI_Finalize hook.
int
freeMpiOpCallback (MPI_Comm, int, void*, void*)
{
  // We don't need any of the arguments to this function, since we're
  // just freeing the singleton.
  if (theMpiOp_ == MPI_OP_NULL) {
    return MPI_SUCCESS;
  } else {
    return freeMpiOp (&theMpiOp_);
  }
}

// Create the MPI_Op singleton that invokes the
// Teuchos_MPI_reduction_op callback.  Assign the MPI_Op to theMpiOp_,
// and set it up with an MPI_Finalize hook so it gets freed
// automatically.
void createReductOp ()
{
#if MPI_VERSION >= 2

  // This function has side effects on the global singleton theMpiOp_.
  // This check ensures that the function is idempotent.  We only need
  // to create the MPI_Op singleton once.
  if (theMpiOp_ != MPI_OP_NULL) {
    return; // We've already called this function; we don't have to again.
  }

  MPI_Op mpi_op = MPI_OP_NULL;

  // FIXME (mfh 23 Nov 2014) I found the following comment here:
  // "Assume op is commutative".  That's what it means to pass 1 as
  // the second argument.  I don't know whether it's a good idea to
  // keep that assumption.
  int err = MPI_Op_create (&Teuchos_MPI_reduction_op, 1, &mpi_op);
  TEUCHOS_TEST_FOR_EXCEPTION(
    err != MPI_SUCCESS, std::runtime_error, "Teuchos::createReductOp: "
    "MPI_Op_create (for custom reduction operator) failed!");

  // Use the standard MPI idiom (attach a (key,value) pair to
  // MPI_COMM_SELF with a "destructor" function) in order that
  // theMpiOp_ gets freed at MPI_Finalize, if necessary.

  // 'key' is an output argument of MPI_Comm_create_keyval.
  int key = MPI_KEYVAL_INVALID;
  err = MPI_Comm_create_keyval (MPI_COMM_NULL_COPY_FN, freeMpiOpCallback,
                                &key, NULL);
  if (err != MPI_SUCCESS) {
    // Attempt to clean up by freeing the newly created MPI_Op.  If
    // cleaning up fails, just let it slide, since we're already in
    // trouble if MPI can't create a (key,value) pair.
    (void) MPI_Op_free (&mpi_op);
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "Teuchos::createReductOp: "
      "MPI_Comm_create_keyval (for custom reduction operator) failed!");
  }
  int val = key; // doesn't matter

  // Attach the attribute to MPI_COMM_SELF.
  err = MPI_Comm_set_attr (MPI_COMM_SELF, key, &val);
  if (err != MPI_SUCCESS) {
    // MPI (versions up to and including 3.0) doesn't promise correct
    // behavior after any function returns something other than
    // MPI_SUCCESS.  Thus, it's not required to try to free the new
    // key via MPI_Comm_free_keyval.  Furthermore, if something went
    // wrong with MPI_Comm_set_attr, it's likely that the attribute
    // mechanism is broken.  Thus, it would be unwise to call
    // MPI_Comm_free_keyval.
    //
    // I optimistically assume that the "rest" of MPI is still
    // working, and attempt to clean up by freeing the newly created
    // MPI_Op.  If cleaning up fails, just let it slide, since we're
    // already in trouble if MPI can't create a (key,value) pair.
    (void) MPI_Op_free (&mpi_op);
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "Teuchos::createReductOp: "
      "MPI_Comm_set_attr (for custom reduction operator) failed!");
  }

  // It looks weird to "free" the key right away.  However, this does
  // not actually cause the "destructor" to be called.  It only gets
  // called at MPI_FINALIZE.  See MPI 3.0 standard, Section 6.7.2,
  // MPI_COMM_FREE_KEYVAL:
  //
  // "Note that it is not erroneous to free an attribute key that is
  // in use, because the actual free does not transpire until after
  // all references (in other communicators on the process) to the key
  // have been freed.  These references need to be explicitly freed by
  // the program, either via calls to MPI_COMM_DELETE_ATTR that free
  // one attribute instance, or by calls to MPI_COMM_FREE that free
  // all attribute instances associated with the freed communicator."
  //
  // We rely here on the latter mechanism.  MPI_FINALIZE calls
  // MPI_COMM_FREE on MPI_COMM_SELF, so we do not need to call it
  // explicitly.
  //
  // It's not clear what to do if the MPI_* calls above succeeded, but
  // this call fails (i.e., returns != MPI_SUCCESS).  We could throw;
  // this would make sense to do, because MPI (versions up to and
  // including 3.0) doesn't promise correct behavior after any MPI
  // function returns something other than MPI_SUCCESS.  We could also
  // be optimistic and just ignore the return value, hoping that if
  // the above calls succeeded, then the communicator will get freed
  // at MPI_FINALIZE, even though the unfreed key may leak memory (see
  // Bug 6338).  I've chosen the latter.
  (void) MPI_Comm_free_keyval (&key);

  // The "transaction" succeeded; save the result.
  theMpiOp_ = mpi_op;

#else // MPI_VERSION < 2
#  error "Sorry, you need an MPI implementation that supports at least MPI 2.0 in order to build this code.  MPI 2.0 came out in 1997.  I wrote this comment in 2017.  If you really _really_ want MPI 1.x support, please file a GitHub issue for this feature request at github.com/trilinos/trilinos/issues with an expression of its priority and we will get to it as soon as we can."
#endif // MPI_VERSION >= 2
}

void
setReductOp (const Teuchos::Details::MpiReductionOpBase* reductOp)
{
  if (theMpiOp_ == MPI_OP_NULL) {
    createReductOp ();
  }
  theReductOp_ = reductOp;
}

} // namespace (anonymous)

extern "C" {

void CALL_API
Teuchos_MPI_reduction_op (void* invec,
                          void* inoutvec,
                          int* len,
                          MPI_Datatype* datatype)
{
  if (theReductOp_ != NULL) {
    theReductOp_->reduce (invec, inoutvec, len, datatype);
  }
}

} // extern "C"

namespace Teuchos {
namespace Details {

MPI_Op setMpiReductionOp (const MpiReductionOpBase& reductOp)
{
  setReductOp (&reductOp);
  TEUCHOS_TEST_FOR_EXCEPTION
    (theMpiOp_ == MPI_OP_NULL, std::logic_error, "Teuchos::Details::"
     "setMpiReductionOp: Failed to create reduction MPI_Op theMpiOp_.  "
     "This should never happen.  "
     "Please report this bug to the Teuchos developers.");
  return theMpiOp_;
}

} // namespace Details
} // namespace Teuchos

#endif // HAVE_MPI
