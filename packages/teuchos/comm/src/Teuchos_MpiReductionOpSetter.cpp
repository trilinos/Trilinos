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
// The two static variables theMpiOp_ and theMpiOpKey_ are persistent
// and initialized lazily.
//

// The MPI_Op singleton that implements the reduction or scan
// operation.  We only need to create the MPI_Op once (lazily, on
// demand).  When we create the MPI_Op, we stash its "destructor" in
// MPI_COMM_SELF so that it gets freed at MPI_Finalize.  (This is a
// standard MPI idiom.)
MPI_Op theMpiOp_ = MPI_OP_NULL;

// When we stash the MPI_Op singleton's "destructor" in MPI_COMM_SELF,
// we get back a key.  This is the key.  If we ever needed to change
// the destructor, we could use the key to look it up and free it.  We
// don't in this case, but knowing that the key value is valid (!=
// MPI_KEYVAL_INVALID) tells us whether the MPI_Op singleton has been
// created yet.
int theMpiOpKey_ = MPI_KEYVAL_INVALID; // flag value

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
// returned by freeMpiOp_.  This is the singleton's "destructor" that
// we attach to MPI_COMM_SELF as an MPI_Finalize hook.
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
// automatically.  Store the hook's key in theMpiOpKey_.
void createReductOp ()
{
  // This function has side effects on global singletons.  This check
  // ensures that the function is idempotent.  We only need to create
  // the MPI_Op singleton once.
  if (theMpiOpKey_ != MPI_KEYVAL_INVALID) {
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
  // theMpiOp_ gets freed at MPI_Finalize, if necessary.  Save the
  // resulting key in theMpiOpKey_.

  // 'key' is an output argument; we don't need to come up with a
  // unique key (thank goodness!).
  int key = 0;
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

  // OpenMPI 1.8.1 man page: "MPI_Comm_set_attr stores the stipulated
  // attribute value attribute_val for subsequent retrieval by
  // MPI_Comm_get_attr. If the value is already present, then the
  // outcome is as if MPI_Comm_delete_attr was first called to delete
  // the previous value (and the callback function delete_fn was
  // executed), and a new value was next stored."
  err = MPI_Comm_set_attr (MPI_COMM_SELF, key, &val);
  if (err != MPI_SUCCESS) {
    // Attempt to clean up by freeing the newly created MPI_Op.  If
    // cleaning up fails, just let it slide, since we're already in
    // trouble if MPI can't create a (key,value) pair.
    (void) MPI_Op_free (&mpi_op);
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "Teuchos::createReductOp: "
      "MPI_Comm_set_attr (for custom reduction operator) failed!");
  }

  // The "transaction" succeeded; save the results.
  theMpiOp_ = mpi_op;
  theMpiOpKey_ = key;
}

void
setReductOp (const Teuchos::Details::MpiReductionOpBase* reductOp)
{
  if (theMpiOpKey_ == MPI_KEYVAL_INVALID) {
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
  TEUCHOS_TEST_FOR_EXCEPTION(
    theMpiOpKey_ == MPI_KEYVAL_INVALID || theMpiOp_ == MPI_OP_NULL,
    std::logic_error, "Teuchos::Details::setMpiReductionOp: Failed to create "
    "theMpiOpKey_ or theMpiOp_.  This should never happen.  "
    "Please report this bug to the Teuchos developers.");
  return theMpiOp_;
}

} // namespace Details
} // namespace Teuchos

#endif // HAVE_MPI
