/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_BlockMultiVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_TypeNameTraits.hpp"

namespace {
  using Kokkos::ALL;
  using Kokkos::subview;
  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::TypeNameTraits;
  using std::endl;

  //
  // UNIT TESTS
  //

  // Test BlockMultiVector::blockWiseMultiply (analog of
  // MultiVector::elementWiseMultiply).
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMultiVector, BlockWiseMultiply, Scalar, LO, GO, Node )
  {
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef typename BMV::device_type device_type;
    typedef typename BMV::impl_scalar_type IST;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    typedef Kokkos::Details::ArithTraits<IST> KAT;
    typedef typename KAT::mag_type MT;
    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    int lclSuccess = 1;
    int gblSuccess = 0;

    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test BlockMultiVector::blockWiseMultiply" << endl;
    Teuchos::OSTab tab1 (out);
    myOut << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
          << "LO: " << TypeNameTraits<LO>::name () << endl
          << "GO: " << TypeNameTraits<LO>::name () << endl
          << "device_type: " << TypeNameTraits<device_type>::name () << endl;

    const LO blockSize = 4;
    const LO numVecs = 3;
    const LO numLocalMeshPoints = 12;
    const GO indexBase = 0;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm ();

    // Creating the Map first initializes the Kokkos execution space,
    // if it hasn't already been initialized.
    map_type meshMap (INVALID, static_cast<size_t> (numLocalMeshPoints),
                      indexBase, comm);
    // Make sure that the execution space actually got initialized.
    TEST_ASSERT( Kokkos::is_initialized () );

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      // mfh 24 Jan 2016: Failure to initialize Kokkos::Threads in a
      // release build before allocating a Kokkos::View may cause the
      // test to hang.
      myOut << "The execution space is not initialized.  No sense in "
        "continuing the test." << endl;
      return;
    }

    myOut << "Make prototype block and test GETF2 and GETRI" << endl;

    // Each block in the block diagonal will be a diagonally dominant
    // but nonsymmetric matrix.  To make sure that the implementation
    // doesn't mix up the blocks, we scale each block by a different
    // factor.  Ditto for the right-hand side.
    const IST zero = KAT::zero ();
    const IST one = KAT::one ();
    const IST two = one + one;
    const IST four = two + two;
    const IST five = four + one;

    typename Kokkos::View<IST**, device_type>::HostMirror prototypeBlock_h ("prototypeBlock", blockSize, blockSize);
    Kokkos::deep_copy (prototypeBlock_h, zero);
    for (LO i = 0; i < blockSize; ++i) {
      prototypeBlock_h(i,i) = four; // diagonally dominant
      if (i >= 1) {
        prototypeBlock_h(i, i-1) = one;
      }
      if (i + 1 < blockSize) {
        prototypeBlock_h(i, i+1) = -one; // nonsymmetric
      }
    }

    // Make a Teuchos copy of the matrix, so that we can use
    // Teuchos::BLAS to factor it.  (Kokkos::View may have a different
    // layout than the BLAS knows how to accept.
    // Teuchos::SerialDenseMatrix is always column major.)
    Teuchos::SerialDenseMatrix<int, Scalar> teuchosBlock (blockSize, blockSize);
    for (LO j = 0; j < blockSize; ++j) {
      for (LO i = 0; i < blockSize; ++i) {
        teuchosBlock(i,j) = prototypeBlock_h(i,j);
      }
    }

    // Use LAPACK (through the BLAS interface) to compute the inverse
    // (in place) of teuchosBlock.  We will use this to check our
    // implementation of the LU factorization and explicit inverse.
    typename Kokkos::View<int*, device_type>::HostMirror ipiv ("ipiv", blockSize);
    Teuchos::LAPACK<int, Scalar> lapack;
    typename Kokkos::View<IST*, device_type>::HostMirror work ("work", blockSize);
    int info = 0;
    {
      lapack.GETRF (blockSize, blockSize, teuchosBlock.values (),
                    teuchosBlock.stride (), ipiv.data (),
                    &info);
      TEST_EQUALITY_CONST( info, 0 );
      if (info != 0) {
        myOut << "GETRF returned info = " << info << " != 0 on the test matrix.  "
          "No point in continuing." << endl;
        return;
      }
      IST workQuery = zero;
      lapack.GETRI (blockSize, teuchosBlock.values (),
                    teuchosBlock.stride (), ipiv.data (),
                    reinterpret_cast<Scalar*> (&workQuery), -1, &info);
      TEST_EQUALITY_CONST( info, 0 );
      if (info != 0) {
        myOut << "GETRI workspace query returned info = " << info << " != 0.  "
          "No point in continuing." << endl;
        return;
      }
      const int lwork = static_cast<int> (KAT::real (workQuery));
      if (work.extent (0) < static_cast<size_t> (lwork)) {
        work = decltype (work) ("work", lwork);
      }
      lapack.GETRI (blockSize, teuchosBlock.values (),
                    teuchosBlock.stride (), ipiv.data (),
                    reinterpret_cast<Scalar*> (work.data ()),
                    lwork, &info);
      TEST_EQUALITY_CONST( info, 0 );
      if (info != 0) {
        myOut << "GETRI returned info = " << info << " != 0.  "
          "No point in continuing." << endl;
        return;
      }
    }

    Tpetra::Experimental::GETF2 (prototypeBlock_h, ipiv, info);
    TEST_EQUALITY_CONST( info, 0 );
    if (info != 0) {
      myOut << "Our GETF2 returned info = " << info << " != 0.  "
        "No point in continuing." << endl;
      return;
    }
    Tpetra::Experimental::GETRI (prototypeBlock_h, ipiv, work, info);
    TEST_EQUALITY_CONST( info, 0 );
    if (info != 0) {
      myOut << "Our GETF2 returned info = " << info << " != 0.  "
        "No point in continuing." << endl;
      return;
    }

    // Check that the explicit matrix inverse computed using our GETF2
    // and GETRI matches that computed by LAPACK.  I'm not sure if
    // it's legit to ask for componentwise accuracy, so I'll check
    // normwise accuracy.
    MT frobNorm = 0.0;
    for (LO i = 0; i < blockSize; ++i) {
      for (LO j = 0; j < blockSize; ++j) {
        frobNorm += KAT::abs (prototypeBlock_h(i,j) - static_cast<IST> (teuchosBlock(i,j)));
      }
    }
    myOut << "KAT::eps() = " << KAT::eps ()
        << " and blockSize = " << blockSize << endl;
    const MT maxMatFrobNorm = KAT::eps () * static_cast<MT> (blockSize) *
      static_cast<MT> (blockSize);
    myOut << "maxMatFrobNorm = " << maxMatFrobNorm << endl;
    TEST_ASSERT( frobNorm <= maxMatFrobNorm );

    if (! success) {
      myOut << "Returning early due to FAILURE." << endl;
      return;
    }

    // Compute the expected solution (little) vector in prototypeY,
    // when applying the explicit inverse to a vector of all ones.
    typename Kokkos::View<IST*, device_type>::HostMirror
      prototypeX ("prototypeX", blockSize);
    Kokkos::deep_copy (prototypeX, one);
    typename Kokkos::View<IST*, device_type>::HostMirror
      prototypeY ("prototypeY", blockSize);
    Teuchos::BLAS<int, Scalar> blas;
    blas.GEMV (Teuchos::NO_TRANS, blockSize, blockSize,
               static_cast<Scalar> (1.0),
               teuchosBlock.values (), teuchosBlock.stride (),
               reinterpret_cast<Scalar*> (prototypeX.data ()), 1,
               static_cast<Scalar> (0.0),
               reinterpret_cast<Scalar*> (prototypeY.data ()), 1);

    myOut << "Constructing block diagonal (as 3-D Kokkos::View)" << endl;

    // Each (little) block of D gets scaled by a scaling factor, and
    // every (little) vector of X gets scaled by the inverse of the
    // same factor.  This means the solution should not change.  The
    // point of this is to catch indexing errors.

    Kokkos::View<IST***, device_type> D ("D", numLocalMeshPoints,
                                         blockSize, blockSize);
    auto D_host = Kokkos::create_mirror_view (D);
    IST curScalingFactor = one;
    for (LO whichBlk = 0; whichBlk < numLocalMeshPoints; ++whichBlk) {
      auto D_cur = subview (D_host, whichBlk, ALL (), ALL ());
      Tpetra::Experimental::COPY (prototypeBlock_h, D_cur); // copy into D_cur
      Tpetra::Experimental::SCAL (curScalingFactor, D_cur);
      curScalingFactor += one;
    }
    Kokkos::deep_copy (D, D_host);

    myOut << "Make and fill BlockMultiVector instances" << endl;

    BMV X (meshMap, blockSize, numVecs);
    map_type pointMap = X.getPointMap ();
    BMV Y (meshMap, pointMap, blockSize, numVecs);

    X.putScalar (one);
    Y.putScalar (zero);
    {
      X.sync_host ();
      X.modify_host();
      curScalingFactor = one;
      for (LO whichBlk = 0; whichBlk < numLocalMeshPoints; ++whichBlk) {
        for (LO whichVec = 0; whichVec < numVecs; ++whichVec) {
          auto X_cur = X.getLocalBlock (whichBlk, whichVec);
          // This doesn't actually assume UVM, since we're modifying
          // the current block of X on the host.  The sync below will
          // sync back to device memory.
          Tpetra::Experimental::SCAL (one / curScalingFactor, X_cur);
        }
        curScalingFactor += one;
      }
      X.template sync<device_type> ();
    }

    myOut << "Call Y.blockWiseMultiply(alpha, D, X)" << endl;

    Y.blockWiseMultiply (one, D, X);

    myOut << "Check results of Y.blockWiseMultiply(alpha, D, X)" << endl;

    const MT maxVecNorm = KAT::eps () * static_cast<MT> (blockSize);
    myOut << "maxVecNorm = " << maxVecNorm << endl;

    {
      Y.sync_host ();
      Y.modify_host ();
      curScalingFactor = one;
      for (LO whichBlk = 0; whichBlk < numLocalMeshPoints; ++whichBlk) {
        for (LO whichVec = 0; whichVec < numVecs; ++whichVec) {
          auto Y_cur = Y.getLocalBlock (whichBlk, whichVec);

          // Compare Y_cur normwise to prototypeY.  This doesn't
          // actually assume UVM, since we're modifying the host
          // version of Y.  We will sync back to device below.
          MT frobNorm2 = 0.0;
          for (LO i = 0; i < blockSize; ++i) {
            frobNorm2 += KAT::abs (prototypeY(i) - Y_cur(i));
          }
          myOut << "frobNorm = " << frobNorm2 << endl;
          TEST_ASSERT( frobNorm2 <= maxVecNorm );
          if (! success) {
            myOut << "Retuning early due to FAILURE." << endl;
            return;
          }
        }
        curScalingFactor += one;
      }
      Y.template sync<device_type> ();
    }

    myOut << "Call Y.blockWiseMultiply(alpha, D, X) again, where Y has nonzero "
      "initial contents (the method should ignore those contents)" << endl;
    Y.putScalar (-five);
    Y.blockWiseMultiply (one, D, X);

    myOut << "Check results of Y.blockWiseMultiply(alpha, D, X)" << endl;

    {
      Y.sync_host ();
      curScalingFactor = one;
      for (LO whichBlk = 0; whichBlk < numLocalMeshPoints; ++whichBlk) {
        for (LO whichVec = 0; whichVec < numVecs; ++whichVec) {
          auto Y_cur = Y.getLocalBlock (whichBlk, whichVec);

          // Compare Y_cur normwise to prototypeY.  This doesn't
          // actually assume UVM, since we're modifying the host
          // version of Y.  We will sync back to device below.
          MT frobNorm2 = 0.0;
          for (LO i = 0; i < blockSize; ++i) {
            frobNorm2 += KAT::abs (prototypeY(i) - Y_cur(i));
          }
          myOut << "frobNorm = " << frobNorm2 << endl;
          TEST_ASSERT( frobNorm2 <= maxVecNorm );
        }
        curScalingFactor += one;
      }
      Y.template sync<device_type> ();
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (gblSuccess == 1) {
      myOut << "Test succeeded" << endl;
    }
    else {
      myOut << "Test FAILED" << endl;
    }
  }

  //
  // Test BlockMultiVector::blockJacobiUpdate.
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( BlockMultiVector, BlockJacobiUpdate, Scalar, LO, GO, Node )
  {
    typedef Tpetra::BlockMultiVector<Scalar, LO, GO, Node> BMV;
    typedef typename BMV::device_type device_type;
    typedef typename BMV::impl_scalar_type IST;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    typedef Tpetra::global_size_t GST;
    typedef Kokkos::Details::ArithTraits<IST> KAT;
    typedef typename KAT::mag_type MT;
    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    int lclSuccess = 1;
    int gblSuccess = 0;

    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test BlockMultiVector::blockJacobiUpdate" << endl;
    Teuchos::OSTab tab1 (out);
    myOut << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
          << "LO: " << TypeNameTraits<LO>::name () << endl
          << "GO: " << TypeNameTraits<LO>::name () << endl
          << "device_type: " << TypeNameTraits<device_type>::name () << endl;

    const LO blockSize = 4;
    const LO numVecs = 3;
    const LO numLocalMeshPoints = 12;
    const GO indexBase = 0;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm ();

    // Creating the Map first initializes the Kokkos execution space,
    // if it hasn't already been initialized.
    map_type meshMap (INVALID, static_cast<size_t> (numLocalMeshPoints),
                      indexBase, comm);
    // Make sure that the execution space actually got initialized.
    TEST_ASSERT( Kokkos::is_initialized () );

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      // mfh 24 Jan 2016: Failure to initialize Kokkos::Threads in a
      // release build before allocating a Kokkos::View may cause the
      // test to hang.
      myOut << "The execution space is not initialized.  No sense in "
        "continuing the test." << endl;
      return;
    }

    myOut << "Make prototype block" << endl;

    // Each block in the block diagonal will be a diagonally dominant
    // but nonsymmetric matrix.  To make sure that the implementation
    // doesn't mix up the blocks, we scale each block by a different
    // factor.  Ditto for the right-hand side.
    const IST zero = KAT::zero ();
    const IST one = KAT::one ();
    const IST two = one + one;
    const IST three = two + one;
    const IST four = two + two;
    const IST five = four + one;

    typename Kokkos::View<IST**, device_type>::HostMirror prototypeBlock_h ("prototypeBlock", blockSize, blockSize);
    Kokkos::deep_copy (prototypeBlock_h, zero);
    for (LO i = 0; i < blockSize; ++i) {
      prototypeBlock_h(i,i) = four; // diagonally dominant
      if (i >= 1) {
        prototypeBlock_h(i, i-1) = one;
      }
      if (i + 1 < blockSize) {
        prototypeBlock_h(i, i+1) = -one; // nonsymmetric
      }
    }

    typename Kokkos::View<int*, device_type>::HostMirror ipiv ("ipiv", blockSize);
    int info = 0;
    Tpetra::Experimental::GETF2 (prototypeBlock_h, ipiv, info);
    TEST_EQUALITY_CONST( info, 0 );
    if (info != 0) {
      myOut << "Our GETF2 returned info = " << info << " != 0.  "
        "No point in continuing." << endl;
      return;
    }
    typename Kokkos::View<IST*, device_type>::HostMirror work ("work", blockSize);
    Tpetra::Experimental::GETRI (prototypeBlock_h, ipiv, work, info);
    TEST_EQUALITY_CONST( info, 0 );
    if (info != 0) {
      myOut << "Our GETF2 returned info = " << info << " != 0.  "
        "No point in continuing." << endl;
      return;
    }

    myOut << "Constructing block diagonal (as 3-D Kokkos::View)" << endl;

    // Each (little) block of D gets scaled by a scaling factor, and
    // every (little) vector of X gets scaled by the inverse of the
    // same factor.  This means the solution should not change.  The
    // point of this is to catch indexing errors.

    Kokkos::View<IST***, device_type> D ("D", numLocalMeshPoints,
                                         blockSize, blockSize);
    auto D_host = Kokkos::create_mirror_view (D);

    IST curScalingFactor = one;
    for (LO whichBlk = 0; whichBlk < numLocalMeshPoints; ++whichBlk) {
      auto D_cur = subview (D_host, whichBlk, ALL (), ALL ());
      Tpetra::Experimental::COPY (prototypeBlock_h, D_cur); // copy into D_cur
      Tpetra::Experimental::SCAL (curScalingFactor, D_cur);
      curScalingFactor += one;
    }
    Kokkos::deep_copy (D, D_host);

    myOut << "Make and fill BlockMultiVector instances" << endl;

    BMV X (meshMap, blockSize, numVecs);
    map_type pointMap = X.getPointMap ();
    BMV Y (meshMap, pointMap, blockSize, numVecs);

    X.putScalar (one);
    Y.putScalar (zero);
    {
      X.sync_host ();
      X.modify_host ();
      curScalingFactor = one;
      for (LO whichBlk = 0; whichBlk < numLocalMeshPoints; ++whichBlk) {
        for (LO whichVec = 0; whichVec < numVecs; ++whichVec) {
          auto X_cur = X.getLocalBlock (whichBlk, whichVec);
          // This doesn't actually assume UVM, since we're modifying
          // the current block of X on the host.  The sync below will
          // sync back to device memory.
          Tpetra::Experimental::SCAL (one / curScalingFactor, X_cur);
        }
        curScalingFactor += one;
      }
      X.template sync<device_type> ();
    }

    // Fill Y with some initial value, so that using beta != 0 gives a
    // nontrivial result.  Fill Y2 to the same value; we will use it
    // to compute the same result in a different way.

    Y.putScalar (-five);
    BMV Y2 (meshMap, pointMap, blockSize, numVecs);
    Y2.putScalar (-five);

    const IST alpha = four / five;
    IST beta = -one / two; // just something different

    // Fill Z with some other value.  Remember that blockJacobiUpdate
    // reserves the right to treat Z as scratch space on output!
    // Thus, we will have to refill Z when computing Y2.
    BMV Z (meshMap, pointMap, blockSize, numVecs);
    Z.putScalar (one / three);

    myOut << "Call Y.blockJacobiUpdate(alpha, D, X, Z, beta) with alpha = "
        << alpha << " and beta = " << beta << endl;

    TEST_NOTHROW( Y.blockJacobiUpdate (alpha, D, X, Z, beta) );
    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      return; // abort the test
    }

    myOut << "Test result of blockJacobiUpdate" << endl;

    myOut << "Refill Z in preparation for the test" << endl;
    Z.putScalar (one / three);

    // We know that blockWiseMultiply works, so we can use it to test
    // blockJacobiUpdate.  Y2.blockJacobiUpdate(alpha, D, X, Z, beta)
    // computes Y2 := beta * Y2 + alpha * D * (X - Z), which we can also
    // compute as:
    //
    // Z := X - Z, that is, Z.update (1.0, X, -1.0);
    // U.blockWiseMultiply (alpha, D, Z);
    // Y2 := beta * Y2 + U, that is, Y2.update (1.0, U, beta);
    //
    // We use Y2 instead of Y here, so that we can compare against the
    // result of Y computed above.

    myOut << "Use blockWiseMultiply to test blockJacobiUpdate" << endl;
    Z.update (one, X, -one);
    BMV U (meshMap, pointMap, blockSize, numVecs);
    TEST_NOTHROW( U.blockWiseMultiply (alpha, D, Z) );
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (gblSuccess != 1) {
        myOut << "U.blockWiseMultiply(alpha,D,Z) threw an exception on some "
          "process!" << endl;
        return; // abort the test
      }
    }
    TEST_NOTHROW( Y2.update (one, U, beta) );
    {
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY_CONST( gblSuccess, 1 );
      if (gblSuccess != 1) {
        myOut << "Y2.update(one,U,beta) threw an exception on some process!"
            << endl;
        return; // abort the test
      }
    }

    // Take the norm of the difference between the two results.  We
    // can use the inf norm as a reasonable approximation for the
    // actual norm that we should use (something like the inf norm
    // between blocks, but 2-norm within a block).

    myOut << "Take the norm of the difference between the two results" << endl;
    Y2.update (one, Y, -one); // Y2 := Y - Y2
    Kokkos::View<MT*, device_type> norms ("norms", numVecs);
    Y2.getMultiVectorView ().normInf (norms);

    auto norms_host = Kokkos::create_mirror_view (norms);
    Kokkos::deep_copy (norms_host, norms);

    const MT maxResultNorm = KAT::eps () * static_cast<MT> (numLocalMeshPoints) *
      static_cast<MT> (blockSize);

    for (LO j = 0; j < numVecs; ++j) {
      myOut << "Norm of vector " << j << ": " << norms_host(j) << endl;
      TEST_ASSERT( norms_host(j) <= maxResultNorm );
    }

    myOut << "Retest Y.blockJacobiUpdate with beta == 0" << endl;
    // Retest Y.blockJacobiUpdate(alpha, D, X, Z, beta) with beta = 0.
    // The implementation should treat this as a special case by
    // ignoring the initial contents of Y.
    beta = zero;
    Y.putScalar (-five); // don't use zero, since we're testing beta = 0
    Y2.putScalar (-five); // same as Y
    Z.putScalar (one / three); // refill Z, since we used it as scratch

    myOut << "Call Y.blockJacobiUpdate(alpha, D, X, Z, beta) with alpha = "
        << alpha << " and beta = " << beta << endl;

    Y.blockJacobiUpdate (alpha, D, X, Z, beta);

    myOut << "Test result of blockJacobiUpdate" << endl;

    // Refill Z in preparation for the test.
    Z.putScalar (one / three);

    // We know that blockWiseMultiply works, so we can use it to test
    // blockJacobiUpdate.  Y2.blockJacobiUpdate(alpha, D, X, Z, 0.0)
    // computes Y2 := alpha * D * (X - Z), which we can also compute
    // as:
    //
    // Z := X - Z, that is, Z.update (1.0, X, -1.0);
    // Y2.blockWiseMultiply (alpha, D, Z);
    //
    // We use Y2 instead of Y here, so that we can compare against the
    // result of Y computed above.

    Z.update (one, X, -one);
    Y2.blockWiseMultiply (alpha, D, Z);

    // Take the norm of the difference between the two results.

    Y2.update (one, Y, -one); // Y2 := Y - Y2
    Y2.getMultiVectorView ().normInf (norms);
    Kokkos::deep_copy (norms_host, norms);

    for (LO j = 0; j < numVecs; ++j) {
      myOut << "Norm of vector " << j << ": " << norms_host(j) << endl;
      TEST_ASSERT( norms_host(j) <= maxResultNorm );
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess == 1) {
      myOut << "Test succeeded" << endl;
    }
    else {
      myOut << "Test FAILED" << endl;
    }
  }


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMultiVector, BlockWiseMultiply, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( BlockMultiVector, BlockJacobiUpdate, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)


