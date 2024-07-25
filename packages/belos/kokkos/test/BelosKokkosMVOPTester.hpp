// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
#ifndef BELOS_KOKKOS_MVOPTESTER_HPP
#define BELOS_KOKKOS_MVOPTESTER_HPP

// Assumptions that I have made:
// * I assume/verify that a multivector must have at least one std::vector. This seems
//   to be consistent with Epetra_MultiVec.
// * I do not assume that an operator is deterministic; I do assume that the
//   operator, applied to 0, will return 0.

/*! \file BelosKokkosMVOPTester.hpp
  \brief Test routines for Kokkos MultiVecTraits and OperatorTraits conformity.
  Note: This tester file is a direct copy of the BelosMVOPTester.hpp, 
  except that the tests for accessing non-contiguous columns of data have
  been removed for the testing of the Kokkos adapter. 
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosOutputManager.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_SetScientific.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

namespace Belos {

  /// \brief Test correctness of a MultiVecTraits specialization and
  ///   multivector implementation.
  ///
  /// \tparam ScalarType The type of the entries in the multivectors;
  ///   the first template parameter of MultiVecTraits.
  /// \tparam MV The multivector type; the second template parameter
  ///   of MultiVecTraits.
  ///
  /// \param om [in/out] A valid OutputManager, for displaying test results.
  ///
  /// \param A [in] An initial multivector, to use for making new
  ///   multivectors.  (Belos doesn't currently have a "vector space"
  ///   abstraction; making a new multivector requires a valid input
  ///   multivector to clone.)
  ///
  /// \return Test result: true if all tests passed, else false.
  template< class ScalarType, class MV >
  bool
  TestKokkosMultiVecTraits (const Teuchos::RCP<OutputManager<ScalarType> > &om,
                      const Teuchos::RCP<const MV> &A)
  {
    using Teuchos::SetScientific;
    using std::endl;
    typedef MultiVecTraits<ScalarType, MV>    MVT;
    typedef Teuchos::ScalarTraits<ScalarType> STS;
    typedef typename STS::magnitudeType       MagType;

    // Make sure that all floating-point numbers are printed with the
    // right precision.
    SetScientific<ScalarType> sci (om->stream (Warnings));

    // FIXME (mfh 09 Jan 2013) Added an arbitrary tolerance in case
    // norms are not computed deterministically (which is possible
    // even with MPI only, and more likely with threads).
    const MagType tol = Teuchos::as<MagType> (100) * STS::eps ();

    /* MVT Contract:

         Clone(MV,int)
         CloneCopy(MV)
         CloneCopy(MV,vector<int>)
           USER: will request positive number of vectors
             MV: will return a multivector with exactly the number of
                   requested vectors.
                 vectors are the same dimension as the cloned MV


         CloneView(MV,vector<int>) [const and non-const]
           USER: There is no assumed communication between creation and
           destruction of a view. I.e., after a view is created, changes to the
           source multivector are not reflected in the view. Likewise, until
           destruction of the view, changes in the view are not reflected in the
           source multivector.

         GetGlobalLength 
             MV: will always be positive (MV cannot have zero vectors)

         GetNumberVecs
             MV: will always be positive (MV cannot have zero vectors)

         MvAddMv
           USER: multivecs will be of the same dimension and same number of vecs
             MV: input vectors will not be modified
                 performing C=0*A+1*B will assign B to C exactly

         MvTimesMatAddMv
           USER: multivecs and serialdensematrix will be of the proper shape
             MV: input arguments will not be modified
                 following arithmetic relations hold exactly:
                   A*I = A
                   0*B = B
                   1*B = B

         MvTransMv
           USER: SerialDenseMatrix will be large enough to hold results.
             MV: SerialDenseMatrix will not be resized.
                 Inner products will satisfy |a'*b| <= |a|*|b|
                 alpha == 0  =>  SerialDenseMatrix == 0

         MvDot
          USER: Results vector will be large enough for results.
                Both multivectors will have the same number of vectors.
                    (Epetra crashes, otherwise.)
            MV: Inner products will satisfy |a'*b| <= |a|*|b|
                Results vector will not be resized.

         MvNorm
             MV: vector norm is always non-negative, and zero
                   only for zero vectors.
                 results vector should not be resized

         SetBlock
          USER: indices will be distinct
            MV: assigns copies of the vectors to the specified
                locations, leaving the other vectors untouched.

         MvRandom
             MV: Generate zero vector with "zero" probability
                 Don't gen the same vectors twice.

         MvInit
             MV: Init(alpha) sets all elements to alpha

         MvScale (two versions)
             MV: scales multivector values

         MvPrint
             MV: routine does not modify vectors (not tested here)
    *********************************************************************/

    const ScalarType one      = STS::one();
    const ScalarType zero     = STS::zero();
    const MagType    zero_mag = Teuchos::ScalarTraits<MagType>::zero();

    // Don't change these two without checking the initialization of ind below
    const int numvecs   = 10;
    const int numvecs_2 = 5;

    std::vector<int> ind(numvecs_2);
    std::vector<int> indContig(numvecs_2);

    /* Initialize indices for selected copies/views
       The MVT specialization should not assume that
       these are ordered or even distinct.
       Also retrieve the edges.

       However, to spice things up, grab the first std::vector,
       last std::vector, and choose the others randomly.
    */
    TEUCHOS_TEST_FOR_EXCEPT(numvecs_2 != 5);
    ind[0] = 0;
    ind[1] = 5;
    ind[2] = 2;
    ind[3] = 2;
    ind[4] = 9;

    indContig[0] = 2;
    indContig[1] = 3;
    indContig[2] = 4;
    indContig[3] = 5;
    indContig[4] = 6;

    /*********** GetNumberVecs() *****************************************
       Verify:
       1) This number should be strictly positive
    *********************************************************************/
    if ( MVT::GetNumberVecs(*A) <= 0 ) {
      om->stream(Warnings)
        << "*** ERROR *** MultiVectorTraits::GetNumberVecs()." << endl
        << "Returned <= 0." << endl;
      return false;
    }


    /*********** GetGlobalLength() ***************************************
       Verify:
       1) This number should be strictly positive
    *********************************************************************/
    if ( MVT::GetGlobalLength(*A) <= 0 ) {
      om->stream(Warnings)
        << "*** ERROR *** MultiVectorTraitsExt::GetGlobalLength()" << endl
        << "Returned <= 0." << endl;
      return false;
    }


    /*********** Clone() and MvNorm() ************************************
       Verify:
       1) Clone() allows us to specify the number of vectors
       2) Clone() returns a multivector of the same dimension
       3) Vector norms shouldn't be negative
       4) MvNorm result std::vector should not be resized
    *********************************************************************/
    {
      Teuchos::RCP<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(2*numvecs);
      bool ResizeWarning = false;
      if ( MVT::GetNumberVecs(*B) != numvecs ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::Clone()." << endl
          << "Did not allocate requested number of vectors." << endl;
        return false;
      }
      if ( MVT::GetGlobalLength(*B) != MVT::GetGlobalLength(*A) ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::Clone()." << endl
          << "Did not allocate requested number of vectors." << endl;
        return false;
      }
      MVT::MvNorm(*B, norms);
      if ( norms.size() != 2*numvecs && ResizeWarning==false ) {
        om->stream(Warnings)
          << "*** WARNING *** MultiVecTraits::MvNorm()." << endl
          << "Method resized the output vector." << endl;
        ResizeWarning = true;
      }
      for (int i=0; i<numvecs; i++) {
        if ( norms[i] < zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::Clone()." << endl
            << "Vector had negative norm." << endl;
          return false;
        }
      }
    }


    /*********** MvRandom() and MvNorm() and MvInit() ********************
       Verify:
       1) Set vectors to zero
       2) Check that norm is zero
       3) Perform MvRandom.
       4) Verify that vectors aren't zero anymore
       5) Perform MvRandom again.
       6) Verify that std::vector norms are different than before

       Without knowing something about the random distribution,
       this is about the best that we can do, to make sure that MvRandom
       did at least *something*.

       Also, make sure std::vector norms aren't negative.
    *********************************************************************/
    {
      Teuchos::RCP<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(numvecs), norms2(numvecs);

      MVT::MvInit(*B);
      MVT::MvNorm(*B, norms);
      for (int i=0; i<numvecs; i++) {
        if ( norms[i] != zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvInit() "
            << "and MultiVecTraits::MvNorm()" << endl
            << "Supposedly zero vector has non-zero norm." << endl;
          return false;
        }
      }
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, norms);
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, norms2);
      for (int i=0; i<numvecs; i++) {
        if ( norms[i] == zero_mag || norms2[i] == zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvRandom()." << endl
            << "Random vector was empty (very unlikely)." << endl;
          return false;
        }
        else if ( norms[i] < zero_mag || norms2[i] < zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvRandom()." << endl
            << "Vector had negative norm." << endl;
          return false;
        }
        else if ( norms[i] == norms2[i] ) {
          om->stream(Warnings)
            << "*** ERROR *** MutliVecTraits::MvRandom()." << endl
            << "Vectors not random enough." << endl;
          return false;
        }
      }
    }


    /*********** MvRandom() and MvNorm() and MvScale() *******************
       Verify:
       1) Perform MvRandom.
       2) Verify that vectors aren't zero
       3) Set vectors to zero via MvScale
       4) Check that norm is zero
    *********************************************************************/
    {
      Teuchos::RCP<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(numvecs);

      MVT::MvRandom(*B);
      MVT::MvScale(*B,STS::zero());
      MVT::MvNorm(*B, norms);
      for (int i=0; i<numvecs; i++) {
        if ( norms[i] != zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvScale(alpha) "
            << "Supposedly zero vector has non-zero norm." << endl;
          return false;
        }
      }

      MVT::MvRandom(*B);
      std::vector<ScalarType> zeros(numvecs,STS::zero());
      MVT::MvScale(*B,zeros);
      MVT::MvNorm(*B, norms);
      for (int i=0; i<numvecs; i++) {
        if ( norms[i] != zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvScale(alphas) "
            << "Supposedly zero vector has non-zero norm." << endl;
          return false;
        }
      }
    }


    /*********** MvInit() and MvNorm() ***********************************
       A std::vector of ones of dimension n should have norm std::sqrt(n)
       1) Init vectors to all ones
       2) Verify that norm is std::sqrt(n)
       3) Verify that norms aren't negative

       Note: I'm not sure that we can expect this to hold in practice.
              Maybe something like std::abs(norm-std::sqrt(n)) < STS::eps()  ???
              The sum of 1^2==1 should be n, but what about std::sqrt(n)?
              They may be using a different square root than ScalartTraits
              On my iBook G4 and on jeter, this test works.
              Right now, this has been demoted to a warning.
    *********************************************************************/
    {
      Teuchos::RCP<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(numvecs);

      MVT::MvInit(*B,one);
      MVT::MvNorm(*B, norms);
      bool BadNormWarning = false;
      for (int i=0; i<numvecs; i++) {
        if ( norms[i] < zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvRandom()." << endl
            << "Vector had negative norm." << endl;
          return false;
        }
        else if ( norms[i] != STS::squareroot(MVT::GetGlobalLength(*B)) && !BadNormWarning ) {
          om->stream(Warnings)
            << endl
            << "Warning testing MultiVecTraits::MvInit()." << endl
            << "Ones std::vector should have norm std::sqrt(dim)." << endl
            << "norms[i]: " << norms[i] << "\tdim: " << MVT::GetGlobalLength(*B) << endl << endl;
          BadNormWarning = true;
        }
      }
    }


    /*********** MvInit() and MvNorm() ***********************************
       A std::vector of zeros of dimension n should have norm 0
       1) Verify that norms aren't negative
       2) Verify that norms are zero

       We must know this works before the next tests.
    *********************************************************************/
    {
      Teuchos::RCP<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(numvecs);
      MVT::MvInit(*B, zero_mag);
      MVT::MvNorm(*B, norms);
      for (int i=0; i<numvecs; i++) {
        if ( norms[i] < zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvInit()." << endl
            << "Vector had negative norm." << endl;
          return false;
        }
        else if ( norms[i] != zero_mag ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvInit()." << endl
            << "Zero std::vector should have norm zero." << endl;
          return false;
        }
      }
    }


    /*********** CloneCopy(MV,std::vector<int>) and MvNorm ********************
       1) Check quantity/length of vectors
       2) Check std::vector norms for agreement
       3) Zero out B and make sure that C norms are not affected
    *********************************************************************/
    {
      Teuchos::RCP<MV> B, C;
      std::vector<MagType> norms(numvecs), norms2(numvecs);

      B = MVT::Clone(*A,numvecs);
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, norms);
      C = MVT::CloneCopy(*B,ind);
      MVT::MvNorm(*C, norms2);
      if ( MVT::GetNumberVecs(*C) != numvecs_2 ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::CloneCopy(ind)." << endl
          << "Wrong number of vectors." << endl;
        return false;
      }
      if ( MVT::GetGlobalLength(*C) != MVT::GetGlobalLength(*B) ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::CloneCopy(ind)." << endl
          << "Vector lengths don't match." << endl;
        return false;
      }
      for (int i=0; i<numvecs_2; i++) {
        if (STS::magnitude (norms2[i] - norms[ind[i]]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::CloneCopy(ind)." << endl
            << "Copied vectors do not agree: "
            << norms2[i] << " != " << norms[ind[i]] << endl
            << "Difference " << STS::magnitude (norms2[i] - norms[ind[i]])
            << " exceeds the tolerance 100*eps = " << tol << endl;
          return false;
        }
      }
      MVT::MvInit(*B,zero);
      MVT::MvNorm(*C, norms);
      for (int i=0; i<numvecs_2; i++) {
        //if ( norms2[i] != norms[i] ) {
        if (STS::magnitude (norms2[i] - norms[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::CloneCopy(ind)." << endl
            << "Copied vectors were not independent." << endl
            << norms2[i] << " != " << norms[i] << endl
            << "Difference " << STS::magnitude (norms2[i] - norms[i])
            << " exceeds the tolerance 100*eps = " << tol << endl;
          return false;
        }
      }
    }

    /*********** CloneCopy(MV) and MvNorm ********************************
       1) Check quantity
       2) Check value of norms
       3) Zero out B and make sure that C is still okay
    *********************************************************************/
    {
      Teuchos::RCP<MV> B, C;
      std::vector<MagType> norms(numvecs), norms2(numvecs);

      B = MVT::Clone(*A,numvecs);
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, norms);
      C = MVT::CloneCopy(*B);
      MVT::MvNorm(*C, norms2);
      if ( MVT::GetNumberVecs(*C) != numvecs ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::CloneCopy()." << endl
          << "Wrong number of vectors." << endl;
        return false;
      }
      for (int i=0; i<numvecs; i++) {
        if (STS::magnitude (norms2[i] - norms[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::CloneCopy()." << endl
            << "Copied vectors do not agree: "
            << norms2[i] << " != " << norms[i] << endl
            << "Difference " << STS::magnitude (norms2[i] - norms[i])
            << " exceeds the tolerance 100*eps = " << tol << endl;
          return false;
        }
      }
      MVT::MvInit(*B,zero);
      MVT::MvNorm(*C, norms);
      for (int i=0; i<numvecs; i++) {
        //if ( norms2[i] != norms[i] ) {
        if (STS::magnitude (norms2[i] - norms[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::CloneCopy()." << endl
            << "Copied vectors were not independent." << endl
            << norms2[i] << " != " << norms[i] << endl
            << "Difference " << STS::magnitude (norms2[i] - norms[i])
            << " exceeds the tolerance 100*eps = " << tol << endl;
          return false;
        }
      }
    }


    /*********** CloneView(MV,std::vector<int>) and MvNorm ********************
       Check that we have a view of the selected vectors
       1) Check quantity
       2) Check value of norms
       3) Zero out B and make sure that C is zero as well
    *********************************************************************/
    {
      Teuchos::RCP<MV> B, C;
      std::vector<MagType> norms(numvecs), norms2(numvecs);

      B = MVT::Clone(*A,numvecs);
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, norms);
      C = MVT::CloneViewNonConst(*B,indContig);
      MVT::MvNorm(*C, norms2);
      if ( MVT::GetNumberVecs(*C) != numvecs_2 ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::CloneView(indContig)." << endl
          << "Wrong number of vectors." << endl;
        return false;
      }
      for (int i=0; i<numvecs_2; i++) {
        //if ( norms2[i] != norms[indContig[i]] ) {
        if (STS::magnitude (norms2[i] - norms[indContig[i]]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::CloneView(indContig)." << endl
            << "Viewed vectors do not agree." << endl;
          return false;
        }
      }
    }


    /*********** CloneView(const MV,std::vector<int>) and MvNorm() ************
       Check that we have a view of the selected vectors.
       1) Check quantity
       2) Check value of norms for agreement
       3) Zero out B and make sure that C is zerod as well
    *********************************************************************/
   {
      Teuchos::RCP<MV> B;
      Teuchos::RCP<const MV> C;
      std::vector<MagType> normsB(numvecs), normsC(numvecs_2);
      std::vector<int> allind(numvecs);
      for (int i=0; i<numvecs; i++) {
        allind[i] = i;
      }

      B = MVT::Clone(*A,numvecs);
      MVT::MvRandom( *B );
      MVT::MvNorm(*B, normsB);
      C = MVT::CloneView(*B,indContig);
      MVT::MvNorm(*C, normsC);
      if ( MVT::GetNumberVecs(*C) != numvecs_2 ) {
        om->stream(Warnings)
          << "*** ERROR *** const MultiVecTraits::CloneView(indContig)." << endl
          << "Wrong number of vectors." << endl;
        return false;
      }
      for (int i=0; i<numvecs_2; i++) {
        //if ( normsC[i] != normsB[ind[i]] ) {
        if (STS::magnitude (normsC[i] - normsB[indContig[i]]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** const MultiVecTraits::CloneView(indContig)." << endl
            << "Viewed vectors do not agree." << endl;
          return false;
        }
      }
    }

    /*********** SetBlock() and MvNorm() *********************************
       SetBlock() will copy the vectors from C into B
       1) Verify that the specified vectors were copied
       2) Verify that the other vectors were not modified
       3) Verify that C was not modified
       4) Change C and then check B to make sure it was not modified

       Use a different index set than has been used so far (distinct entries).
       This is because duplicate entries will cause the std::vector to be
       overwritten, making it more difficult to test.
    *********************************************************************/
    {
      Teuchos::RCP<MV> B, C; //numvecs = 10, numvecs_2 = 5
      std::vector<MagType> normsB1(numvecs), normsB2(numvecs),
                           normsC1(numvecs_2), normsC2(numvecs_2);

      B = MVT::Clone(*A,numvecs);
      C = MVT::Clone(*A,numvecs_2);
      // Just do every other one, interleaving the vectors of C into B
      ind.resize(numvecs_2);
      for (int i=0; i<numvecs_2; i++) {
        ind[i] = 2*i; // 0,2,4,6,..
      }
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);

      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      MVT::SetBlock(*C,ind,*B);// Copy vectors of C to vectors specified by ind in B. 
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);

      // check that C was not changed by SetBlock
      for (int i=0; i<numvecs_2; i++) {
        //if ( normsC1[i] != normsC2[i] ) {
        if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
            << "Operation modified source vectors." << endl;
          return false;
        }
      }
      // check that the correct vectors of B were modified
      // and the others were not
      for (int i=0; i<numvecs; i++) {
        if (i % 2 == 0) {
          // should be a vector from C
          if (STS::magnitude (normsB2[i] - normsC1[i/2]) > tol) {
            om->stream(Warnings)
              << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
              << "Copied vectors do not agree: " << endl
              << normsB2[i] << " != " << normsC1[i/2] << endl
              << "Difference " << STS::magnitude (normsB2[i] - normsC1[i/2])
              << " exceeds the tolerance 100*eps = " << tol << endl;
            return false;
          }
        }
        else {
          // should be an original vector
          if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
            om->stream(Warnings)
              << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
              << "Incorrect vectors were modified." << endl
              << normsB1[i] << " != " << normsB2[i] << endl
              << "Difference " << STS::magnitude (normsB2[i] - normsB2[i])
              << " exceeds the tolerance 100*eps = " << tol << endl;
            return false;
          }
        }
      }
      MVT::MvInit(*C,zero);
      MVT::MvNorm(*B,normsB1);
      // verify that we copied and didn't reference
      for (int i=0; i<numvecs; i++) {
        //if ( normsB1[i] != normsB2[i] ) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
            << "Copied vectors were not independent." << endl
            << normsB1[i] << " != " << normsB2[i] << endl
            << "Difference " << STS::magnitude (normsB1[i] - normsB2[i])
            << " exceeds the tolerance 100*eps = " << tol << endl;
          return false;
        }
      }
    }


    /*********** SetBlock() and MvNorm() *********************************
       SetBlock() will copy the vectors from C into B
       1) Verify that the specified vectors were copied
       2) Verify that the other vectors were not modified
       3) Verify that C was not modified
       4) Change C and then check B to make sure it was not modified

       Use a different index set than has been used so far (distinct entries).
       This is because duplicate entries will cause the std::vector to be
       overwritten, making it more difficult to test.

       These tests are the same as the ones above, except that the
       number of indices (to be copied into B) is less than the number
       of vectors in C, so that not all of C is put into B.
    *********************************************************************/
    {
      Teuchos::RCP<MV> B, C;
      // set these: we assume below that setSize*2=BSize
      const int CSize   = 6,
                setSize = 5,
                BSize   = 2*setSize;
      std::vector<MagType> normsB1(BSize), normsB2(BSize),
                           normsC1(CSize), normsC2(CSize);

      B = MVT::Clone(*A,BSize);
      C = MVT::Clone(*A,CSize);
      // Just do every other one, interleaving the vectors of C into B
      ind.resize(setSize);
      for (int i=0; i<setSize; i++) {
        ind[i] = 2*i;
      }
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);

      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      MVT::SetBlock(*C,ind,*B);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);

      // check that C was not changed by SetBlock
      for (int i=0; i<CSize; i++) {
        //if ( normsC1[i] != normsC2[i] ) {
        if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
            << "Operation modified source vectors." << endl;
          return false;
        }
      }
      // check that the correct vectors of B were modified
      // and the others were not
      for (int i=0; i<BSize; i++) {
        if (i % 2 == 0) {
          // should be a vector from C
          const MagType diff = STS::magnitude (normsB2[i] - normsC1[i/2]);
          if (diff > tol) {
            om->stream(Warnings)
              << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
              << "Copied vectors do not agree: " << endl
              << normsB2[i] << " != " << normsC1[i/2] << endl
              << "Difference " << diff << " exceeds the tolerance 100*eps = "
              << tol << endl;
            return false;
          }
        }
        else {
          // should be an original vector
          const MagType diff = STS::magnitude (normsB1[i] - normsB2[i]);
          //if ( normsB1[i] != normsB2[i] ) {
          if (diff > tol) {
            om->stream(Warnings)
              << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
              << "Incorrect vectors were modified." << endl
              << normsB1[i] << " != " << normsB2[i] << endl
              << "Difference " << diff << " exceeds the tolerance 100*eps = "
              << tol << endl;
            return false;
          }
        }
      }
      MVT::MvInit(*C,zero);
      MVT::MvNorm(*B,normsB1);
      // verify that we copied and didn't reference
      for (int i=0; i<numvecs; i++) {
        //if ( normsB1[i] != normsB2[i] ) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
            << "Copied vectors were not independent." << endl
            << normsB1[i] << " != " << normsB2[i] << endl
            << "Difference " << STS::magnitude (normsB1[i] - normsB2[i])
            << " exceeds the tolerance 100*eps = " << tol << endl;
          return false;
        }
      }
    }


    /*********** MvTransMv() *********************************************
      Performs C = alpha * A^H * B, where
        alpha is type ScalarType
        A,B are type MV with p and q vectors, respectively
        C is a SerialDenseMatrix<int,ScalarType> ALREADY sized to p by q

        Verify:
        1) C is not resized by the routine
        3) Check that zero*(A^H B) == zero
        3) Check inner product inequality:
                                        [ |a1|*|b1|    ...    |ap|*|b1| ]
           [a1 ... ap]^H [b1 ... bq] <= [    ...    |ai|*|bj|    ...    ]
                                        [ |ap|*|b1|    ...    |ap|*|bq| ]
        4) Zero B and check that C is zero
        5) Zero A and check that C is zero

        Note: Should we really require that C is correctly sized already?
        Epetra does (and crashes if it isn't.)
    *********************************************************************/
    {
      const int p = 7;
      const int q = 9;
      Teuchos::RCP<MV> B, C;
      std::vector<MagType> normsB(p), normsC(q);
      Teuchos::SerialDenseMatrix<int,ScalarType> SDM(p,q);

      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,q);

      // randomize the multivectors
      MVT::MvRandom(*B);
      MVT::MvNorm(*B,normsB);
      MVT::MvRandom(*C);
      MVT::MvNorm(*C,normsC);

      // perform SDM  = zero() * B^H * C
      MVT::MvTransMv( zero, *B, *C, SDM );

      // check the sizes: not allowed to have shrunk
      if ( SDM.numRows() != p || SDM.numCols() != q ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
          << "Routine resized SerialDenseMatrix." << endl;
        return false;
      }

      // check that zero**A^H*B == zero
      if ( SDM.normOne() != zero ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
          << "Scalar argument processed incorrectly." << endl;
        return false;
      }

      // perform SDM  = one * B^H * C
      MVT::MvTransMv( one, *B, *C, SDM );

      // check the norms: a^H b = |a| |b| cos(theta) <= |a| |b|
      // with equality only when a and b are colinear
      for (int i=0; i<p; i++) {
        for (int j=0; j<q; j++) {
          if (   STS::magnitude(SDM(i,j))
               > STS::magnitude(normsB[i]*normsC[j]) ) {
            om->stream(Warnings)
              << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
              << "Triangle inequality did not hold: "
              << STS::magnitude(SDM(i,j))
              << " > "
              << STS::magnitude(normsB[i]*normsC[j])
              << endl;
            return false;
          }
        }
      }
      MVT::MvInit(*C);
      MVT::MvRandom(*B);
      MVT::MvTransMv( one, *B, *C, SDM );
      for (int i=0; i<p; i++) {
        for (int j=0; j<q; j++) {
          if ( SDM(i,j) != zero ) {
            om->stream(Warnings)
              << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
              << "Inner products not zero for C==0." << endl;
            return false;
          }
        }
      }
      MVT::MvInit(*B);
      MVT::MvRandom(*C);
      MVT::MvTransMv( one, *B, *C, SDM );
      for (int i=0; i<p; i++) {
        for (int j=0; j<q; j++) {
          if ( SDM(i,j) != zero ) {
            om->stream(Warnings)
              << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
              << "Inner products not zero for B==0." << endl;
            return false;
          }
        }
      }
    }


    /*********** MvDot() *************************************************
        Verify:
        1) Results std::vector not resized
        2) Inner product inequalities are satisfied
        3) Zero vectors give zero inner products
    *********************************************************************/
    {
      const int p = 7;
      const int q = 9;
      Teuchos::RCP<MV> B, C;
      std::vector<ScalarType> iprods(p+q);
      std::vector<MagType> normsB(numvecs), normsC(numvecs);

      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,p);

      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB);
      MVT::MvNorm(*C,normsC);
      MVT::MvDot( *B, *C, iprods );
      if ( iprods.size() != p+q ) {
        om->stream(Warnings)
          << "*** ERROR *** MultiVecTraits::MvDot." << endl
          << "Routine resized results std::vector." << endl;
        return false;
      }
      for (int i=0; i<BELOS_MIN(p,q); i++) {
        if ( STS::magnitude(iprods[i])
             > STS::magnitude(normsB[i]*normsC[i]) ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvDot()." << endl
            << "Inner products not valid." << endl;
          return false;
        }
      }
      MVT::MvInit(*B);
      MVT::MvRandom(*C);
      MVT::MvDot( *B, *C, iprods );
      for (int i=0; i<p; i++) {
        if ( iprods[i] != zero ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvDot()." << endl
            << "Inner products not zero for B==0." << endl;
          return false;
        }
      }
      MVT::MvInit(*C);
      MVT::MvRandom(*B);
      MVT::MvDot( *B, *C, iprods );
      for (int i=0; i<p; i++) {
        if ( iprods[i] != zero ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvDot()." << endl
            << "Inner products not zero for C==0." << endl;
          return false;
        }
      }
    }


    /*********** MvAddMv() ***********************************************
       D = alpha*B + beta*C
       1) Use alpha==0,beta==1 and check that D == C
       2) Use alpha==1,beta==0 and check that D == B
       3) Use D==0 and D!=0 and check that result is the same
       4) Check that input arguments are not modified
    *********************************************************************/
    {
      const int p = 7;
      Teuchos::RCP<MV> B, C, D;
      std::vector<MagType> normsB1(p), normsB2(p),
                           normsC1(p), normsC2(p),
                           normsD1(p), normsD2(p);

      Teuchos::SerialDenseMatrix<int,ScalarType> Alpha(1,1), Beta(1,1);
      Teuchos::randomSyncedMatrix( Alpha );
      Teuchos::randomSyncedMatrix( Beta );
      ScalarType alpha = Alpha(0,0),
                  beta = Beta(0,0);

      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,p);
      D = MVT::Clone(*A,p);

      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);

      // check that 0*B+1*C == C
      MVT::MvAddMv(zero,*B,one,*C,*D);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      MVT::MvNorm(*D,normsD1);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Input arguments were modified." << endl;
          return false;
        }
        else if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Input arguments were modified." << endl;
          return false;
        }
        else if (STS::magnitude (normsC1[i] - normsD1[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Assignment did not work." << endl;
          return false;
        }
      }

      // check that 1*B+0*C == B
      MVT::MvAddMv(one,*B,zero,*C,*D);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      MVT::MvNorm(*D,normsD1);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Input arguments were modified." << endl;
          return false;
        }
        else if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Input arguments were modified." << endl;
          return false;
        }
        else if (STS::magnitude (normsB1[i] - normsD1[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Assignment did not work." << endl;
          return false;
        }
      }

      // check that alpha*B+beta*C -> D is invariant under initial D
      // first, try random D
      MVT::MvRandom(*D);
      MVT::MvAddMv(alpha,*B,beta,*C,*D);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      MVT::MvNorm(*D,normsD1);
      // check that input args are not modified
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Input arguments were modified." << endl;
          return false;
        }
        else if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Input arguments were modified." << endl;
          return false;
        }
      }
      // next, try zero D
      MVT::MvInit(*D);
      MVT::MvAddMv(alpha,*B,beta,*C,*D);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      MVT::MvNorm(*D,normsD2);
      // check that input args are not modified and that D is the same
      // as the above test
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Input arguments were modified." << endl;
          return false;
        }
        else if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Input arguments were modified." << endl;
          return false;
        }
        else if (STS::magnitude (normsD1[i] - normsD2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
            << "Results varies depending on initial state of dest vectors." << endl;
          return false;
        }
      }
    }

    /*********** MvAddMv() ***********************************************
       Similar to above, but where B or C are potentially the same
       object as D. This case is commonly used, for example, to affect
       A <- alpha*A
       via
       MvAddMv(alpha,A,zero,A,A)
          ** OR **
       MvAddMv(zero,A,alpha,A,A)

       The result is that the operation has to be "atomic". That is,
       B and C are no longer reliable after D is modified, so that
       the assignment to D must be the last thing to occur.

       D = alpha*B + beta*C

       1) Use alpha==0,beta==1 and check that D == C
       2) Use alpha==1,beta==0 and check that D == B
    *********************************************************************/
    {
      const int p = 7;
      Teuchos::RCP<MV> B, D;
      Teuchos::RCP<const MV> C;
      std::vector<MagType> normsB(p),
                           normsD(p);
      std::vector<int> lclindex(p);
      for (int i=0; i<p; i++) lclindex[i] = i;

      B = MVT::Clone(*A,p);
      C = MVT::CloneView(*B,lclindex);
      D = MVT::CloneViewNonConst(*B,lclindex);

      MVT::MvRandom(*B);
      MVT::MvNorm(*B,normsB);

      // check that 0*B+1*C == C
      MVT::MvAddMv(zero,*B,one,*C,*D);
      MVT::MvNorm(*D,normsD);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB[i] - normsD[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv() #2" << endl
            << "Assignment did not work." << endl;
          return false;
        }
      }

      // check that 1*B+0*C == B
      MVT::MvAddMv(one,*B,zero,*C,*D);
      MVT::MvNorm(*D,normsD);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB[i] - normsD[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvAddMv() #2" << endl
            << "Assignment did not work." << endl;
          return false;
        }
      }

    }

    /*********** MvTimesMatAddMv() 7 by 5 ********************************
       C = alpha*B*SDM + beta*C
       1) Use alpha==0, SDM!=0, beta==1 and check that C is unchanged
       2) Use alpha==0, SDM!=0, beta==0 and check that C is set to zero
       3) Use alpha==1, SDM==I, beta==0 and check that C is set to B
       4) Use alpha==1, SDM==0, beta==1 and check that C is unchanged
       5) Test with non-square matrices
       6) Always check that input arguments are not modified
    *********************************************************************/
    {
      const int p = 7, q = 5;
      Teuchos::RCP<MV> B, C;
      Teuchos::SerialDenseMatrix<int,ScalarType> SDM(p,q);
      std::vector<MagType> normsC1(q), normsC2(q),
                           normsB1(p), normsB2(p);

      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,q);

      // Test 1: alpha==0, SDM!=0, beta==1 and check that C is unchanged
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      Teuchos::randomSyncedMatrix(SDM);
      MVT::MvTimesMatAddMv(zero,*B,SDM,one,*C);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Input vectors were modified." << endl;
          return false;
        }
      }
      for (int i=0; i<q; i++) {
        if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 1 failed." << endl;
          return false;
        }
      }

      // Test 2: alpha==0, SDM!=0, beta==0 and check that C is set to zero
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      Teuchos::randomSyncedMatrix(SDM);
      MVT::MvTimesMatAddMv(zero,*B,SDM,zero,*C);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Input vectors were modified." << endl;
          return false;
        }
      }
      for (int i=0; i<q; i++) {
        if ( normsC2[i] != zero ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 2 failed: "
            << normsC2[i]
            << " != "
            << zero
            << endl;
          return false;
        }
      }

      // Test 3: alpha==1, SDM==|I|, beta==0 and check that C is set to B
      //                        |0|
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      SDM.scale(zero);
      for (int i=0; i<q; i++) {
        SDM(i,i) = one;
      }
      MVT::MvTimesMatAddMv(one,*B,SDM,zero,*C);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Input vectors were modified." << endl;
          return false;
        }
      }
      for (int i=0; i<q; i++) {
        if (STS::magnitude (normsB1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 3 failed: "
            << normsB1[i]
            << " != "
            << normsC2[i]
            << endl;
          return false;
        }
      }

      // Test 4: alpha==1, SDM==0, beta==1 and check that C is unchanged
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      SDM.scale(zero);
      MVT::MvTimesMatAddMv(one,*B,SDM,one,*C);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Input vectors were modified." << endl;
          return false;
        }
      }
      for (int i=0; i<q; i++) {
        if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 4 failed." << endl;
          return false;
        }
      }
    }

    /*********** MvTimesMatAddMv() 5 by 7 ********************************
       C = alpha*B*SDM + beta*C
       1) Use alpha==0, SDM!=0, beta==1 and check that C is unchanged
       2) Use alpha==0, SDM!=0, beta==0 and check that C is set to zero
       3) Use alpha==1, SDM==I, beta==0 and check that C is set to B
       4) Use alpha==1, SDM==0, beta==1 and check that C is unchanged
       5) Test with non-square matrices
       6) Always check that input arguments are not modified
    *********************************************************************/
    {
      const int p = 5, q = 7;
      Teuchos::RCP<MV> B, C;
      Teuchos::SerialDenseMatrix<int,ScalarType> SDM(p,q);
      std::vector<MagType> normsC1(q), normsC2(q),
                           normsB1(p), normsB2(p);

      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,q);

      // Test 5: alpha==0, SDM!=0, beta==1 and check that C is unchanged
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      Teuchos::randomSyncedMatrix(SDM);
      MVT::MvTimesMatAddMv(zero,*B,SDM,one,*C);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Input vectors were modified." << endl;
          return false;
        }
      }
      for (int i=0; i<q; i++) {
        if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 5 failed." << endl;
          return false;
        }
      }

      // Test 6: alpha==0, SDM!=0, beta==0 and check that C is set to zero
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      Teuchos::randomSyncedMatrix(SDM);
      MVT::MvTimesMatAddMv(zero,*B,SDM,zero,*C);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Input vectors were modified." << endl;
          return false;
        }
      }
      for (int i=0; i<q; i++) {
        if ( normsC2[i] != zero ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 6 failed: "
            << normsC2[i]
            << " != "
            << zero
            << endl;
          return false;
        }
      }

      // Test 7: alpha==1, SDM==[I 0], beta==0 and check that C is set to B
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      SDM.scale(zero);
      for (int i=0; i<p; i++) {
        SDM(i,i) = one;
      }
      MVT::MvTimesMatAddMv(one,*B,SDM,zero,*C);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Input vectors were modified." << endl;
          return false;
        }
      }
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 7 failed." << endl;
          return false;
        }
      }
      for (int i=p; i<q; i++) {
        if ( normsC2[i] != zero ) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 7 failed." << endl;
          return false;
        }
      }

      // Test 8: alpha==1, SDM==0, beta==1 and check that C is unchanged
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,normsB1);
      MVT::MvNorm(*C,normsC1);
      SDM.scale(zero);
      MVT::MvTimesMatAddMv(one,*B,SDM,one,*C);
      MVT::MvNorm(*B,normsB2);
      MVT::MvNorm(*C,normsC2);
      for (int i=0; i<p; i++) {
        if (STS::magnitude (normsB1[i] - normsB2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Input vectors were modified." << endl;
          return false;
        }
      }
      for (int i=0; i<q; i++) {
        if (STS::magnitude (normsC1[i] - normsC2[i]) > tol) {
          om->stream(Warnings)
            << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
            << "Arithmetic test 8 failed." << endl;
          return false;
        }
      }
    }

    return true;

  }



  /// \brief Test correctness of OperatorTraits specialization and its
  ///   operator implementation.
  ///
  /// \tparam ScalarType The type of the entries in the multivectors;
  ///   the first template parameter of MultiVecTraits and
  ///   OperatorTraits.
  /// \tparam MV The multivector type; the second template parameter
  ///   of MultiVecTraits and OperatorTraits.
  /// \tparam OP The operator type; the third template parameter
  ///   of OperatorTraits.
  ///
  /// \param om [in/out] A valid OutputManager, for displaying test results.
  ///
  /// \param A [in] An initial multivector, to use for making new
  ///   multivectors.  (Belos doesn't currently have a "vector space"
  ///   abstraction; making a new multivector requires a valid input
  ///   multivector to clone.)
  ///
  /// \param M [in] The operator to test.
  ///
  /// \return Test result: true if all tests passed, else false.
  template< class ScalarType, class MV, class OP>
  bool
  TestKokkosOperatorTraits (const Teuchos::RCP<OutputManager<ScalarType> > &om,
                      const Teuchos::RCP<const MV> &A,
                      const Teuchos::RCP<const OP> &M)
  {
    using Teuchos::SetScientific;
    using std::endl;
    typedef MultiVecTraits<ScalarType, MV>    MVT;
    typedef Teuchos::ScalarTraits<ScalarType> STS;
    typedef typename STS::magnitudeType       MagType;

    // Make sure that all floating-point numbers are printed with the
    // right precision.
    SetScientific<ScalarType> sci (om->stream (Warnings));

    // FIXME (mfh 09 Jan 2013) Added an arbitrary tolerance in case
    // norms are not computed deterministically (which is possible
    // even with MPI only, and more likely with threads).
    const MagType tol = Teuchos::as<MagType> (100) * STS::eps ();

    /* OPT Contract:
       Apply()
         MV: OP*zero == zero
             Warn if OP is not deterministic (OP*A != OP*A)
             Does not modify input arguments
    *********************************************************************/

    typedef MultiVecTraits<ScalarType, MV>     MVT;
    typedef Teuchos::ScalarTraits<ScalarType>  STS;
    typedef OperatorTraits<ScalarType, MV, OP> OPT;
    typedef typename STS::magnitudeType        MagType;

    const int numvecs = 10;

    Teuchos::RCP<MV> B = MVT::Clone(*A,numvecs),
                             C = MVT::Clone(*A,numvecs);

    std::vector<MagType> normsB1(numvecs), normsB2(numvecs),
                         normsC1(numvecs), normsC2(numvecs);
    bool NonDeterministicWarning;

    /*********** Apply() *************************************************
        Verify:
        1) OP*B == OP*B; OP is deterministic (just warn on this)
        2) OP*zero == 0
        3) OP*B doesn't modify B
        4) OP*B is invariant under initial state of destination vectors
    *********************************************************************/
    MVT::MvInit(*B);
    MVT::MvRandom(*C);
    MVT::MvNorm(*B,normsB1);
    OPT::Apply(*M,*B,*C);
    MVT::MvNorm(*B,normsB2);
    MVT::MvNorm(*C,normsC2);
    for (int i=0; i<numvecs; i++) {
      if (STS::magnitude (normsB2[i] - normsB1[i]) > tol) {
        om->stream(Warnings)
          << "*** ERROR *** OperatorTraits::Apply() [1]" << endl
          << "Apply() modified the input vectors." << endl
          << "Original: " << normsB1[i] << "; After: " << normsB2[i] << endl
          << "Difference " << STS::magnitude (normsB2[i] - normsB1[i])
          << " exceeds the tolerance 100*eps = " << tol << endl;
        return false;
      }
      if (normsC2[i] != STS::zero()) {
        om->stream(Warnings)
          << "*** ERROR *** OperatorTraits::Apply() [1]" << endl
          << "Operator applied to zero did not return zero." << endl;
        return false;
      }
    }

    // If we send in a random matrix, we should not get a zero return
    MVT::MvRandom(*B);
    MVT::MvNorm(*B,normsB1);
    OPT::Apply(*M,*B,*C);
    MVT::MvNorm(*B,normsB2);
    MVT::MvNorm(*C,normsC2);
    bool ZeroWarning = false;
    for (int i=0; i<numvecs; i++) {
      if (STS::magnitude (normsB2[i] - normsB1[i]) > tol) {
        om->stream(Warnings)
          << "*** ERROR *** OperatorTraits::Apply() [2]" << endl
          << "Apply() modified the input vectors." << endl
          << "Original: " << normsB1[i] << "; After: " << normsB2[i] << endl
          << "Difference " << STS::magnitude (normsB2[i] - normsB1[i])
          << " exceeds the tolerance 100*eps = " << tol << endl;
        return false;
      }
      if (normsC2[i] == STS::zero() && ZeroWarning==false ) {
        om->stream(Warnings)
          << "*** ERROR *** OperatorTraits::Apply() [2]" << endl
          << "Operator applied to random vectors returned zero." << endl;
        ZeroWarning = true;
      }
    }

    // Apply operator with C init'd to zero
    MVT::MvRandom(*B);
    MVT::MvNorm(*B,normsB1);
    MVT::MvInit(*C);
    OPT::Apply(*M,*B,*C);
    MVT::MvNorm(*B,normsB2);
    MVT::MvNorm(*C,normsC1);
    for (int i=0; i<numvecs; i++) {
      if (STS::magnitude (normsB2[i] - normsB1[i]) > tol) {
        om->stream(Warnings)
          << "*** ERROR *** OperatorTraits::Apply() [3]" << endl
          << "Apply() modified the input vectors." << endl
          << "Original: " << normsB1[i] << "; After: " << normsB2[i] << endl
          << "Difference " << STS::magnitude (normsB2[i] - normsB1[i])
          << " exceeds the tolerance 100*eps = " << tol << endl;
        return false;
      }
    }

    // Apply operator with C init'd to random
    // Check that result is the same as before; warn if not.
    // This could be a result of a bug, or a stochastic
    //   operator. We do not want to prejudice against a
    //   stochastic operator.
    MVT::MvRandom(*C);
    OPT::Apply(*M,*B,*C);
    MVT::MvNorm(*B,normsB2);
    MVT::MvNorm(*C,normsC2);
    NonDeterministicWarning = false;
    for (int i=0; i<numvecs; i++) {
      if (STS::magnitude (normsB2[i] - normsB1[i]) > tol) {
        om->stream(Warnings)
          << "*** ERROR *** OperatorTraits::Apply() [4]" << endl
          << "Apply() modified the input vectors." << endl
          << "Original: " << normsB1[i] << "; After: " << normsB2[i] << endl
          << "Difference " << STS::magnitude (normsB2[i] - normsB1[i])
          << " exceeds the tolerance 100*eps = " << tol << endl;
        return false;
      }
      if (normsC1[i] != normsC2[i] && !NonDeterministicWarning) {
        om->stream(Warnings)
          << endl
          << "*** WARNING *** OperatorTraits::Apply() [4]" << endl
          << "Apply() returned two different results." << endl << endl;
        NonDeterministicWarning = true;
      }
    }

    return true;

  }

}

#endif
