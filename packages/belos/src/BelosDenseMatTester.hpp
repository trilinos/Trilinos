//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
//@HEADER
//
#ifndef BELOS_DENSETRAITS_TESTER_HPP
#define BELOS_DENSETRAITS_TESTER_HPP

/*! \file BelosDenseTraitsTester.hpp
  \brief Test routines for MultiVecTraits and OperatorTraits conformity.
*/

#include "BelosConfigDefs.hpp"//TODO is this needed?
#include "BelosTypes.hpp"//TODO is this needed?

#include "BelosOutputManager.hpp"
#include "BelosDenseMatTraits.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_SetScientific.hpp"//TODO is this needed?
#include "Teuchos_SerialDenseMatrix.hpp"

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
  template< class ScalarType, class DM = Teuchos::SerialDenseMatrix<int,ScalarType>>
  bool
  TestDenseMatTraits (const Teuchos::RCP<OutputManager<ScalarType> > &om)
  {
    using Teuchos::SetScientific;
    using Teuchos::RCP;
    using std::endl;
    typedef DenseMatTraits<ScalarType, DM>    DMT;
    typedef Teuchos::ScalarTraits<ScalarType> STS;
    typedef typename STS::magnitudeType       MagType;

    // Make sure that all floating-point numbers are printed with the
    // right precision.
    SetScientific<ScalarType> sci (om->stream (Warnings));

    // Arbitrary tolerance in case
    // norms are not computed deterministically (which is possible
    // even with MPI only, and more likely with threads).
    const MagType tol = Teuchos::as<MagType> (100) * STS::eps ();

    /* Dense Traits Contract:

         Clone(MV,int)
         CloneCopy(MV)
         CloneCopy(MV,vector<int>)
           USER: will request positive number of vectors
             MV: will return a multivector with exactly the number of
                   requested vectors.
                 vectors are the same dimension as the cloned MV


    *********************************************************************/

    const ScalarType one      = STS::one();
    const ScalarType zero     = STS::zero();
    const MagType    zero_mag = Teuchos::ScalarTraits<MagType>::zero();

    /*********** Basic Functions: *****************************************
       Verify:
       1) Can call DM constructor.
       2) Number of rows and cols matches requested.
       3) Init to zero if requested. 
       4) Call constructor with no init to zero.
       5) Can randomize matrix values. 
       6) Can change a matrix value.
       7) Const value access is const.
       8) 
       Basic tester: Just call all the functions. 
    *********************************************************************/
    { //begin test scope
      int numrows = 4;
      int numcols = 3;
      RCP<DM> dm1 = DMT::Create(numrows,numcols);

      //Check number of rows and cols.
      int checknumrows = DMT::GetNumRows(*dm1);
      int checknumcols = DMT::GetNumCols(*dm1);
      if( checknumrows != numrows || checknumcols != numcols){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::GetNumRows or GetNumCols" << endl
          << "Returned incorrect value." << endl;
        return false;
      }
      
      //Check init to zero on create.
      for(int i = 0; i<numrows; i++){
        for(int j = 0; j<numcols; j++){
          ScalarType val = DMT::Value(*dm1,i,j);
          if(val != zero){
            om->stream(Warnings)
              << "*** ERROR *** DenseMatTraits::Create(-,-,true)" << endl
              << "Did not initialize to zero!" << endl;
            return false;
          }
        }
      }

      //Try to change a value.
      DMT::Value(*dm1,0,0) = (ScalarType)5.0; //TODO: Does this compile? Is an lvalue?
      if(DMT::Value(*dm1,0,0) != (ScalarType)5.0){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::Value" << endl
          << "Does not give write access to matrix values!!" << endl;
        return false;
      }

      //Try to sync to host and vice-versa
      DMT::SyncHostToDevice(*dm1);
      DMT::SyncDeviceToHost(*dm1);
      
      //Call create with non-default third arg.
      RCP<DM> dm2 = DMT::Create(numrows, numcols, false);
      
      //Try calling const version? TODO need to check const-ness??
      const ScalarType constval = DMT::Value(*dm1,1,1);
      
      //get stride
      int stride = DMT::GetStride(*dm2);

      //randomize and add
      DMT::Randomize(*dm2);
      ScalarType tmpVal = DMT::Value(*dm1,numrows-1,numcols-1) + DMT::Value(*dm2,numrows-1,numcols-1);
      ScalarType tmpVal2 = DMT::Value(*dm1,0,0) + DMT::Value(*dm2,0,0);
      DMT::Add(*dm1,*dm2);
      if(DMT::Value(*dm1,numrows-1,numcols-1) != tmpVal ||
      DMT::Value(*dm1,0,0) != tmpVal2){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Add failed" << endl;
        return false;
      }

      //Test assign and scale: 
      DMT::Assign(*dm1,*dm2);
      DMT::Scale(*dm1,2.0);
      if(DMT::Value(*dm1,1,1) != (ScalarType)DMT::Value(*dm2,1,1)*2.0){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Assign or scale failed" << endl;
        return false;
      }

      //Test syncs again here? 
      //Really should throw an error because syncing now
      //would overwrite the changes we just made on host.
      //DMT::SyncDeviceToHost(*dm1);
      //
      
      RCP<DM> dm3 = DMT::Subview(*dm1, numrows-1, numcols-1, 1, 1);
      if(DMT::Value(*dm3,0,0) != DMT::Value(*dm1,1,1)){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Subview failed" << endl;
        return false;
      }
      if(DMT::GetNumRows(*dm3) != numrows-1 ||
         DMT::GetNumCols(*dm3) != numcols-1){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Subview gives wrong dimensions." << endl;
        return false;
      }
      // Try to change a value in the subview.
      ScalarType testVal = 237.1;
      DMT::Value(*dm3,0,0) = testVal;
      if(DMT::Value(*dm3,0,0) != DMT::Value(*dm1,1,1) ||
         DMT::Value(*dm3,0,0) != testVal){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Subview did not edit value or did not edit original." << endl;
        return false;
      }

      RCP<DM> dm4 = DMT::SubviewCopy(*dm1, numrows-2, numcols-2, 2, 2);
      if(DMT::Value(*dm4,0,0) != DMT::Value(*dm1,2,2)){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "SubviewCopy failed" << endl;
        return false;
      }
      if(DMT::GetNumRows(*dm4) != numrows-2 ||
         DMT::GetNumCols(*dm4) != numcols-2){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "SubviewCopy gives wrong dimensions." << endl;
        return false;
      }
      // Try to change a value in the subview.
      DMT::Value(*dm4,0,0) = testVal-5;
      if(DMT::Value(*dm4,0,0) == DMT::Value(*dm1,2,2) ||
         DMT::Value(*dm4,0,0) != testVal-5){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "SubviewCopy is incorrect. Possible view but not copy." << endl;
        return false;
      }
      

      //TODO: Try to add matrices of mismatched size (eg dm3, dm4) and make sure it throws error. 
      //
      //Try reshaping:
      //TODO: What should happen if reshape is called on a subview?
      //Increase Dimensions
      DMT::Reshape(*dm2, numrows+5, numcols+5);
      if(DMT::GetNumRows(*dm2) != numrows+5 ||
         DMT::GetNumCols(*dm2) != numcols+5){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Reshape test 1 gives wrong dimensions." << endl;
        return false;
      }
      //Decrease dimensions. 
      DMT::Reshape(*dm2, numrows-2, numcols-2);
      if(DMT::GetNumRows(*dm2) != numrows-2 ||
         DMT::GetNumCols(*dm2) != numcols-2){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Reshape test 2 gives wrong dimensions." << endl;
        return false;
      }

      ScalarType * testPtr = DMT::GetRawHostPtr(*dm2);
      //TODO: Test handing this to lapack function? 
      //TODO: Should something err if view is modified on device and 
      //try to access on host without sync? 


      return true;
    } //end test scope

#if 0
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

#endif
  }

} //namespace Belos

#endif
