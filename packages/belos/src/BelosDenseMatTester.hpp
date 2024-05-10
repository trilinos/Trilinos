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
  template< class ScalarType, class DM > 
  bool
  TestDenseMatTraits (const Teuchos::RCP<OutputManager<ScalarType> > &om)
  {
    using Teuchos::SetScientific;
    using Teuchos::RCP;
    using std::endl;
    typedef DenseMatTraits<ScalarType, DM>    DMT;
    typedef Teuchos::ScalarTraits<ScalarType> STS;

    // Make sure that all floating-point numbers are printed with the
    // right precision.
    SetScientific<ScalarType> sci (om->stream (Warnings));

    // Arbitrary tolerance in case
    // norms are not computed deterministically (which is possible
    // even with MPI only, and more likely with threads).

    /* Dense Traits Contract:
        
        //Example: 
         CloneCopy(MV,vector<int>)
           USER: will request positive number of vectors
             MV: will return a multivector with exactly the number of
                   requested vectors.
                 vectors are the same dimension as the cloned MV


    *********************************************************************/

    //TODO const ScalarType one      = STS::one();
    const ScalarType zero     = STS::zero();

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
     
      DMT::SyncDeviceToHost(*dm1); 
      //Check init to zero on create.
      for(int i = 0; i<numrows; i++){
        for(int j = 0; j<numcols; j++){
          ScalarType val = DMT::ValueConst(*dm1,i,j);
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
      if(DMT::ValueConst(*dm1,0,0) != (ScalarType)5.0){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::Value" << endl
          << "Does not give write access to matrix values!!" << endl;
        return false;
      }

      //Try to sync to host and vice-versa
      DMT::SyncHostToDevice(*dm1);
      
      //Call create with non-default third arg.
      RCP<DM> dm2 = DMT::Create(numrows, numcols, false);
      DMT::PutScalar(*dm2,(ScalarType)47.2);
      DMT::SyncDeviceToHost(*dm2);
      if( DMT::ValueConst(*dm2,0,0) != (ScalarType)47.2){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::PutScalar " << endl
          << "Value at dm2 0,0 is: " << DMT::ValueConst(*dm2,0,0) << endl;
        return false;
      }
      DMT::PutScalar(*dm2);
      DMT::SyncDeviceToHost(*dm2);
      if( DMT::ValueConst(*dm2,0,0) != zero){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "PutScalar default args failed" << endl;
        return false;
      }
      
      //get stride
      int stride = DMT::GetStride(*dm2);
      if (stride != numrows) {
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Stride failed" << endl;
        return false;
      }

      //randomize and add
      DMT::Randomize(*dm2);
      DMT::SyncDeviceToHost(*dm2);
      ScalarType tmpVal = DMT::ValueConst(*dm1,numrows-1,numcols-1) + DMT::ValueConst(*dm2,numrows-1,numcols-1);
      ScalarType tmpVal2 = DMT::ValueConst(*dm1,0,0) + DMT::ValueConst(*dm2,0,0);
      DMT::Add(*dm1,*dm2);
      DMT::SyncDeviceToHost(*dm1);
      if(DMT::ValueConst(*dm1,numrows-1,numcols-1) != tmpVal ||
      DMT::ValueConst(*dm1,0,0) != tmpVal2){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Add failed" << endl;
        return false;
      }

      //Test assign and scale: 
      DMT::Assign(*dm1,*dm2);
      DMT::Scale(*dm1,2.0);
      DMT::SyncDeviceToHost(*dm1);
      if(DMT::ValueConst(*dm1,1,1) != (ScalarType)DMT::ValueConst(*dm2,1,1)*(ScalarType)2.0){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Assign or scale failed" << endl;
        return false;
      }

      RCP<DM> dm3 = DMT::Subview(*dm1, numrows-1, numcols-1, 1, 1);
      DMT::SyncDeviceToHost(*dm3);
      if(DMT::ValueConst(*dm3,0,0) != DMT::ValueConst(*dm1,1,1)){
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
      if(DMT::ValueConst(*dm3,0,0) != DMT::ValueConst(*dm1,1,1) ||
         DMT::ValueConst(*dm3,0,0) != testVal){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Subview did not edit value or did not edit original." << endl;
        return false;
      }

      RCP<DM> dm4 = DMT::SubviewCopy(*dm1, numrows-2, numcols-2, 2, 2);
      if(DMT::ValueConst(*dm4,0,0) != DMT::ValueConst(*dm1,2,2)){
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
      DMT::Value(*dm4,0,0) = testVal-(ScalarType)5;
      if(DMT::ValueConst(*dm4,0,0) == DMT::ValueConst(*dm1,2,2) ||
         DMT::ValueConst(*dm4,0,0) != testVal-(ScalarType)5){
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
      //Increase Dimensions and init to zeros.
      DMT::Reshape(*dm2, numrows+5, numcols+5, true);
      if(DMT::GetNumRows(*dm2) != numrows+5 ||
         DMT::GetNumCols(*dm2) != numcols+5){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Reshape test 3 gives wrong dimensions." << endl;
        return false;
      }
      if(DMT::ValueConst(*dm2,0,0) != zero || DMT::ValueConst(*dm2,numrows+4,numcols+4) != zero){
        om->stream(Warnings)
          << "*** ERROR *** DenseMatTraits::" << endl
          << "Reshape test 3 does not init to zero." << endl;
        return false;
      }

      //ScalarType * testPtr = DMT::GetRawHostPtr(*dm2);
      //DMT::RawPtrDataModified(*dm2);
      //ScalarType const * testPtr2 = DMT::GetConstRawHostPtr(*dm2);
      //TODO: Test handing this to lapack function? 
      //TODO: Should something err if view is modified on device and 
      //try to access on host without sync? 
      //
      //Compute Frobenius norm. //TODO: Do Frob norm of a matrix we know the answer for and check it. 
      DMT::NormFrobenius(*dm2);

      //TODO: Definitely need testing to check host/device sync semantics. 
      return true;
    } //end test scope
  } //end test function
} //namespace Belos

#endif
