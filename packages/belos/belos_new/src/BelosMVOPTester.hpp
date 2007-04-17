// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
#ifndef _BELOS_MVOPTESTER_HPP
#define _BELOS_MVOPTESTER_HPP

// Assumptions that I have made:
// * I assume/verify that a multivector must have at least one vector. This seems 
//   to be consistent with Epetra_MultiVec.
// * I do not assume that an operator is deterministic; I do assume that the
//   operator, applied to 0, will return 0.

/*! \file BelosMVOPTester.hpp
  \brief Test routines for MultiVecTraits and OperatorTraits conformity.
*/

#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosTypes.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "BelosOutputManager.hpp"

namespace Belos {

/*! \fn TestMultiVecTraits
    \brief This is a function to test the correctness of a MultiVecTraits 
    specialization and multivector implementation.
    \return Status of the test: OK or Error
*/
  template< class ScalarType, class MV >
  ReturnType TestMultiVecTraits( 
                const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                const Teuchos::RefCountPtr<const MV> &A ) {

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

         GetVecLength 
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

         MvPrint
             MV: routine does not modify vectors (not tested here)
    *********************************************************************/

    typedef MultiVecTraits<ScalarType, MV>    MVT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType       MagType;
    typedef Teuchos::ScalarTraits<MagType>    MT;
    
    const ScalarType one      = SCT::one();
    const ScalarType zero     = SCT::zero();
    const MagType    zero_mag = Teuchos::ScalarTraits<MagType>::zero();

    // Don't change these two without checking the initialization of ind below
    const int numvecs   = 10;
    const int numvecs_2 = 5;

    int i,j;
    std::vector<int> ind(numvecs_2);
    std::ostream &out = *om->GetOStream();

    /* Initialize indices for selected copies/views
       The MVT specialization should not assume that 
       these are ordered or even distinct.
       Also retrieve the edges.

       However, to spice things up, grab the first vector,
       last vector, and choose the others randomly.
    */
    ind[0] = 0;
    ind[numvecs_2-1] = numvecs-1;
    for (i=1; i<numvecs_2-1; i++) {
      ind[i] = rand() % numvecs;
    }


    /*********** GetNumberVecs() *****************************************
       Verify:
       1) This number should be strictly positive
    *********************************************************************/
    if ( MVT::GetNumberVecs(*A) <= 0 ) {
      if ( om->isVerbosityAndPrint( Errors  ) ) {
        out << "*** ERROR *** MultiVectorTraits::GetNumberVecs()." << endl
            << "Returned <= 0." << endl;
      }
      return Error;
    }


    /*********** GetVecLength() ******************************************
       Verify:
       1) This number should be strictly positive
    *********************************************************************/
    if ( MVT::GetVecLength(*A) <= 0 ) {
      if ( om->isVerbosityAndPrint( Errors ) ) {
        out << "*** ERROR *** MultiVectorTraits::GetVecLength()" << endl
            << "Returned <= 0." << endl;
      }
      return Error;
    }


    /*********** Clone() and MvNorm() ************************************
       Verify:
       1) Clone() allows us to specify the number of vectors
       2) Clone() returns a multivector of the same dimension
       3) Vector norms shouldn't be negative
       4) MvNorm result vector should not be resized
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(2*numvecs);
      if ( MVT::GetNumberVecs(*B) != numvecs ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::Clone()." << endl
              << "Did not allocate requested number of vectors." << endl;
        }
        return Error;
      }
      if ( MVT::GetVecLength(*B) != MVT::GetVecLength(*A) ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::Clone()." << endl
              << "Did not allocate requested number of vectors." << endl;
        }
        return Error;
      }
      MVT::MvNorm(*B, &norms);
      if ( norms.size() != 2*numvecs ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::MvNorm()." << endl
              << "Method resized the output vector." << endl;
        }
        return Error;
      }
      for (i=0; i<numvecs; i++) {
        if ( norms[i] < zero_mag ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::Clone()." << endl
                << "Vector had negative norm." << endl;
          }
          return Error;
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
       6) Verify that vector norms are different than before
       
       Without knowing something about the random distribution, 
       this is about the best that we can do, to make sure that MvRandom 
       did at least *something*.
       
       Also, make sure vector norms aren't negative.
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(numvecs), norms2(numvecs);

      MVT::MvInit(*B);
      MVT::MvNorm(*B, &norms);
      for (i=0; i<numvecs; i++) {
        if ( norms[i] != zero_mag ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvInit() "
                << "and MultiVecTraits::MvNorm()" << endl
                << "Supposedly zero vector has non-zero norm." << endl;
          }
          return Error;
        }
      }
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, &norms);
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, &norms2);
      for (i=0; i<numvecs; i++) {
        if ( norms[i] == zero_mag || norms2[i] == zero_mag ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvRandom()." << endl
                << "Random vector was empty (very unlikely)." << endl;
          }
          return Error;
        }
        else if ( norms[i] < zero_mag || norms2[i] < zero_mag ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvRandom()." << endl
                << "Vector had negative norm." << endl;
          }
          return Error;
        }
        else if ( norms[i] == norms2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MutliVecTraits::MvRandom()." << endl
                << "Vectors not random enough." << endl;
          }
          return Error;
        }
      }
    }


    /*********** MvInit() and MvNorm() ***********************************
       A vector of ones of dimension n should have norm sqrt(n)
       1) Init vectors to all ones
       2) Verify that norm is sqrt(n)
       3) Verify that norms aren't negative

       Note: I'm not sure that we can expect this to hold in practice.
              Maybe something like abs(norm-sqrt(n)) < SCT::eps()  ???
              The sum of 1^2==1 should be n, but what about sqrt(n)?
              They may be using a different square root than ScalartTraits
              On my iBook G4 and on jeter, this test works.
              Right now, this has been demoted to a warning.
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(numvecs);

      MVT::MvInit(*B,one);
      MVT::MvNorm(*B, &norms);
      bool BadNormWarning = false;
      for (i=0; i<numvecs; i++) {
        if ( norms[i] < zero_mag ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvRandom()." << endl
                << "Vector had negative norm." << endl;
          }
          return Error;
        }
        else if ( norms[i] != MT::squareroot(MVT::GetVecLength(*B)) && !BadNormWarning ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << endl;
            out << "Warning testing MultiVecTraits::MvInit()." << endl
                << "Ones vector should have norm sqrt(dim)." << endl
                << "norms[i]: " << norms[i] << "\tdim: " << MVT::GetVecLength(*B) << endl << endl;
          }
          BadNormWarning = true;
        }
      }
    }


    /*********** MvInit() and MvNorm() ***********************************
       A vector of zeros of dimension n should have norm 0
       1) Verify that norms aren't negative
       2) Verify that norms are zero

       We must know this works before the next tests.
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B = MVT::Clone(*A,numvecs);
      std::vector<MagType> norms(numvecs);
      MVT::MvInit(*B, zero_mag);
      MVT::MvNorm(*B, &norms);
      for (i=0; i<numvecs; i++) {
        if ( norms[i] < zero_mag ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvInit()." << endl
                << "Vector had negative norm." << endl;
          }
          return Error;
        }
        else if ( norms[i] != zero_mag ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvInit()." << endl
                << "Zero vector should have norm zero." << endl;
          }
          return Error;
        }
      }
    }


    /*********** CloneCopy(MV,vector<int>) and MvNorm ********************
       1) Check quantity/length of vectors
       2) Check vector norms for agreement
       3) Zero out B and make sure that C norms are not affected
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B, C;
      std::vector<MagType> norms(numvecs), norms2(numvecs);

      B = MVT::Clone(*A,numvecs);
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, &norms);
      C = MVT::CloneCopy(*B,ind);
      MVT::MvNorm(*C, &norms2);
      if ( MVT::GetNumberVecs(*C) != numvecs_2 ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::CloneCopy(ind)." << endl
              << "Wrong number of vectors." << endl;
        }
        return Error;
      }
      if ( MVT::GetVecLength(*C) != MVT::GetVecLength(*B) ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::CloneCopy(ind)." << endl
              << "Vector lengths don't match." << endl;
        }
        return Error;
      }
      for (i=0; i<numvecs_2; i++) {
        if ( norms2[i] != norms[ind[i]] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::CloneCopy(ind)." << endl
                << "Copied vectors do not agree:" 
                << norms2[i] << " != " << norms[ind[i]] << endl;
          }
          return Error;
        }
      }
      MVT::MvInit(*B,zero);
      MVT::MvNorm(*C, &norms); 
      for (i=0; i<numvecs_2; i++) {
        if ( norms2[i] != norms[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::CloneCopy(ind)." << endl
                << "Copied vectors were not independent." << endl;
          }
          return Error;
        }
      }
    }    


    /*********** CloneCopy(MV) and MvNorm ********************************
       1) Check quantity
       2) Check value of norms
       3) Zero out B and make sure that C is still okay
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B, C;
      std::vector<MagType> norms(numvecs), norms2(numvecs);

      B = MVT::Clone(*A,numvecs);
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, &norms);
      C = MVT::CloneCopy(*B);
      MVT::MvNorm(*C, &norms2);
      if ( MVT::GetNumberVecs(*C) != numvecs ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::CloneCopy()." << endl
              << "Wrong number of vectors." << endl;
        }
        return Error;
      }
      for (i=0; i<numvecs; i++) {
        if ( norms2[i] != norms[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::CloneCopy()." << endl
                << "Copied vectors do not agree." << endl;
          }
          return Error;
        }
      }
      MVT::MvInit(*B,zero);
      MVT::MvNorm(*C, &norms); 
      for (i=0; i<numvecs; i++) {
        if ( norms2[i] != norms[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::CloneCopy()." << endl
                << "Copied vectors were not independent." << endl;
          }
          return Error;
        }
      }
    }


    /*********** CloneView(MV,vector<int>) and MvNorm ********************
       Check that we have a view of the selected vectors
       1) Check quantity
       2) Check value of norms
       3) Zero out B and make sure that C is zero as well
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B, C;
      std::vector<MagType> norms(numvecs), norms2(numvecs);

      B = MVT::Clone(*A,numvecs); 
      MVT::MvRandom(*B);
      MVT::MvNorm(*B, &norms);
      C = MVT::CloneView(*B,ind);
      MVT::MvNorm(*C, &norms2);
      if ( MVT::GetNumberVecs(*C) != numvecs_2 ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::CloneView(ind)." << endl
              << "Wrong number of vectors." << endl;
        }
        return Error;
      }
      for (i=0; i<numvecs_2; i++) {
        if ( norms2[i] != norms[ind[i]] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::CloneView(ind)." << endl
                << "Viewed vectors do not agree." << endl;
          }
          return Error;
        }
      }
      /*
      MVT::MvInit(*B,zero);
      MVT::MvNorm(*C, &norms2); 
      for (i=0; i<numvecs_2; i++) {
        if ( norms2[i] != zero ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::CloneView(ind)." << endl
                << "Copied vectors were not dependent." << endl;
          }
          return Error;
        }
      }
      */
    }


    /*********** const CloneView(MV,vector<int>) and MvNorm() ************
       Check that we have a view of the selected vectors.
       1) Check quantity
       2) Check value of norms for agreement
       3) Zero out B and make sure that C is zerod as well
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B;
      Teuchos::RefCountPtr<const MV> constB, C;
      std::vector<MagType> normsB(numvecs), normsC(numvecs_2);
      std::vector<int> allind(numvecs);
      for (i=0; i<numvecs; i++) {
        allind[i] = i;
      }

      B = MVT::Clone(*A,numvecs);
      MVT::MvRandom( *B );
      // need a const MV to test const CloneView
      constB = MVT::CloneView(*B,allind);
      MVT::MvNorm(*constB, &normsB);
      C = MVT::CloneView(*constB,ind);
      MVT::MvNorm(*C, &normsC);
      if ( MVT::GetNumberVecs(*C) != numvecs_2 ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** const MultiVecTraits::CloneView(ind)." << endl
              << "Wrong number of vectors." << endl;
        }
        return Error;
      }
      for (i=0; i<numvecs_2; i++) {
        if ( normsC[i] != normsB[ind[i]] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** const MultiVecTraits::CloneView(ind)." << endl
                << "Viewed vectors do not agree." << endl;
          }
          return Error;
        }
      }
      /*
      MVT::MvInit(const_cast<MV&>(*C),zero);
      MVT::MvNorm(*constB, &normsB); 
      for (i=0; i<numvecs_2; i++) {
        if ( normsB[ind[i]] != SCT::zero() ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** const MultiVecTraits::CloneView(ind)." << endl
                << "Copied vectors were not dependent." << endl;
          }
          return Error;
        }
      }
      */
    }


    /*********** SetBlock() and MvNorm() *********************************
       SetBlock() will copy the vectors from C into B 
       1) Verify that the specified vectors were copied
       2) Verify that the other vectors were not modified
       3) Verify that C was not modified
       4) Change C and then check B to make sure it was not modified
      
       Use a different index set than has been used so far (distinct entries).
       This is because duplicate entries will cause the vector to be
       overwritten, making it more difficult to test.
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B, C;
      std::vector<MagType> normsB1(numvecs), normsB2(numvecs),
                           normsC1(numvecs_2), normsC2(numvecs_2);

      B = MVT::Clone(*A,numvecs);
      C = MVT::Clone(*A,numvecs_2);
      // Just do every other one, interleaving the vectors of C into B
      ind.resize(numvecs_2);
      for (i=0; i<numvecs_2; i++) {
        ind[i] = 2*i;
      }
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);

      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      MVT::SetBlock(*C,ind,*B);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);

      // check that C was not changed by SetBlock
      for (i=0; i<numvecs_2; i++) {
        if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
                << "Operation modified source vectors." << endl;
          }
          return Error;
        }
      }
      // check that the correct vectors of B were modified
      // and the others were not
      for (i=0; i<numvecs; i++) {
        if (i % 2 == 0) {
          // should be a vector from C
          if ( normsB2[i] != normsC1[i/2] ) {
            if ( om->isVerbosityAndPrint( Errors ) ) {
              out << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
                  << "Copied vectors do not agree." << endl;
            }
            return Error;
          }
        }
        else {
          // should be an original vector
          if ( normsB1[i] != normsB2[i] ) {
            if ( om->isVerbosityAndPrint( Errors ) ) {
              out << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
                  << "Incorrect vectors were modified." << endl;
            }
            return Error;
          }
        }
      }
      MVT::MvInit(*C,zero);
      MVT::MvNorm(*B,&normsB1);
      // verify that we copied and didn't reference
      for (i=0; i<numvecs; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
                << "Copied vectors were not independent." << endl;
          }
          return Error;
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
       This is because duplicate entries will cause the vector to be
       overwritten, making it more difficult to test.

       These tests are the same as the ones above, except that the
       number of indices (to be copied into B) is less than the number
       of vectors in C, so that not all of C is put into B.
    *********************************************************************/
    {
      Teuchos::RefCountPtr<MV> B, C;
      // set these: we assume below that setSize*2=BSize
      const int BSize   = 10, 
                CSize   = 6,
                setSize = 5;
      std::vector<MagType> normsB1(BSize), normsB2(BSize),
                           normsC1(CSize), normsC2(CSize);

      B = MVT::Clone(*A,BSize);
      C = MVT::Clone(*A,CSize);
      // Just do every other one, interleaving the vectors of C into B
      ind.resize(setSize);
      for (i=0; i<setSize; i++) {
        ind[i] = 2*i;
      }
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);

      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      MVT::SetBlock(*C,ind,*B);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);

      // check that C was not changed by SetBlock
      for (i=0; i<CSize; i++) {
        if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
                << "Operation modified source vectors." << endl;
          }
          return Error;
        }
      }
      // check that the correct vectors of B were modified
      // and the others were not
      for (i=0; i<BSize; i++) {
        if (i % 2 == 0) {
          // should be a vector from C
          if ( normsB2[i] != normsC1[i/2] ) {
            if ( om->isVerbosityAndPrint( Errors ) ) {
              out << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
                  << "Copied vectors do not agree." << endl;
            }
            return Error;
          }
        }
        else {
          // should be an original vector
          if ( normsB1[i] != normsB2[i] ) {
            if ( om->isVerbosityAndPrint( Errors ) ) {
              out << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
                  << "Incorrect vectors were modified." << endl;
            }
            return Error;
          }
        }
      }
      MVT::MvInit(*C,zero);
      MVT::MvNorm(*B,&normsB1);
      // verify that we copied and didn't reference
      for (i=0; i<numvecs; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::SetBlock()." << endl
                << "Copied vectors were not independent." << endl;
          }
          return Error;
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
      Teuchos::RefCountPtr<MV> B, C;
      std::vector<MagType> normsB(p), normsC(q);
      Teuchos::SerialDenseMatrix<int,ScalarType> SDM(p,q);

      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,q);
   
      // randomize the multivectors
      MVT::MvRandom(*B);
      MVT::MvNorm(*B,&normsB);
      MVT::MvRandom(*C);
      MVT::MvNorm(*C,&normsC);
   
      // perform SDM  = zero() * B^H * C
      MVT::MvTransMv( zero, *B, *C, SDM );
   
      // check the sizes: not allowed to have shrunk
      if ( SDM.numRows() != p || SDM.numCols() != q ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
              << "Routine resized SerialDenseMatrix." << endl;
        }
        return Error;
      }
   
      // check that zero**A^H*B == zero
      if ( SDM.normOne() != zero ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
              << "Scalar argument processed incorrectly." << endl;
        }
        return Error;
      }
   
      // perform SDM  = one * B^H * C
      MVT::MvTransMv( one, *B, *C, SDM );
   
      // check the norms: a^H b = |a| |b| cos(theta) <= |a| |b|
      // with equality only when a and b are colinear
      for (i=0; i<p; i++) {
        for (j=0; j<q; j++) {
          if (   SCT::magnitude(SDM(i,j)) 
               > SCT::magnitude(normsB[i]*normsC[j]) ) {
            if ( om->isVerbosityAndPrint( Errors ) ) {
              out << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
                  << "Triangle inequality did not hold: " 
                  << SCT::magnitude(SDM(i,j)) 
                  << " > " 
                  << SCT::magnitude(normsB[i]*normsC[j]) 
                  << endl;
            }
            return Error;
          }
        }
      }
      MVT::MvInit(*C);
      MVT::MvRandom(*B);
      MVT::MvTransMv( one, *B, *C, SDM );
      for (i=0; i<p; i++) {
        for (j=0; j<q; j++) {
          if ( SDM(i,j) != zero ) {
            if ( om->isVerbosityAndPrint( Errors ) ) {
              out << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
                  << "Inner products not zero for C==0." << endl;
            }
            return Error;
          }
        }
      }
      MVT::MvInit(*B);
      MVT::MvRandom(*C);
      MVT::MvTransMv( one, *B, *C, SDM );
      for (i=0; i<p; i++) {
        for (j=0; j<q; j++) {
          if ( SDM(i,j) != zero ) {
            if ( om->isVerbosityAndPrint( Errors ) ) {
              out << "*** ERROR *** MultiVecTraits::MvTransMv()." << endl
                  << "Inner products not zero for B==0." << endl;
            }
            return Error;
          }
        }
      }
    } 


    /*********** MvDot() *************************************************
        Verify:
        1) Results vector not resized
        2) Inner product inequalities are satisfied
        3) Zero vectors give zero inner products
    *********************************************************************/
    {
      const int p = 7;
      const int q = 9;
      Teuchos::RefCountPtr<MV> B, C;
      vector<ScalarType> iprods(p+q);
      std::vector<MagType> normsB(numvecs), normsC(numvecs);

      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,p);

      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB);
      MVT::MvNorm(*C,&normsC);
      MVT::MvDot( *B, *C, &iprods );
      if ( iprods.size() != p+q ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** MultiVecTraits::MvDot." << endl
              << "Routine resized results vector." << endl;
        }
        return Error;
      }
      for (i=0; i<BELOS_MIN(p,q); i++) {
        if ( SCT::magnitude(iprods[i]) 
             > SCT::magnitude(normsB[i]*normsC[i]) ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvDot()." << endl
                << "Inner products not valid." << endl;
          }
          return Error;
        }
      }
      MVT::MvInit(*B);
      MVT::MvRandom(*C);
      MVT::MvDot( *B, *C, &iprods );
      for (i=0; i<p; i++) {
        if ( iprods[i] != zero ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvDot()." << endl
                << "Inner products not zero for B==0." << endl;
          }
          return Error;
        }
      }
      MVT::MvInit(*C);
      MVT::MvRandom(*B);
      MVT::MvDot( *B, *C, &iprods );
      for (i=0; i<p; i++) {
        if ( iprods[i] != zero ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvDot()." << endl
                << "Inner products not zero for C==0." << endl;
          }
          return Error;
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
      Teuchos::RefCountPtr<MV> B, C, D;
      std::vector<MagType> normsB1(p), normsB2(p),
                           normsC1(p), normsC2(p),
                           normsD1(p), normsD2(p);
      ScalarType alpha = SCT::random(),
                  beta = SCT::random();

      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,p);
      D = MVT::Clone(*A,p);

      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
   
      // check that 0*B+1*C == C
      MVT::MvAddMv(zero,*B,one,*C,*D);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      MVT::MvNorm(*D,&normsD1);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Input arguments were modified." << endl;
          }
          return Error;
        }
        else if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Input arguments were modified." << endl;
          }
          return Error;
        }
        else if ( normsC1[i] != normsD1[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Assignment did not work." << endl;
          }
          return Error;
        }
      }

      // check that 1*B+0*C == B
      MVT::MvAddMv(one,*B,zero,*C,*D);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      MVT::MvNorm(*D,&normsD1);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Input arguments were modified." << endl;
          }
          return Error;
        }
        else if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Input arguments were modified." << endl;
          }
          return Error;
        }
        else if ( normsB1[i] != normsD1[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Assignment did not work." << endl;
          }
          return Error;
        }
      }

      // check that alpha*B+beta*C -> D is invariant under initial D
      // first, try random D
      MVT::MvRandom(*D);
      MVT::MvAddMv(alpha,*B,beta,*C,*D);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      MVT::MvNorm(*D,&normsD1);
      // check that input args are not modified
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Input arguments were modified." << endl;
          }
          return Error;
        }
        else if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Input arguments were modified." << endl;
          }
          return Error;
        }
      }
      // next, try zero D
      MVT::MvInit(*D);
      MVT::MvAddMv(alpha,*B,beta,*C,*D);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      MVT::MvNorm(*D,&normsD2);
      // check that input args are not modified and that D is the same
      // as the above test
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Input arguments were modified." << endl;
          }
          return Error;
        }
        else if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Input arguments were modified." << endl;
          }
          return Error;
        }
        else if ( normsD1[i] != normsD2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv()." << endl
                << "Results varies depending on initial state of dest vectors." << endl;
          }
          return Error;
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
      Teuchos::RefCountPtr<MV> B, C, D;
      std::vector<MagType> normsB(p),
                           normsD(p);
      std::vector<int> lclindex(p);
      for (i=0; i<p; i++) lclindex[i] = i;

      B = MVT::Clone(*A,p);
      C = MVT::CloneView(*B,lclindex);
      D = MVT::CloneView(*B,lclindex);

      MVT::MvRandom(*B);
      MVT::MvNorm(*B,&normsB);
   
      // check that 0*B+1*C == C
      MVT::MvAddMv(zero,*B,one,*C,*D);
      MVT::MvNorm(*D,&normsD);
      for (i=0; i<p; i++) {
        if ( normsB[i] != normsD[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv() #2" << endl
                << "Assignment did not work." << endl;
          }
          return Error;
        }
      }

      // check that 1*B+0*C == B
      MVT::MvAddMv(one,*B,zero,*C,*D);
      MVT::MvNorm(*D,&normsD);
      for (i=0; i<p; i++) {
        if ( normsB[i] != normsD[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvAddMv() #2" << endl
                << "Assignment did not work." << endl;
          }
          return Error;
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
      Teuchos::RefCountPtr<MV> B, C;
      Teuchos::SerialDenseMatrix<int,ScalarType> SDM(p,q);
      std::vector<MagType> normsC1(q), normsC2(q),
                           normsB1(p), normsB2(p);
      
      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,q);

      // Test 1: alpha==0, SDM!=0, beta==1 and check that C is unchanged
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      SDM.random();
      MVT::MvTimesMatAddMv(zero,*B,SDM,one,*C);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Input vectors were modified." << endl;
          }
          return Error;
        }
      }
      for (i=0; i<q; i++) {
        if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 1 failed." << endl;
          }
          return Error;
        }
      }

      // Test 2: alpha==0, SDM!=0, beta==0 and check that C is set to zero
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      SDM.random();
      MVT::MvTimesMatAddMv(zero,*B,SDM,zero,*C);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Input vectors were modified." << endl;
          }
          return Error;
        }
      }
      for (i=0; i<q; i++) {
        if ( normsC2[i] != zero ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 2 failed: " 
                << normsC2[i] 
                << " != " 
                << zero 
                << endl;
          }
          return Error;
        }
      }

      // Test 3: alpha==1, SDM==|I|, beta==0 and check that C is set to B
      //                        |0|
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      SDM.scale(zero);
      for (i=0; i<q; i++) {
        SDM(i,i) = one;
      }
      MVT::MvTimesMatAddMv(one,*B,SDM,zero,*C);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Input vectors were modified." << endl;
          }
          return Error;
        }
      }
      for (i=0; i<q; i++) {
        if ( normsB1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 3 failed: "
                << normsB1[i] 
                << " != "
                << normsC2[i]
                << endl;
          }
          return Error;
        }
      }

      // Test 4: alpha==1, SDM==0, beta==1 and check that C is unchanged
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      SDM.scale(zero);
      MVT::MvTimesMatAddMv(one,*B,SDM,one,*C);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Input vectors were modified." << endl;
          }
          return Error;
        }
      }
      for (i=0; i<q; i++) {
        if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 4 failed." << endl;
          }
          return Error;
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
      Teuchos::RefCountPtr<MV> B, C;
      Teuchos::SerialDenseMatrix<int,ScalarType> SDM(p,q);
      std::vector<MagType> normsC1(q), normsC2(q),
                           normsB1(p), normsB2(p);
      
      B = MVT::Clone(*A,p);
      C = MVT::Clone(*A,q);

      // Test 5: alpha==0, SDM!=0, beta==1 and check that C is unchanged
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      SDM.random();
      MVT::MvTimesMatAddMv(zero,*B,SDM,one,*C);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Input vectors were modified." << endl;
          }
          return Error;
        }
      }
      for (i=0; i<q; i++) {
        if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 5 failed." << endl;
          }
          return Error;
        }
      }

      // Test 6: alpha==0, SDM!=0, beta==0 and check that C is set to zero
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      SDM.random();
      MVT::MvTimesMatAddMv(zero,*B,SDM,zero,*C);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Input vectors were modified." << endl;
          }
          return Error;
        }
      }
      for (i=0; i<q; i++) {
        if ( normsC2[i] != zero ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 6 failed: " 
                << normsC2[i] 
                << " != " 
                << zero 
                << endl;
          }
          return Error;
        }
      }

      // Test 7: alpha==1, SDM==[I 0], beta==0 and check that C is set to B
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      SDM.scale(zero);
      for (i=0; i<p; i++) {
        SDM(i,i) = one;
      }
      MVT::MvTimesMatAddMv(one,*B,SDM,zero,*C);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Input vectors were modified." << endl;
          }
          return Error;
        }
      }
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 7 failed." << endl;
          }
          return Error;
        }
      }
      for (i=p; i<q; i++) {
        if ( normsC2[i] != zero ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 7 failed." << endl;
          }
          return Error;
        }
      }

      // Test 8: alpha==1, SDM==0, beta==1 and check that C is unchanged
      MVT::MvRandom(*B);
      MVT::MvRandom(*C);
      MVT::MvNorm(*B,&normsB1);
      MVT::MvNorm(*C,&normsC1);
      SDM.scale(zero);
      MVT::MvTimesMatAddMv(one,*B,SDM,one,*C);
      MVT::MvNorm(*B,&normsB2);
      MVT::MvNorm(*C,&normsC2);
      for (i=0; i<p; i++) {
        if ( normsB1[i] != normsB2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Input vectors were modified." << endl;
          }
          return Error;
        }
      }
      for (i=0; i<q; i++) {
        if ( normsC1[i] != normsC2[i] ) {
          if ( om->isVerbosityAndPrint( Errors ) ) {
            out << "*** ERROR *** MultiVecTraits::MvTimesMatAddMv()." << endl
                << "Arithmetic test 8 failed." << endl;
          }
          return Error;
        }
      }
    }

    return Ok;

  }



/*! \fn TestOperatorTraits
    \brief This function tests the correctness of an operator 
    implementation with respect to an OperatorTraits specialization 
    \return Status of the test: OK or Error
*/
  template< class ScalarType, class MV, class OP>
  ReturnType TestOperatorTraits( 
                const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                const Teuchos::RefCountPtr<const MV> &A,
                const Teuchos::RefCountPtr<const OP> &M) {

    /* OPT Contract:
       Apply()
         MV: OP*zero == zero
             Warn if OP is not deterministic (OP*A != OP*A)
             Does not modify input arguments
    *********************************************************************/

    typedef MultiVecTraits<ScalarType, MV>     MVT;
    typedef Teuchos::ScalarTraits<ScalarType>  SCT;
    typedef OperatorTraits<ScalarType, MV, OP> OPT;
    typedef typename SCT::magnitudeType        MagType;

    const int numvecs = 10;

    Teuchos::RefCountPtr<MV> B = MVT::Clone(*A,numvecs), 
                             C = MVT::Clone(*A,numvecs);

    std::vector<MagType> normsB1(numvecs), normsB2(numvecs),
                         normsC1(numvecs), normsC2(numvecs);
    ReturnType ret;
    std::ostream &out = *om->GetOStream();
    bool NonDeterministicWarning;
    int i;


    /*********** Apply() *************************************************
        Verify:
        1) OP*B == OP*B; OP is deterministic (just warn on this)
        2) OP*zero == 0
        3) OP*B doesn't modify B
        4) OP*B is invariant under initial state of destination vectors
    *********************************************************************/
    MVT::MvInit(*B);
    MVT::MvRandom(*C);
    MVT::MvNorm(*B,&normsB1);
    ret = OPT::Apply(*M,*B,*C);
    MVT::MvNorm(*B,&normsB2);
    MVT::MvNorm(*C,&normsC2);
    if (ret != Ok) {
      if (om->isVerbosityAndPrint( Errors )) {
        out << "*** ERROR *** OperatorTraits::Apply() [1]" << endl
            << "Apply() returned an error." << endl;
      }
      return Error;
    }
    for (i=0; i<numvecs; i++) {
      if (normsB2[i] != normsB1[i]) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** OperatorTraits::Apply() [1]" << endl
              << "Apply() modified the input vectors." << endl;
        }
        return Error;
      }
      if (normsC2[i] != SCT::zero()) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** OperatorTraits::Apply() [1]" << endl
              << "Operator applied to zero did not return zero." << endl;
        }
        return Error;
      }
    }

    // If we send in a random matrix, we should not get a zero return
    MVT::MvRandom(*B);
    MVT::MvNorm(*B,&normsB1);
    ret = OPT::Apply(*M,*B,*C);
    MVT::MvNorm(*B,&normsB2);
    MVT::MvNorm(*C,&normsC2);
    if (ret != Ok) {
      if (om->isVerbosityAndPrint( Errors )) {
        out << "*** ERROR *** OperatorTraits::Apply() [2]" << endl
            << "Apply() returned an error." << endl;
      }
      return Error;
    }
    bool ZeroWarning = false;
    for (i=0; i<numvecs; i++) {
      if (normsB2[i] != normsB1[i]) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** OperatorTraits::Apply() [2]" << endl
              << "Apply() modified the input vectors." << endl;
        }
        return Error;
      }
      if (normsC2[i] == SCT::zero() && ZeroWarning==false ) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** OperatorTraits::Apply() [2]" << endl
              << "Operator applied to random vectors returned zero." << endl;
          ZeroWarning = true;
        }
      }
    }

    // Apply operator with C init'd to zero
    MVT::MvRandom(*B);
    MVT::MvNorm(*B,&normsB1);
    MVT::MvInit(*C);
    ret = OPT::Apply(*M,*B,*C);
    MVT::MvNorm(*B,&normsB2);
    MVT::MvNorm(*C,&normsC1);
    if (ret != Ok) {
      if ( om->isVerbosityAndPrint( Errors ) ) {
        out << "*** ERROR *** OperatorTraits::Apply() [3]" << endl
            << "Apply() returned an error." << endl;
      }
      return Error;
    }
    for (i=0; i<numvecs; i++) {
      if (normsB2[i] != normsB1[i]) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** OperatorTraits::Apply() [3]" << endl
              << "Apply() modified the input vectors." << endl;
        }
        return Error;
      }
    }

    // Apply operator with C init'd to random
    // Check that result is the same as before; warn if not.
    // This could be a result of a bug, or a stochastic
    //   operator. We do not want to prejudice against a 
    //   stochastic operator.
    MVT::MvRandom(*C);
    ret = OPT::Apply(*M,*B,*C);
    MVT::MvNorm(*B,&normsB2);
    MVT::MvNorm(*C,&normsC2);
    NonDeterministicWarning = false;
    if (ret != Ok) {
      if (om->isVerbosityAndPrint( Errors )) {
        out << "*** ERROR *** OperatorTraits::Apply() [4]" << endl
            << "Apply() returned an error." << endl;
      }
      return Error;
    }
    for (i=0; i<numvecs; i++) {
      if (normsB2[i] != normsB1[i]) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << "*** ERROR *** OperatorTraits::Apply() [4]" << endl
              << "Apply() modified the input vectors." << endl;
        }
        return Error;
      }
      if (normsC1[i] != normsC2[i] && !NonDeterministicWarning) {
        if ( om->isVerbosityAndPrint( Errors ) ) {
          out << endl;
          out << "*** WARNING *** OperatorTraits::Apply() [4]" << endl
              << "Apply() returned two different results." << endl << endl;
        }
        NonDeterministicWarning = true;
      }
    }

    return Ok;

  }

}

#endif
