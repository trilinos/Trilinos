// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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


/*! \file AnasaziBasicOrthoManager.hpp
  \brief Basic implementation of the Anasazi::OrthoManager class
*/

#ifndef ANASAZI_BASIC_ORTHOMANAGER_HPP
#define ANASAZI_BASIC_ORTHOMANAGER_HPP

/*!   \class Anasazi::BasicOrthoManager
      \brief An implementation of the Anasazi::MatOrthoManager that performs orthogonalization
      using (potentially) multiple steps of classical Gram-Schmidt.
      
      \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/


// #define ORTHO_DEBUG

#include "AnasaziMatOrthoManager.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"

namespace Anasazi {

  template<class ScalarType, class MV, class OP>
  class BasicOrthoManager : public MatOrthoManager<ScalarType,MV,OP> {

  private:
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<ScalarType>  SCT;
    typedef MultiVecTraits<ScalarType,MV>      MVT;
    typedef OperatorTraits<ScalarType,MV,OP>   OPT;

  public:
    
    //@{ \name Constructor/Destructor.
    //! Constructor specifying re-orthogonalization tolerance.
    BasicOrthoManager( Teuchos::RefCountPtr<const OP> Op = Teuchos::null,
                       const MagnitudeType kappa = SCT::magnitude(1.5625) ) : MatOrthoManager<ScalarType,MV,OP>(Op), _kappa(kappa) {};

    //! Destructor
    ~BasicOrthoManager() {};
    //@}


    //@{ \name Accessor routines

    //! Set parameter for re-orthogonalization threshhold.
    void setKappa( const MagnitudeType kappa ) { _kappa = kappa; };

    //! Return parameter for re-orthogonalization threshhold.
    MagnitudeType getKappa() const { return _kappa; } 

    //@} 


    //@{ \name Orthogonalization methods.

    /*! \brief This method takes a multivector and projects it onto the space orthogonal to 
     *  another given multivector, in a specified inner product. The method has the option of
     *  exploiting a caller-provided \c MX, and returning updated information to the caller.
     *  
     @param X [in/out] The multivector to the modified. 
       On output, \c X will be orthogonal to \c Q with respect to \c innerProd().
      
     @param MX [in/out] The image of \c X under the operator \c Op. 
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.
            
     @param Q [in] A multivector specifying the space to be orthogonalized against. 
            
     @return Code specifying failure of the routine, as defined by the implementation.
    */
    ReturnType project ( MV &X, Teuchos::RefCountPtr<MV> MX, const MV &Q ) const;


    /*! \brief This method takes a multivector and orthonormalizes the columns, with respect to \c innerProd().
     *  The method has the option of
     *  exploiting a caller-provided \c MX, and returning updated information to the caller.
     *  
     @param X [in/out] The multivector to the modified. 
       On output, the columns are Op-orthonormal
    
     @param MX [in/out] The image of \c X under the operator \c Op. 
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.
      
     @param rank [out] Rank of the basis computed by this method.
    
     @return Code specifying failure of the routine, as defined by the implementation.
    */
    ReturnType normalize ( MV &X, Teuchos::RefCountPtr<MV> MX, int &rank ) const;


    /*! \brief This method takes a multivector and projects it onto the space orthogonal to 
     *  another given multivector.  It also orthonormalizes the 
     *  columns of the resulting multivector. Both of these operations are conducted 
     *  with respect to \c innerProd().
     *  The method has the option of
     *  exploiting a caller-provided \c MX, and returning updated information to the caller.
     *  
     @param X [in/out] The multivector to the modified. 
       On output, the columns of X are Op-orthonormal and Op-orthogonal to Q.
      
     @param MX [in/out] The image of \c X under the operator \c Op. 
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.
      
     @param Q [in] A multivector specifying the space to be orthogonalized against. \c Q is assumed to have orthonormal
     columns with respect to \c innerProd().
      
     @param rank [out] Rank of the basis computed by this method.
      
     @return Code specifying failure of the routine, as defined by the implementation.
    */
    ReturnType projectAndNormalize ( MV &X, Teuchos::RefCountPtr<MV> MX, const MV &Q, int &rank ) const;

    //@}

    //@{ \name Error methods.

    /*! \brief This method computes the error in orthonormality of a multivector, measured via \f$ \|X^T Op X - I\|_F \f$.
     *  The method has the option of
     *  exploiting a caller-provided \c MX.
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthonormError(const MV &X, Teuchos::RefCountPtr<const MV> MX) const;

    /*! \brief This method computes the error in orthogonality of two multivectors, measure via \f$ \|Q^T Op X \|_F \f$.
     *  The method has the option of
     *  exploiting a caller-provided \c MX.
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthogError(const MV &X1, Teuchos::RefCountPtr<const MV> MX1, const MV &X2) const;

    //@}

  private:
    
    //! Parameter for re-orthogonalization.
    MagnitudeType _kappa;
  
    // ! Routine to find an orthonormal basis for the 
    ReturnType findBasis(MV &X, Teuchos::RefCountPtr<MV> MX, int &rank, bool completeBasis ) const;
    
  };


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthonormality
  template<class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
  BasicOrthoManager<ScalarType,MV,OP>::orthonormError(const MV &X, Teuchos::RefCountPtr<const MV> MX) const {
    const ScalarType ONE = SCT::one();
    int rank = MVT::GetNumberVecs(X);
    Teuchos::SerialDenseMatrix<int,double> xTx(rank,rank);
    ReturnType ret = innerProd(X,X,MX,xTx);
    TEST_FOR_EXCEPTION( ret != Ok, std::runtime_error, "BasicOrthoManager::orthonormError(): innerProd returned error!" );
    for (int i=0; i<rank; i++) {
      xTx(i,i) -= ONE;
    }
    return xTx.normFrobenius();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthogonality
  template<class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
  BasicOrthoManager<ScalarType,MV,OP>::orthogError(const MV &X1, Teuchos::RefCountPtr<const MV> MX1, const MV &X2) const {
    int r1 = MVT::GetNumberVecs(X1);
    int r2  = MVT::GetNumberVecs(X2);
    Teuchos::SerialDenseMatrix<int,double> xTx(r2,r1);
    ReturnType ret = innerProd(X2,X1,MX1,xTx);
    TEST_FOR_EXCEPTION( ret != Ok, std::runtime_error, "BasicOrthoManager::orthogError(): innerProd returned error!" );
    return xTx.normFrobenius();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X) - span(W)
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::projectAndNormalize(
            MV &X, Teuchos::RefCountPtr<MV> MX, const MV &Q, int &rank ) const {

    int qc = MVT::GetNumberVecs( Q );
    int qr = MVT::GetVecLength( Q );
    int xc = MVT::GetNumberVecs( X );
    int xr = MVT::GetVecLength( X );
    ReturnType ret;

    /******   DO NO MODIFY *MX IF _hasOp == false   ******/

    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,MVT::GetNumberVecs(X));
        ret = OPT::Apply(*(this->_Op),X,*MX);
        if (ret != Ok) return Failed;
      }
    }
    else {
      // Op == I  -->  MX = X (ignore it if the user passed it in)
      MX = Teuchos::rcp( &X, false );
    }

    int mxc = MVT::GetNumberVecs( *MX );
    int mxr = MVT::GetVecLength( *MX );

    // short-circuit
    if (xc == 0 || xr == 0) {
      return Ok;
    }

    // check size of X and Q w.r.t. common sense
    if (xc<0 || xr<0 || mxc<0 || mxr<0) {
      return Failed;
    }
    // check size of X w.r.t. MX 
    else if ( xc!=mxc || xr!=mxr ) {
      return Failed;
    }
    // check Q
    else if ( qr!=0 && qc!=0 && qr<qc ) {
      return Failed;
    }
    // check size of Q
    else if ( qr!=0 && qc!=0 && qr!=xr ) {
      return Failed;
    }
    // check feasibility
    else if ( qc+xc > xr ) {
      return Failed;
    }

    // start working
    rank = 0;
    int numTries = 10;   // each vector in X gets 10 random chances to escape degeneracy
    // this is not an infinite loop;
    // either 
    do {

      ReturnType ret;
      int curxsize = xc - rank;
      std::vector<int> curind(curxsize);
      for (int i=0; i<curxsize; i++) {
        curind[i] = rank+i;
      }
      Teuchos::RefCountPtr<MV> curX, curMX;
      curX  = MVT::CloneView(X,curind);
      if (this->_hasOp) {
        curMX = MVT::CloneView(*MX,curind);
      }
      else {
        curMX = curX;
      }

      // orthogonalize current X against Q
      ret = project(*curX,curMX,Q);
      if (ret == Failed) return Failed;

      // orthonormalize X, but quit if it is rank deficient
      ret = findBasis(*curX,curMX,curxsize,false);
      rank += curxsize;
      if (ret == Failed) {
        return Failed;
      }
      else if (ret == Ok) {
        // we are done
        break;
      }
      else {
        TEST_FOR_EXCEPTION( ret != Undefined, std::logic_error, "BasicOrthoManager::projectAndNormalize(): findBasis returned impossible value" );

        if (curxsize > 0) {
          // we finished a vector; reset the chance counter
          numTries = 10;
        }

        // has this vector run out of chances to escape degeneracy?
        if (numTries <= 0) {
          break;
        }
        // use one of this vector's chances
        numTries--;

        // FINISH: CGB: 05/09/2006: have Heidi look at these releases,to reassure me they are sufficient
        curX.release();
        curMX.release();

        // randomize troubled direction
        std::vector<int> ind(1);
        ind[0] = rank;
#ifdef ORTHO_DEBUG
        cout << "Random for column " << rank << endl;
#endif
        curX  = MVT::CloneView(X,ind);
        if (this->_hasOp) {
          curMX = MVT::CloneView(*MX,curind);
        }
        else {
          curMX = curX;
        }
        MVT::MvRandom(*curX);
        if (this->_hasOp) {
          if ( OPT::Apply( *(this->_Op), *curX, *curMX ) != Ok ) return Failed;
        }

        // orthogonalize against Q and prevX
        std::vector<int> prevind(rank);
        for (int i=0; i<rank; i++) {
          prevind[i] = i;
        }
        Teuchos::RefCountPtr<MV> prevX, prevMX;
        prevX  = MVT::CloneView(X,prevind);
        if (this->_hasOp) {
          prevMX = MVT::CloneView(*MX,prevind);
        }
        else {
          prevMX = prevX;
        }
        // orthogonalize curMX against Q
        if ( project(*curX,curMX,Q) != Ok ) return Failed;
        
        // FINISH: CGB: 05/09/2006: have Heidi look at these releases,to reassure me they are sufficient
        curX.release();
        curMX.release();

        // orthogonalize the rest of X against prevX
        curind.resize(xc-rank);
        for (int i=0; i<xc-rank; i++) {
          curind[i] = rank+i;
        }
        curX  = MVT::CloneView(X,curind);
        if (this->_hasOp) {
          curMX = MVT::CloneView(*MX,curind);
        }
        else {
          curMX = curX;
        }
        if ( project(*curX,curMX,*prevX) != Ok ) return Failed;

        // FINISH: CGB: 05/09/2006: have Heidi look at these releases,to reassure me they are sufficient
        curX.release();
        curMX.release();

        continue;
      }

    } while (1);

    if (rank < xc) {
      return Undefined;
    }

    // this should not raise an exception; but our post-conditions oblige us to check
    TEST_FOR_EXCEPTION( rank != xc, std::logic_error, "BasicOrthoManager::projectAndNormalize(): error in rank variable" );
    return Ok; 
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an M-orthonormal basis for span(X), with rank numvectors(X)
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::normalize(MV &X, Teuchos::RefCountPtr<MV> MX, int &rank ) const {
    // call findBasis, with the instruction to try to generate a basis of rank numvecs(X)
    return findBasis(X, MX, rank, true );
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::project(MV &X, Teuchos::RefCountPtr<MV> MX, const MV &Q ) const {
    // For the inner product defined by the operator M or the identity (M == 0)
    //   -> Orthogonalize X against Q
    // Modify MX accordingly
    //
    // Note that when M is 0, MX is not referenced
    //
    // Parameter variables
    //
    // X  : Vectors to be transformed
    //
    // MX : Image of the block vector X by the mass matrix
    //
    // Q  : Vectors to orthogonalize against. These are assumed orthonormal.
    //
    // Return the status of the computation
    //
    // Ok     >> Success
    // Failed >> Failure >> Length of vectors in X != length of vectors in MX
    //                   >> Number of columns of X != number of columns of MX
    //                   >> Length of vectors in X != length of vectors in Q
    //                   >> Failure applying operator M

    ScalarType    ONE  = SCT::one();
    MagnitudeType ZERO = SCT::magnitude(SCT::zero());
    ScalarType    EPS  = SCT::eps();

    int qc = MVT::GetNumberVecs( Q );
    int qr = MVT::GetVecLength( Q );
    int xc = MVT::GetNumberVecs( X );
    int xr = MVT::GetVecLength( X );

    /******   DO NO MODIFY *MX IF _hasOp == false   ******/

    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,MVT::GetNumberVecs(X));
        ReturnType ret = OPT::Apply(*(this->_Op),X,*MX);
        if (ret != Ok) return Failed;
      }
    }
    else {
      // Op == I  -->  MX = X (ignore it if the user passed it in)
      MX = Teuchos::rcp( &X, false );
    }

    int mxc = MVT::GetNumberVecs( *MX );
    int mxr = MVT::GetVecLength( *MX );

    // short-circuit
    if (xc == 0 || xr == 0 || qc == 0 || qr == 0 ) {
      return Ok;
    }

    // check size of X and Q w.r.t. common sense
    if (xc<0 || xr<0 || mxc<0 || mxr<0) {
      return Failed;
    }
    // check size of X w.r.t. MX and Q
    else if ( xc!=mxc || xr!=mxr || xr!=qr ) {
      return Failed;
    }
    // check Q
    else if ( qr<qc ) {
      return Failed;
    }

    std::vector<int> index(xc);
    for (int i=0; i<xc; i++) {
      index[i] = i;
    }

    // Perform the Gram-Schmidt transformation for a block of vectors

    // Compute the initial M-norms
    std::vector<ScalarType> oldDot( xc );
    MVT::MvDot( X, *MX, &oldDot );

    // Define the product Q^T * (M*X)
    // Multiply Q' with MX
    Teuchos::SerialDenseMatrix<int,ScalarType> qTmx( qc, xc );
    innerProd(Q,X,MX,qTmx);
    // Multiply by Q and substract the result in X
    MVT::MvTimesMatAddMv( -ONE, Q, qTmx, ONE, X );

    // Update MX, with the least number of applications of M as possible
    Teuchos::RefCountPtr<MV> MQ;
    if (this->_hasOp) {
      if (xc <= qc) {
        if ( OPT::Apply( *(this->_Op), X, *MX) != Ok ) {
          return Failed;
        }
      }
      else {
        // this will possibly be used again below; don't delete it
        MQ = MVT::Clone( Q, qc );
        if ( OPT::Apply( *(this->_Op), Q, *MQ ) != Ok ) {
          return Failed;
        }
        MVT::MvTimesMatAddMv( -ONE, *MQ, qTmx, ONE, *MX );
      }
    }

    // Compute new M-norms
    std::vector<ScalarType> newDot(xc);
    MVT::MvDot( X, *MX, &newDot );

    // determine (individually) whether to do another step of classical Gram-Schmidt
    for (int j = 0; j < xc; ++j) {
      
      if ( SCT::magnitude(_kappa*newDot[j]) < SCT::magnitude(oldDot[j]) ) {
        
        // Apply another step of classical Gram-Schmidt
        innerProd(Q,X,MX,qTmx);
        MVT::MvTimesMatAddMv( -ONE, Q, qTmx, ONE, X );
        
        // Update MX, with the least number of applications of M as possible
        if (this->_hasOp) {
          if (MQ.get()) {
            // MQ was allocated and computed above; use it
            MVT::MvTimesMatAddMv( -ONE, *MQ, qTmx, ONE, *MX );
          }
          else if (xc <= qc) {
            // MQ was not allocated and computed above; it was cheaper to use X before and it still is
            if ( OPT::Apply( *(this->_Op), X, *MX) != Ok ) { 
              return Failed;
            }
          }
        }
        
        break;
      } // if (_kappa*newDot[j] < oldDot[j])
    } // for (int j = 0; j < xc; ++j)

    return Ok;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an M-orthonormal basis for span(X), with the option of extending the subspace so that 
  // the rank is numvectors(X)
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::findBasis(MV &X, Teuchos::RefCountPtr<MV> MX, int &rank, bool completeBasis ) const {
    // For the inner product defined by the operator M or the identity (M == 0)
    //   -> Orthonormalize X 
    // Modify MX accordingly
    //
    // Note that when M is 0, MX is not referenced
    //
    // Parameter variables
    //
    // X  : Vectors to be orthonormalized
    //
    // MX : Image of the multivector X under the operator M
    //
    // M  : Pointer to the operator for the inner product
    //
    // TODO: add reference
    // kappa= Coefficient determining when to perform a second Gram-Schmidt step
    //        Default value = 1.5625 = (1.25)^2 (as suggested in Parlett's book)
    //
    // Return the status of the computation
    // Ok        >> Success
    // Undefined >> X is rank deficient and completeBasis == false
    // Failed    >> Failure >> X has zero columns or zero rows
    //                      >> Length of vectors in X != length of vectors in MX
    //                      >> Number of columns of X != number of columns of MX
    //                      >> Number of columns of X  > number of rows of X
    //                      >> X is rank deficient and completeBasis == true and we couldn't fix it
    //                      >> Failure applying operator M

    const ScalarType ONE  = SCT::one();
    const ScalarType ZERO = SCT::zero();
    const ScalarType EPS  = SCT::eps();

    int xc = MVT::GetNumberVecs( X );
    int xr = MVT::GetVecLength( X );

    /******   DO NO MODIFY *MX IF _hasOp == false   ******/

    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,MVT::GetNumberVecs(X));
        ReturnType ret = OPT::Apply(*(this->_Op),X,*MX);
        if (ret != Ok) return Failed;
      }
    }
    else {
      // Op == I  -->  MX = X (ignore it if the user passed it in)
      MX = Teuchos::rcp( &X, false );
    }

    int mxc = MVT::GetNumberVecs( *MX );
    int mxr = MVT::GetVecLength( *MX );

    // short-circuit
    if (xc != mxc) {
      return Failed;
    }
    else if (xc > xr) {
      return Failed;
    }
    else if (xr != mxr) {
      return Failed;
    }
    else if (xc == 0 || xr == 0) {
      return Failed;
    }

    for (int j = 0; j < xc; j++) {

      // numX is 
      // * number of currently orthonormal columns of X
      // * the index of the current column of X
      int numX = j;

      // Get a view of the vector currently being worked on.
      std::vector<int> index(1);
      index[0] = numX;
      Teuchos::RefCountPtr<MV> Xj = MVT::CloneView( X, index );
      Teuchos::RefCountPtr<MV> MXj;
      if ((this->_hasOp)) {
        // MXj is a view of the current vector in MX
        MXj = MVT::CloneView( *MX, index );
      }
      else {
        // MXj is a pointer to Xj, and MUST NOT be modified
        MXj = Xj;
      }

      // Get a view of the previous vectors.
      std::vector<int> prev_idx( numX );
      Teuchos::RefCountPtr<MV> prevX;

      if (numX > 0) {
        for (int i=0; i<numX; i++) {
          prev_idx[i] = i;
        }
        prevX = MVT::CloneView( X, prev_idx );
      } 

      // Make storage for these Gram-Schmidt iterations.
      Teuchos::SerialDenseMatrix<int,ScalarType> product( numX, 1 );

      bool rankDef = true;
      for (int numTrials = 0; numTrials < 10; numTrials++) {

        std::vector<ScalarType> oldDot( 1 ), newDot( 1 );

        //
        // Save old MXj vector and compute M-norm
        //
        Teuchos::RefCountPtr<MV> oldMXj = MVT::CloneCopy( *MXj ); 
        MVT::MvDot( *Xj, *oldMXj, &oldDot );
        // Xj^H M Xj should be real and positive, by the hermitian positive definiteness of M
        if ( SCT::real(oldDot[0]) < ZERO ) {
          return Failed;
        }

        if (numX > 0) {
          // Apply the first step of Gram-Schmidt

          // product <- prevX^T MXj
          innerProd(*prevX,*MXj,product);

          // Xj <- Xj - prevX prevX^T MXj   
          //     = Xj - prevX product
          MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *Xj );

          // Update MXj
          Teuchos::RefCountPtr<const MV> prevMX;
          if (this->_hasOp) {
            // prevMX may be used below
            prevMX = MVT::CloneView( *MX, prev_idx );
            // MXj <- M*Xj_new
            //      = M*(Xj_old - prevX prevX^T MXj)
            //      = MXj - prevMX product
            MVT::MvTimesMatAddMv( -ONE, *prevMX, product, ONE, *MXj );
          }

          // Compute new M-norm
          MVT::MvDot( *Xj, *MXj, &newDot );

          // Check if a correction is needed.
          if ( SCT::magnitude(_kappa*newDot[0]) < SCT::magnitude(oldDot[0]) ) {
            // Apply the second step of Gram-Schmidt
            // This is the same as above
            innerProd(*prevX,*MXj,product);
            MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *Xj );
            if ((this->_hasOp)) {
              // prevMX was initialized above
              MVT::MvTimesMatAddMv( -ONE, *prevMX, product, ONE, *MXj );
            }
          } // if (_kappa*newDot[0] < oldDot[0])

        } // if (numX > 0)

        // Compute M-norm with old MXj
        MVT::MvDot( *Xj, *oldMXj, &newDot );

        // Check if Xj has any directional information left after the orthogonalization.
#ifdef ORTHO_DEBUG
        cout << "olddot: " << SCT::magnitude(oldDot[0]) << "    newdot: " << SCT::magnitude(newDot[0]);
#endif
        if ( SCT::magnitude(newDot[0]) > SCT::magnitude(oldDot[0]*EPS*EPS) && SCT::real(newDot[0]) > ZERO ) {
#ifdef ORTHO_DEBUG
          cout << " ACCEPTED" << endl;
#endif
          // Normalize Xj.
          // Xj <- Xj / sqrt(newDot*EPS*EPS)
          MVT::MvAddMv( ONE/SCT::squareroot(SCT::magnitude(newDot[0])), *Xj, ZERO, *Xj, *Xj );
          if (this->_hasOp) {
            // Update MXj.
            MVT::MvAddMv( ONE/SCT::squareroot(SCT::magnitude(newDot[0])), *MXj, ZERO, *MXj, *MXj );
          }
          // We are not rank deficient in this vector. Move on to the next vector in X.
          rankDef = false;
          break;
        }
        else {
#ifdef ORTHO_DEBUG
          cout << " REJECTED" << endl;
#endif
          // There was nothing left in Xj after orthogonalizing against previous columns in X.
          // X is rank deficient.

          if (completeBasis) {
            // Nothing left in Xj. Fill it with random information and keep going.
#ifdef ORTHO_DEBUG
            cout << "Random for column " << j << endl;
#endif
            MVT::MvRandom( *Xj );
            if (this->_hasOp) {
              if ( OPT::Apply( *(this->_Op), *Xj, *MXj ) != Ok ) return Failed;
            }
          }
          else {
            rankDef = true;
            break;
          }

        } // if (norm > oldDot*EPS*EPS)

      }  // for (numTrials = 0; numTrials < 10; ++numTrials)

      // if rankDef == true, then quit and notify user of rank obtained
      if (rankDef == true) {
        rank = j;
        MVT::MvInit( *Xj, ZERO );
        if (this->_hasOp) {
          MVT::MvInit( *MXj, ZERO );
        }
        if (completeBasis) {
          return Failed;
        }
        return Undefined;
      }

    } // for (j = 0; j < xc; ++j)

    rank = xc;
    return Ok;
  }

} // namespace Anasazi

#endif // ANASAZI_BASIC_ORTHOMANAGER_HPP

