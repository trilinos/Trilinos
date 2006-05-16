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


// FINISH: using findBasis may have been more trouble than it is worth; consider rewriting?

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
    ReturnType project ( MV &X, Teuchos::RefCountPtr<MV> MX, 
                         Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > C, 
                         const MV &Q ) const;


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
    ReturnType normalize ( MV &X, Teuchos::RefCountPtr<MV> MX, 
                           Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > R, 
                           int &rank ) const;


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
    ReturnType projectAndNormalize ( MV &X, Teuchos::RefCountPtr<MV> MX, 
                                     Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > C, 
                                     Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > R, 
                                     const MV &Q, int &rank ) const;

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
    ReturnType findBasis(MV &X, Teuchos::RefCountPtr<MV> MX, 
                         Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > C, 
                         int &rank, bool completeBasis, int howMany = -1 ) const;
    
  };


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthonormality
  template<class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
  BasicOrthoManager<ScalarType,MV,OP>::orthonormError(const MV &X, Teuchos::RefCountPtr<const MV> MX) const {
    const ScalarType ONE = SCT::one();
    int rank = MVT::GetNumberVecs(X);
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(rank,rank);
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
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(r2,r1);
    ReturnType ret = innerProd(X2,X1,MX1,xTx);
    TEST_FOR_EXCEPTION( ret != Ok, std::runtime_error, "BasicOrthoManager::orthogError(): innerProd returned error!" );
    return xTx.normFrobenius();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X) - span(W)
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::projectAndNormalize(
                                    MV &X, Teuchos::RefCountPtr<MV> MX, 
                                    Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > C, 
                                    Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > R, 
                                    const MV &Q, int &rank ) const {

    int qc = MVT::GetNumberVecs( Q );
    int qr = MVT::GetVecLength( Q );
    int xc = MVT::GetNumberVecs( X );
    int xr = MVT::GetVecLength( X );
    ReturnType ret;


    /* if the user doesn't want to store the coefficienets, 
     * allocate some local memory for them 
     */
    if ( C == Teuchos::null ) {
      C = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(qc,xc) );
    }
    if ( R == Teuchos::null ) {
      R = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xc,xc) );
    }

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

    // check size of C, R
    if ( C->numRows() != qc || C->numCols() != xc || 
         R->numRows() != xc || R->numCols() != xc ) {
      return Failed;
    }
    // check size of X and Q w.r.t. common sense
    else if (xc<0 || xr<0 || mxc<0 || mxr<0) {
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

    // orthogonalize all of X against Q
    ret = project(X,MX,C,Q);
    if (ret == Failed) return Failed;

    Teuchos::SerialDenseMatrix<int,ScalarType> oldCoeff(xc,1);

    // start working
    rank = 0;
    int numTries = 10;   // each vector in X gets 10 random chances to escape degeneracy
    int oldrank = -1;
    do {

      int curxsize = xc - rank;

      // orthonormalize X, but quit if it is rank deficient
      // we can't let findBasis generated random vectors to complete the basis,
      // because it doesn't know about Q; we will do this ourselves below
      ret = findBasis(X,MX,R,rank,false,curxsize);

      if (rank < xc && numTries == 10) {
        // we quit on this vector, and for the first time;
        // save the coefficient information, because findBasis will overwrite it
        for (int i=0; i<xc; i++) {
          oldCoeff(i,0) = (*R)(i,rank);
        }
      }

      if (oldrank != -1 && rank != oldrank) {
        // we moved on; restore the previous coefficients
        for (int i=0; i<xc; i++) {
          (*R)(i,oldrank) = oldCoeff(i,0);
        }
      }

      if (ret == Failed) {
        return Failed;
      }
      else if (ret == Ok) {
        // we are done
        break;
      }
      else {
        TEST_FOR_EXCEPTION( ret != Undefined, std::logic_error, "BasicOrthoManager::projectAndNormalize(): findBasis returned impossible value" );
        TEST_FOR_EXCEPTION( rank < oldrank, std::logic_error,   "BasicOrthoManager::projectAndNormalize(): basis lost rank; this shouldn't happen");

        if (rank != oldrank) {
          // we added a basis vector from random info; reset the chance counter
          numTries = 10;
        }

        // store old rank
        oldrank = rank;

        // has this vector run out of chances to escape degeneracy?
        if (numTries <= 0) {
          break;
        }
        // use one of this vector's chances
        numTries--;

        // randomize troubled direction
#ifdef ORTHO_DEBUG
        cout << "Random for column " << rank << endl;
#endif
        Teuchos::RefCountPtr<MV> curX, curMX;
        std::vector<int> ind(1);
        ind[0] = rank;
        curX  = MVT::CloneView(X,ind);
        MVT::MvRandom(*curX);
        if (this->_hasOp) {
          curMX = MVT::CloneView(*MX,ind);
          if ( OPT::Apply( *(this->_Op), *curX, *curMX ) != Ok ) return Failed;
        }

        // orthogonalize against Q
        // if !this->_hasOp, the curMX will be ignored.
        // we don't care about these coefficients; in fact, we need to preserve the previous coeffs
        if ( project(*curX,curMX,Teuchos::null,Q) != Ok ) return Failed;
        
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
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::normalize(
                                MV &X, Teuchos::RefCountPtr<MV> MX, 
                                Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > R, 
                                int &rank ) const {
    // call findBasis, with the instruction to try to generate a basis of rank numvecs(X)
    return findBasis(X, MX, R, rank, true );
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::project(
                                MV &X, Teuchos::RefCountPtr<MV> MX, 
                                Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > C, 
                                const MV &Q ) const {
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

    /* if the user doesn't want to store the coefficienets, 
     * allocate some local memory for them 
     */
    if ( C == Teuchos::null ) {
      C = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(qc,xc) );
    }

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

    // check size of C, R
    if ( C->numRows() != qc || C->numCols() != xc ) {
      return Failed;
    }
    // check size of X and Q w.r.t. common sense
    else if (xc<0 || xr<0 || mxc<0 || mxr<0) {
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
    innerProd(Q,X,MX,*C);
    // Multiply by Q and subtract the result in X
    MVT::MvTimesMatAddMv( -ONE, Q, *C, ONE, X );

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
        MVT::MvTimesMatAddMv( -ONE, *MQ, *C, ONE, *MX );
      }
    }

    // Compute new M-norms
    std::vector<ScalarType> newDot(xc);
    MVT::MvDot( X, *MX, &newDot );

    // determine (individually) whether to do another step of classical Gram-Schmidt
    for (int j = 0; j < xc; ++j) {
      
      if ( SCT::magnitude(_kappa*newDot[j]) < SCT::magnitude(oldDot[j]) ) {

        Teuchos::SerialDenseMatrix<int,ScalarType> C2(*C);
        
        // Apply another step of classical Gram-Schmidt
        innerProd(Q,X,MX,C2);
        *C += C2;
        MVT::MvTimesMatAddMv( -ONE, Q, C2, ONE, X );
        
        // Update MX, with the least number of applications of M as possible
        if (this->_hasOp) {
          if (MQ.get()) {
            // MQ was allocated and computed above; use it
            MVT::MvTimesMatAddMv( -ONE, *MQ, C2, ONE, *MX );
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
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::findBasis(
                MV &X, Teuchos::RefCountPtr<MV> MX, 
                Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > R,
                int &rank, bool completeBasis, int howMany ) const {
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
    const MagnitudeType ZERO = SCT::magnitude(SCT::zero());
    const ScalarType EPS  = SCT::eps();

    int xc = MVT::GetNumberVecs( X );
    int xr = MVT::GetVecLength( X );

    if (howMany == -1) {
      howMany = xc;
    }

    /*******************************************************
     *  If _hasOp == false, we will not reference MX below *
     *******************************************************/

    // if Op==null, MX == X (via pointer)
    // Otherwise, either the user passed in MX or we will allocated and compute it
    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,xc);
        ReturnType ret = OPT::Apply(*(this->_Op),X,*MX);
        if (ret != Ok) return Failed;
      }
    }

    /* if the user doesn't want to store the coefficienets, 
     * allocate some local memory for them 
     */
    if ( R == Teuchos::null ) {
      R = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xc,xc) );
    }

    int mxc = (this->_hasOp) ? MVT::GetNumberVecs( *MX ) : xc;
    int mxr = (this->_hasOp) ? MVT::GetVecLength( *MX )  : xr;

    // check size of C, R
    if ( R->numRows() != xc || R->numCols() != xc ) {
      return Failed;
    }
    else if (xc != mxc) {
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
    else if (howMany < 0 || howMany > xc) {
      return Failed;
    }

    /* xstart is which column we are starting the process with, based on howMany
     * columns before xstart are assumed to be M-orthonormal already
     */
    int xstart = xc - howMany;

    for (int j = xstart; j < xc; j++) {

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
      Teuchos::RefCountPtr<const MV> prevX, prevMX;

      if (numX > 0) {
        for (int i=0; i<numX; i++) {
          prev_idx[i] = i;
        }
        prevX = MVT::CloneView( X, prev_idx );
        if (this->_hasOp) {
          prevMX = MVT::CloneView( *MX, prev_idx );
        }
      } 

      bool rankDef = true;
      /* numTrials>0 will denote that the current vector was randomized for the purpose
       * of finding a basis vector, and that the coefficients of that vector should
       * not be stored in R
       */
      for (int numTrials = 0; numTrials < 10; numTrials++) {

        // Make storage for these Gram-Schmidt iterations.
        Teuchos::SerialDenseMatrix<int,ScalarType> product(numX, 1);
        std::vector<ScalarType> oldDot( 1 ), newDot( 1 );

        //
        // Save old MXj vector and compute M-norm
        //
        Teuchos::RefCountPtr<MV> oldMXj = MVT::CloneCopy( *MXj ); 
        MVT::MvDot( *Xj, *MXj, &oldDot );
        // Xj^H M Xj should be real and positive, by the hermitian positive definiteness of M
        if ( SCT::real(oldDot[0]) < ZERO ) {
          return Failed;
        }

        if (numX > 0) {
          // Apply the first step of Gram-Schmidt

          // product <- prevX^T MXj
          innerProd(*prevX,*prevX,MXj,product);

          // Xj <- Xj - prevX prevX^T MXj   
          //     = Xj - prevX product
          MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *Xj );

          // Update MXj
          if (this->_hasOp) {
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
            Teuchos::SerialDenseMatrix<int,ScalarType> P2(numX,1);

            innerProd(*prevX,*prevX,MXj,P2);
            product += P2;
            MVT::MvTimesMatAddMv( -ONE, *prevX, P2, ONE, *Xj );
            if ((this->_hasOp)) {
              MVT::MvTimesMatAddMv( -ONE, *prevMX, P2, ONE, *MXj );
            }
          } // if (_kappa*newDot[0] < oldDot[0])

        } // if (numX > 0)

        // Compute M-norm with old MXj
        MVT::MvDot( *Xj, *oldMXj, &newDot );

        // save the coefficients, if we are working on the original vector and not a randomly generated one
        if (numTrials == 0) {
          for (int i=0; i<numX; i++) {
            (*R)(i,j) = product(i,0);
          }
        }

        // Check if Xj has any directional information left after the orthogonalization.
#ifdef ORTHO_DEBUG
        cout << "olddot: " << SCT::magnitude(oldDot[0]) << "    newdot: " << SCT::magnitude(newDot[0]);
#endif
        if ( SCT::magnitude(newDot[0]) > SCT::magnitude(oldDot[0]*EPS*EPS) && SCT::real(newDot[0]) > ZERO ) {
#ifdef ORTHO_DEBUG
          cout << " ACCEPTED" << endl;
#endif
          // Normalize Xj.
          // Xj <- Xj / sqrt(newDot)
          ScalarType diag = SCT::squareroot(SCT::magnitude(newDot[0]));

          MVT::MvAddMv( ONE/diag, *Xj, ZERO, *Xj, *Xj );
          if (this->_hasOp) {
            // Update MXj.
            MVT::MvAddMv( ONE/diag, *MXj, ZERO, *MXj, *MXj );
          }

          // save it, if it corresponds to the original vector and not a randomly generated one
          if (numTrials == 0) {
            (*R)(j,j) = diag;
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
          // reflect this in the coefficients
          (*R)(j,j) = ZERO;

          if (completeBasis) {
            // Fill it with random information and keep going.
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
        (*R)(j,j) = ZERO;
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

