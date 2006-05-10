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
      \brief An implementation of the Anasazi::OrthoManager that performs orthogonalization
      using (potentially) multiple steps of classical Gram-Schmidt.
      
      \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziOrthoManager.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Anasazi {

  template<class ScalarType, class MV, class OP>
  class BasicOrthoManager : public OrthoManager<ScalarType,MV,OP> {

  private:
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef typename Teuchos::ScalarTraits<ScalarType>  SCT;
    typedef typename MultiVecTraits<ScalarType,MV>      MVT;
    typedef typename OperatorTraits<ScalarType,MV,OP>   OPT;

  public:
    
    //@{ \name Constructor/Destructor.
    //! Constructor specifying re-orthogonalization tolerance.
    BasicOrthoManager( const MagnitudeType kappa = SCT::magnitude(1.5625) ) { _kappa = kappa; };

    //! Destructor
    ~BasicOrthoManager() {};
    //@}


    //@{ \name Accessor routines

    //! Set parameter for re-orthogonalization threshhold.
    void setKappa( const MagnitudeType kappa ) { _kappa = kappa; };

    //! Return parameter for re-orthogonalization threshhold.
    void getKappa() { return kappa; }

    //@} 


    //@{ \name Orthogonalization methods.

    /*! \brief This method takes a multivector and projects it onto the space orthogonal to 
     *  another given multivector, in a specified inner product.
     *  
     @param X [in/out] The multivector to the modified. 
       On output, \f$X^T M Q = 0\f$ if \c M!=0 or \f$X^T Q = 0\f$ if \c M==0.

     @param MX [in/out] The image of \c X under the operator \c M. 
       On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \c M == 0, \c MX is not referenced.

     @param M [in] Pointer to the operator used for the inner product. If \c M == 0, the Euclidean inner product is used.

     @param Q [in] A multivector specifying the space to be orthogonalized against. 

     @return Code specifying failure of the routine, as defined by the implementation.
    */
    ReturnType project ( MV &X, MV &MX, const OP *M, const MV &Q ) const;



    /*! \brief This method takes a multivector and orthonormalizes the columns using a specified inner product. It extends heroic measures 
     *  to find a basis with as many vectors as the input multivector, in order to avoid returning \c Anasaiz::Undefined.
     *  
     @param X [in/out] The multivector to the modified. 
       On output, \f$X^T M X = I\f$ if \c M!=0 or \f$X^T X = I\f$ if \c M==0.

     @param MX [in/out] The image of \c X under the operator \c M. 
       On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \c M == 0, \c MX is not referenced.

     @param M [in] Pointer to the operator used for the inner product. If \c M == 0, the Euclidean inner product is used.

     @param rank [out] Rank of the basis computed by this method.

     @return Code specifying failure of the routine, as defined by the implementation.
    */
    ReturnType normalize ( MV &X, MV &MX, const OP *M, int &rank ) const;



    /*! \brief This method takes a multivector and projects it onto the space orthogonal to 
     *  another given multivector, in a specified inner product. It also orthonormalizes the 
     *  columns of the resulting multivector. It extends heroic measures to find a basis with as many vectors as the input multivector, in
     *  order to avoid returning \c Anasazi::Undefined.
     *  
     @param X [in/out] The multivector to the modified. 
       On output, \f$X^T M Q = 0\f$ and \f$X^T M X = I\f$ if \c M!=0 or \f$X^T Q = 0\f$ and \f$X^T X = I\f$ if \c M==0.

     @param MX [in/out] The image of \c X under the operator \c M. 
       On input, this is expected to be consistent with \c X. On output, this is updated consistent with updates to \c X.
       If \c M == 0, \c MX is not referenced.

     @param M [in] Pointer to the operator used for the inner product. If \c M == 0, the Euclidean inner product is used.

     @param Q [in] A multivector specifying the space to be orthogonalized against. 

     @param rank [out] Rank of the basis computed by this method.

     @return Code specifying failure of the routine, as defined by the implementation.
    */
    ReturnType projectAndNormalize ( MV &X, MV &MX, const OP *M, const MV &Q, int &rank ) const;

    //@}

  private:
    
    //! Parameter for re-orthogonalization.
    MagnitudeType _kappa;
  
    // ! Routine to find an orthonormal basis for the 
    static ReturnType BasicOrthoManager<ScalarType, MV, OP>::findBasis(MV &X, MV &MX, const OP *M, int &rank, bool completeBasis ) const;
    
  };



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an M-orthonormal basis for span(X) - span(W)
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::projectAndNormalize(MV &X, MV &MX, const OP *M, const MV &Q, int &rank ) const {

    int qc = MVT::GetNumberVecs( Q );
    int qr = MVT::GetVecLength( Q );
    int xc = MVT::GetNumberVecs( X );
    int xr = MVT::GetVecLength( X );
    int mxc = (M) ? MVT::GetNumberVecs( MX ) : xc;
    int mxr = (M) ? MVT::GetVecLength( MX ) : xr;

    // short-circuit
    if (xc == 0 || xr == 0) {
      return Ok;
    }

    // check size of X and Q w.r.t. common sense
    if (xc<0 || xr<0 || mc<0 || mr<0) {
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
    do {

      ReturnType ret;
      int curxsize = xc - rank;
      std::vector<int> curind(curxsize);
      for (int i=0; i<curxsize; i++) {
        curind[i] = rank+i;
      }
      Teuchos::RefCountPtr<MV> curX, curMX;
      curX = MVT::CloneView(X,curind);
      if (M) {
        curMX = MVT::CloneView(MX,curind);
      }
      else {
        curMX = curX; // this won't be referenced
      }

      // orthogonalize current X against Q
      ret = project(*curX,*curMX,M,Q);
      if (ret == Failed) return Failed;

      // orthonormalize X, but quit if it is rank deficient
      ret = findBasis(*curX,*curMX,M,rank,false);
      if (ret == Failed) {
        return Failed;
      }
      else if (ret == Ok) {
        // reset the random-chance counter
        numTries = 10;
        continue;
      }
      else {
        TEST_FOR_EXCEPTION( ret != Undefined, std::logic_error, "BasicOrthoManager::projectAndNormalize(): findBasis returned impossible value" );
        // has this vector run out of chances to escape degeneracy?
        if (numTries <= 0) {
          break;
        }
        // steal one
        numTries--;

        // FINISH: CGB: 05/09/2006: have Heidi look at these releases,to reassure me they are sufficient
        curX.release();
        curMX.release();

        // randomize troubled direction
        std::vector<int> ind(1);
        ind[0] = rank;
        curX = MVT::CloneView(X,ind);
        MVT::MvRandom(*curX);
        if (M) {
          curMX = MVT::CloneView(MX,ind);
          if ( OPT::Apply( *M, *curX, *curMX ) != Ok ) return Failed;
        }
        else {
          curMX = curX;
        }

        // orthogonalize against Q and prevX
        std::vector<int> prevind(rank);
        for (int i=0; i<rank; i++) {
          prevind[i] = i;
        }
        Teuchos::RefCountPtr<MV> prevX = MVT::CloneView(X,prevind);
        if ( project(*curX,*curMX,M,Q) != Ok ) return Failed;
        if ( project(*curX,*curMX,M,*prevX) != Ok ) return Failed;

        continue;
      }

    } while (rank < xc);

    if (cdone < xc) {
      rank = cdone;
      return Undefined;
    }
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an M-orthonormal basis for span(X), with rank numvectors(X)
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::normalize(MV &X, MV &MX, const OP *M, int &rank ) const {
    // call findBasis, with the instruction to try to generate a basis of rank numvecs(X)
    return findBasis(X, MX, M, rank, true );
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::project(MV &X, MV &MX, const OP *M, const MV &Q ) const {
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
    // M  : Pointer to the operator used for the inner product
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

    ScalarType ONE  = SCT::one();
    ScalarType ZERO = SCT::zero();
    ScalarType EPS  = SCT::eps();

    int qc = MVT::GetNumberVecs( Q );
    int qr = MVT::GetVecLength( Q );
    int xc = MVT::GetNumberVecs( X );
    int xr = MVT::GetVecLength( X );
    int mxc = (M) ? MVT::GetNumberVecs( MX ) : xc;
    int mxr = (M) ? MVT::GetVecLength( MX ) : xr;

    // short-circuit
    if (xc == 0 || xr == 0 || qc == 0 || qr == 0 ) {
      return Ok;
    }

    // check size of X and Q w.r.t. common sense
    if (xc<0 || xr<0 || mc<0 || mr<0) {
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
    if (M) {
      MVT::MvDot( X, MX, &oldDot );
    }
    else {
      MVT::MvDot( X, X, &oldDot );
    }

    // Define the product Q^T * (M*X)
    // Multiply Q' with MX
    Teuchos::SerialDenseMatrix<int,ScalarType> qTmx( qc, xc );
    if (M) {
      MVT::MvTransMv( ONE, Q, MX, qTmx );
    }
    else {
      MVT::MvTransMv( ONE, Q, X, qTmx );
    }
    // Multiply by Q and substract the result in X
    MVT::MvTimesMatAddMv( -ONE, Q, qTmx, ONE, *X );

    // Update MX, with the least number of applications of M as possible
    Teuchos::RefCountPtr<MV> MQ;
    if (M) {
      if (xc <= qc) {
        if ( OPT::Apply( *M, X, MX) != Ok ) {
          return Failed;
        }
      }
      else {
        // this will possibly be used again below; don't delete it
        MQ = MVT::Clone( Q, qc );
        if ( OPT::Apply( *M, Q, *MQ ) != Ok ) {
          return Failed;
        }
        MVT::MvTimesMatAddMv( -ONE, *MQ, qTmx, ONE, *MXX );
      }
    }

    // Compute new M-norms
    std::vector<ScalarType> newDot(xc);
    MVT::MvDot( *X, *MXX, &newDot );
    
    // determine (individually) whether to do another step of classical Gram-Schmidt
    for (int j = 0; j < xc; ++j) {
      
      if ( SCT::magnitude(kappa*newDot[j]) < SCT::magnitude(oldDot[j]) ) {
        
        // Apply another step of classical Gram-Schmidt
        if (M) {
          MVT::MvTransMv( ONE, Q, MX, qTmx );
        }
        else {
          MVT::MvTransMv( ONE, Q, X, qTmx );
        }
        MVT::MvTimesMatAddMv( -ONE, Q, qTmx, ONE, X );
        
        // Update MX, with the least number of applications of M as possible
        if (M) {
          if (MQ.get()) {
            // MQ was allocated and computed above; use it
            MVT::MvTimesMatAddMv( -ONE, *MQ, qTmx, ONE, MX );
          }
          else if (xc <= qc) {
            // MQ was not allocated and computed above; it was cheaper to use X before and it still is
            if ( OPT::Apply( *M, X, MX) != Ok ) { 
              return Failed;
            }
          }
        }
        
        break;
      } // if (kappa*newDot[j] < oldDot[j])
    } // for (int j = 0; j < xc; ++j)

    return Ok;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an M-orthonormal basis for span(X), with the option of extending the subspace so that 
  // the rank is numvectors(X)
  template<class ScalarType, class MV, class OP>
  ReturnType BasicOrthoManager<ScalarType, MV, OP>::findBasis(MV &X, MV &MX, const OP *M, int &rank, bool completeBasis ) const {
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

    ScalarType ONE  = SCT::one();
    ScalarType ZERO = SCT::zero();
    ScalarType EPS  = SCT::eps();

    int xc = MVT::GetNumberVecs( X );
    int xr = MVT::GetVecLength( X );
    int mxc = (M) ? MVT::GetNumberVecs( MX ) : xc;
    int mxr = (M) ? MVT::GetVecLength( MX ) : xr;

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
      if (M) {
        // MXj is a view of the current vector in MX
        MXj = MVT::CloneView( MX, index );
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

        if (numX > 0) {
          // Apply the first step of Gram-Schmidt

          // product <- prevX^T MXj
          MVT::MvTransMv( ONE, *prevX, *MXj, product );

          // Xj <- Xj - prevX prevX^T MXj   
          //     = Xj - prevX product
          MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *Xj );

          // Update MXj
          Teuchos::RefCountPtr<MV> prevMX;
          if (M) {
            // prevMX may be used below
            prevMX = MVT::CloneView( MX, prev_idx );
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
            MVT::MvTransMv( ONE, *prevX, *MXj, product );
            MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *Xj );
            if (M) {
              // prevMX was initialized above
              MVT::MvTimesMatAddMv( -ONE, *prevMX, product, ONE, *MXj );
            }
          } // if (_kappa*newDot[0] < oldDot[0])

        } // if (numX > 0)

        // Compute M-norm with old MXj
        // NOTE:  Re-using newDot vector
        MVT::MvDot( *Xj, *oldMXj, &newDot );

        // Check if Xj has any directional information left after the orthogonalization.
        if ( SCT::magnitude(newDot[0]) > SCT::magnitude(oldDot[0]*EPS*EPS) ) {
          // Normalize Xj.
          // Xj <- Xj / sqrt(newDot*EPS*EPS)
          MVT::MvAddMv( ONE/SCT::squareroot(newDot[0]), *Xj, ZERO, *Xj, *Xj );
          if (M) {
            // Update MXj.
            MVT::MvAddMv( ONE/SCT::squareroot(newDot[0]), *MXj, ZERO, *MXj, *MXj );
          }
          // We are not rank deficient in this vector. Move on to the next vector in X.
          rankDef = false;
          break;
        }
        else {
          // There was nothing left in Xj after orthogonalizing against previous columns in X.
          // X is rank deficient.

          if (completeBasis) {
            // Nothing left in Xj. Fill it with random information and keep going.
            MVT::MvRandom( *Xj );
            if (M) {
              if ( OPT::Apply( *M, *Xj, *MXj ) != Ok ) return Failed;
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
        if (M) {
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

