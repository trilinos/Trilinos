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

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziMatOrthoManager.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
#  include <Teuchos_FancyOStream.hpp>
#endif

namespace Anasazi {

  template<class ScalarType, class MV, class OP>
  class BasicOrthoManager : public MatOrthoManager<ScalarType,MV,OP> {

  private:
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<ScalarType>  SCT;
    typedef MultiVecTraits<ScalarType,MV>      MVT;
    typedef MultiVecTraitsExt<ScalarType,MV>   MVText;
    typedef OperatorTraits<ScalarType,MV,OP>   OPT;

  public:

    //! @name Constructor/Destructor
    //@{ 
    //! Constructor specifying re-orthogonalization tolerance.
    BasicOrthoManager( Teuchos::RCP<const OP> Op = Teuchos::null, 
                       typename Teuchos::ScalarTraits<ScalarType>::magnitudeType kappa = 1.41421356 /* sqrt(2) */,
                       typename Teuchos::ScalarTraits<ScalarType>::magnitudeType eps = 0.0,
                       typename Teuchos::ScalarTraits<ScalarType>::magnitudeType tol = 0.20 );


    //! Destructor
    ~BasicOrthoManager() {}
    //@}


    //! @name Methods implementing Anasazi::MatOrthoManager
    //@{ 


    /*! \brief Given a list of mutually orthogonal and internally orthonormal bases \c Q, this method
     * projects a multivector \c X onto the space orthogonal to the individual <tt>Q[i]</tt>, 
     * optionally returning the coefficients of \c X for the individual <tt>Q[i]</tt>. All of this is done with respect
     * to the inner product innerProd().
     *
     * After calling this routine, \c X will be orthogonal to each of the <tt>Q[i]</tt>.
     *
     @param X [in/out] The multivector to be modified.<br>
       On output, the columns of \c X will be orthogonal to each <tt>Q[i]</tt>, satisfying
       \f[
       X_{out} = X_{in} - \sum_i Q[i] \langle Q[i], X_{in} \rangle
       \f]

     @param MX [in/out] The image of \c X under the inner product operator \c Op. 
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c Op \cdot X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.

     @param C [out] The coefficients of \c X in the bases <tt>Q[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c X and <tt>Q[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>, similar to calling
       \code
          innerProd( Q[i], X, C[i] );
       \endcode
       If <tt>C[i]</tt> points to a Teuchos::SerialDenseMatrix with size
       inconsistent with \c X and \c <tt>Q[i]</tt>, then a std::invalid_argument
       exception will be thrown. Otherwise, if <tt>C.size() < i</tt> or
       <tt>C[i]</tt> is a null pointer, the caller will not have access to the
       computed coefficients.

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against, satisfying 
     \f[
        \langle Q[i], Q[j] \rangle = I \quad\textrm{if}\quad i=j
     \f]
     and
     \f[
        \langle Q[i], Q[j] \rangle = 0 \quad\textrm{if}\quad i \neq j\ .
     \f]
    */
    void projectMat ( 
          MV &X, 
          Teuchos::Array<Teuchos::RCP<const MV> > Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
          Teuchos::RCP<MV> MX                        = Teuchos::null,
          Teuchos::Array<Teuchos::RCP<const MV> > MQ = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null))
        ) const;


    /*! \brief This method takes a multivector \c X and attempts to compute an orthonormal basis for \f$colspan(X)\f$, with respect to innerProd().
     *
     * The method uses classical Gram-Schmidt with selective reorthogonalization. As a result, the coefficient matrix \c B is upper triangular.
     *
     * This routine returns an integer \c rank stating the rank of the computed basis. If \c X does not have full rank and the normalize() routine does 
     * not attempt to augment the subspace, then \c rank may be smaller than the number of columns in \c X. In this case, only the first \c rank columns of 
     * output \c X and first \c rank rows of \c B will be valid.
     *  
     * The method attempts to find a basis with dimension equal to the number of columns in \c X. It does this by augmenting linearly dependent 
     * vectors in \c X with random directions. A finite number of these attempts will be made; therefore, it is possible that the dimension of the 
     * computed basis is less than the number of vectors in \c X.
     *
     @param X [in/out] The multivector to be modified.<br>
       On output, the first \c rank columns of \c X satisfy
       \f[
          \langle X[i], X[j] \rangle = \delta_{ij}\ .
       \f]
       Also, 
       \f[
          X_{in}(1:m,1:n) = X_{out}(1:m,1:rank) B(1:rank,1:n)
       \f]
       where \c m is the number of rows in \c X and \c n is the number of columns in \c X.

     @param MX [in/out] The image of \c X under the inner product operator \c Op. 
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c Op \cdot X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.

     @param B [out] The coefficients of the original \c X with respect to the computed basis. If \c B is a non-null pointer and \c B matches the dimensions of \c B, then the
     coefficients computed during the orthogonalization routine will be stored in \c B, similar to calling 
       \code
          innerProd( Xout, Xin, B );
       \endcode
     If \c B points to a Teuchos::SerialDenseMatrix with size inconsistent with \c X, then a std::invalid_argument exception will be thrown. Otherwise, if \c B is null, the caller will not have
     access to the computed coefficients. This matrix is not necessarily triangular (as in a QR factorization); see the documentation of specific orthogonalization managers.<br>
     The first rows in \c B corresponding to the valid columns in \c X will be upper triangular.

     @return Rank of the basis computed by this method, less than or equal to the number of columns in \c X. This specifies how many columns in the returned \c X and rows in the returned \c B are valid.
    */
    int normalizeMat ( 
          MV &X, 
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null,
          Teuchos::RCP<MV> MX                                         = Teuchos::null
        ) const;


    /*! \brief Given a set of bases <tt>Q[i]</tt> and a multivector \c X, this method computes an orthonormal basis for \f$colspan(X) - \sum_i colspan(Q[i])\f$.
     *
     *  This routine returns an integer \c rank stating the rank of the computed basis. If the subspace \f$colspan(X) - \sum_i colspan(Q[i])\f$ does not 
     *  have dimension as large as the number of columns of \c X and the orthogonalization manager doe not attempt to augment the subspace, then \c rank 
     *  may be smaller than the number of columns of \c X. In this case, only the first \c rank columns of output \c X and first \c rank rows of \c B will 
     *  be valid.
     *
     * The method attempts to find a basis with dimension the same as the number of columns in \c X. It does this by augmenting linearly dependent 
     * vectors with random directions. A finite number of these attempts will be made; therefore, it is possible that the dimension of the 
     * computed basis is less than the number of vectors in \c X.
     *
     @param X [in/out] The multivector to be modified.<br>
       On output, the first \c rank columns of \c X satisfy
       \f[
            \langle X[i], X[j] \rangle = \delta_{ij} \quad \textrm{and} \quad \langle X, Q[i] \rangle = 0\ .
       \f]
       Also, 
       \f[
          X_{in}(1:m,1:n) = X_{out}(1:m,1:rank) B(1:rank,1:n) + \sum_i Q[i] C[i]
       \f]
       where \c m is the number of rows in \c X and \c n is the number of columns in \c X.

     @param MX [in/out] The image of \c X under the inner product operator \c Op. 
       If \f$ MX != 0\f$: On input, this is expected to be consistent with \c Op \cdot X. On output, this is updated consistent with updates to \c X.
       If \f$ MX == 0\f$ or \f$ Op == 0\f$: \c MX is not referenced.

     @param C [out] The coefficients of \c X in the <tt>Q[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c X and <tt>Q[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>, similar to calling
       \code
          innerProd( Q[i], X, C[i] );
       \endcode
       If <tt>C[i]</tt> points to a Teuchos::SerialDenseMatrix with size
       inconsistent with \c X and \c <tt>Q[i]</tt>, then a std::invalid_argument
       exception will be thrown. Otherwise, if <tt>C.size() < i</tt> or
       <tt>C[i]</tt> is a null pointer, the caller will not have access to the
       computed coefficients.

     @param B [out] The coefficients of the original \c X with respect to the computed basis. If \c B is a non-null pointer and \c B matches the dimensions of \c B, then the
     coefficients computed during the orthogonalization routine will be stored in \c B, similar to calling 
       \code
          innerProd( Xout, Xin, B );
       \endcode
     If \c B points to a Teuchos::SerialDenseMatrix with size inconsistent with \c X, then a std::invalid_argument exception will be thrown. Otherwise, if \c B is null, the caller will not have
     access to the computed coefficients. This matrix is not necessarily triangular (as in a QR factorization); see the documentation of specific orthogonalization managers.<br>
     The first rows in \c B corresponding to the valid columns in \c X will be upper triangular.

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against, satisfying 
     \f[
        \langle Q[i], Q[j] \rangle = I \quad\textrm{if}\quad i=j
     \f]
     and
     \f[
        \langle Q[i], Q[j] \rangle = 0 \quad\textrm{if}\quad i \neq j\ .
     \f]

     @return Rank of the basis computed by this method, less than or equal to the number of columns in \c X. This specifies how many columns in the returned \c X and rows in the returned \c B are valid.

    */
    int projectAndNormalizeMat ( 
          MV &X, 
          Teuchos::Array<Teuchos::RCP<const MV> >  Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B                  = Teuchos::null, 
          Teuchos::RCP<MV> MX                        = Teuchos::null,
          Teuchos::Array<Teuchos::RCP<const MV> > MQ = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null))
        ) const;

    //@}

    //! @name Error methods
    //@{ 

    /*! \brief This method computes the error in orthonormality of a multivector, measured
     * as the Frobenius norm of the difference <tt>innerProd(X,Y) - I</tt>.
     *  The method has the option of exploiting a caller-provided \c MX.
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthonormErrorMat(const MV &X, Teuchos::RCP<const MV> MX = Teuchos::null) const;

    /*! \brief This method computes the error in orthogonality of two multivectors, measured
     * as the Frobenius norm of <tt>innerProd(X,Y)</tt>.
     *  The method has the option of exploiting a caller-provided \c MX.
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthogErrorMat(const MV &X1, const MV &X2, Teuchos::RCP<const MV> MX1, Teuchos::RCP<const MV> MX2) const;

    //@}

    //! @name Accessor routines
    //@{ 

    //! Set parameter for re-orthogonalization threshold.
    void setKappa( typename Teuchos::ScalarTraits<ScalarType>::magnitudeType kappa ) { kappa_ = kappa; }

    //! Return parameter for re-orthogonalization threshold.
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType getKappa() const { return kappa_; } 

    //@} 

  private:
    //! Parameter for re-orthogonalization.
    MagnitudeType kappa_;
    MagnitudeType eps_;
    MagnitudeType tol_;
  
    // ! Routine to find an orthonormal basis
    int findBasis(MV &X, Teuchos::RCP<MV> MX, 
                         Teuchos::SerialDenseMatrix<int,ScalarType> &B, 
                         bool completeBasis, int howMany = -1 ) const;

    //
    // Internal timers
    //
    Teuchos::RCP<Teuchos::Time> timerReortho_;
    
  };


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template<class ScalarType, class MV, class OP>
  BasicOrthoManager<ScalarType,MV,OP>::BasicOrthoManager( Teuchos::RCP<const OP> Op,
                                                          typename Teuchos::ScalarTraits<ScalarType>::magnitudeType kappa,
                                                          typename Teuchos::ScalarTraits<ScalarType>::magnitudeType eps,
                                                          typename Teuchos::ScalarTraits<ScalarType>::magnitudeType tol ) : 
    MatOrthoManager<ScalarType,MV,OP>(Op), kappa_(kappa), eps_(eps), tol_(tol)
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
    , timerReortho_(Teuchos::TimeMonitor::getNewTimer("Anasazi::BasicOrthoManager::Re-orthogonalization"))
#endif
  {
    TEUCHOS_TEST_FOR_EXCEPTION(eps_ < SCT::magnitude(SCT::zero()),std::invalid_argument,
        "Anasazi::BasicOrthoManager::BasicOrthoManager(): argument \"eps\" must be non-negative.");
    if (eps_ == 0) {
      Teuchos::LAPACK<int,MagnitudeType> lapack;
      eps_ = lapack.LAMCH('E');
      eps_ = Teuchos::ScalarTraits<MagnitudeType>::pow(eps_,.75);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
        tol_ < SCT::magnitude(SCT::zero()) || tol_ > SCT::magnitude(SCT::one()),
        std::invalid_argument,
        "Anasazi::BasicOrthoManager::BasicOrthoManager(): argument \"tol\" must be in [0,1].");
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthonormality
  template<class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
  BasicOrthoManager<ScalarType,MV,OP>::orthonormErrorMat(const MV &X, Teuchos::RCP<const MV> MX) const {
    const ScalarType ONE = SCT::one();
    int rank = MVT::GetNumberVecs(X);
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(rank,rank);
    MatOrthoManager<ScalarType,MV,OP>::innerProdMat(X,X,xTx,MX,MX);
    for (int i=0; i<rank; i++) {
      xTx(i,i) -= ONE;
    }
    return xTx.normFrobenius();
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthogonality
  template<class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
  BasicOrthoManager<ScalarType,MV,OP>::orthogErrorMat(const MV &X1, const MV &X2, Teuchos::RCP<const MV> MX1, Teuchos::RCP<const MV> MX2) const {
    int r1 = MVT::GetNumberVecs(X1);
    int r2  = MVT::GetNumberVecs(X2);
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(r1,r2);
    MatOrthoManager<ScalarType,MV,OP>::innerProdMat(X1,X2,xTx,MX1,MX2);
    return xTx.normFrobenius();
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  template<class ScalarType, class MV, class OP>
  void BasicOrthoManager<ScalarType, MV, OP>::projectMat(
          MV &X, 
          Teuchos::Array<Teuchos::RCP<const MV> > Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
          Teuchos::RCP<MV> MX,
          Teuchos::Array<Teuchos::RCP<const MV> > MQ
      ) const {
    // For the inner product defined by the operator Op or the identity (Op == 0)
    //   -> Orthogonalize X against each Q[i]
    // Modify MX accordingly
    //
    // Note that when Op is 0, MX is not referenced
    //
    // Parameter variables
    //
    // X  : Vectors to be transformed
    //
    // MX : Image of the block vector X by the mass matrix
    //
    // Q  : Bases to orthogonalize against. These are assumed orthonormal, mutually and independently.
    //

#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    // Get a FancyOStream from out_arg or create a new one ...
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    *out << "Entering Anasazi::BasicOrthoManager::projectMat(...)\n";
#endif

    ScalarType ONE  = SCT::one();

    int xc = MVT::GetNumberVecs( X );
    ptrdiff_t xr = MVText::GetGlobalLength( X );
    int nq = Q.length();
    std::vector<int> qcs(nq);
    // short-circuit
    if (nq == 0 || xc == 0 || xr == 0) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
      *out << "Leaving Anasazi::BasicOrthoManager::projectMat(...)\n";
#endif
      return;
    }
    ptrdiff_t qr = MVText::GetGlobalLength ( *Q[0] );
    // if we don't have enough C, expand it with null references
    // if we have too many, resize to throw away the latter ones
    // if we have exactly as many as we have Q, this call has no effect
    C.resize(nq);


    /******   DO NO MODIFY *MX IF _hasOp == false   ******/
    if (this->_hasOp) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
      *out << "Allocating MX...\n";
#endif
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,MVT::GetNumberVecs(X));
        OPT::Apply(*(this->_Op),X,*MX);
        this->_OpCounter += MVT::GetNumberVecs(X);
      }
    }
    else {
      // Op == I  -->  MX = X (ignore it if the user passed it in)
      MX = Teuchos::rcpFromRef(X);
    }
    int mxc = MVT::GetNumberVecs( *MX );
    ptrdiff_t mxr = MVText::GetGlobalLength( *MX );

    // check size of X and Q w.r.t. common sense
    TEUCHOS_TEST_FOR_EXCEPTION( xc<0 || xr<0 || mxc<0 || mxr<0, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::projectMat(): MVT returned negative dimensions for X,MX" );
    // check size of X w.r.t. MX and Q
    TEUCHOS_TEST_FOR_EXCEPTION( xc!=mxc || xr!=mxr || xr!=qr, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::projectMat(): Size of X not consistent with MX,Q" );

    // tally up size of all Q and check/allocate C
    int baslen = 0;
    for (int i=0; i<nq; i++) {
      TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength( *Q[i] ) != qr, std::invalid_argument, 
                          "Anasazi::BasicOrthoManager::projectMat(): Q lengths not mutually consistent" );
      qcs[i] = MVT::GetNumberVecs( *Q[i] );
      TEUCHOS_TEST_FOR_EXCEPTION( qr < static_cast<ptrdiff_t>(qcs[i]), std::invalid_argument, 
                          "Anasazi::BasicOrthoManager::projectMat(): Q has less rows than columns" );
      baslen += qcs[i];

      // check size of C[i]
      if ( C[i] == Teuchos::null ) {
        C[i] = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(qcs[i],xc) );
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION( C[i]->numRows() != qcs[i] || C[i]->numCols() != xc , std::invalid_argument, 
                           "Anasazi::BasicOrthoManager::projectMat(): Size of Q not consistent with size of C" );
      }
    }

    // Perform the Gram-Schmidt transformation for a block of vectors

    // Compute the initial Op-norms
    std::vector<ScalarType> oldDot( xc );
    MVT::MvDot( X, *MX, oldDot );
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    *out << "oldDot = { ";
    std::copy(oldDot.begin(), oldDot.end(), std::ostream_iterator<ScalarType>(*out, " "));
    *out << "}\n";
#endif

    MQ.resize(nq);
    // Define the product Q^T * (Op*X)
    for (int i=0; i<nq; i++) {
      // Multiply Q' with MX
      MatOrthoManager<ScalarType,MV,OP>::innerProdMat(*Q[i],X,*C[i],MQ[i],MX);
      // Multiply by Q and subtract the result in X
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
      *out << "Applying projector P_Q[" << i << "]...\n";
#endif
      MVT::MvTimesMatAddMv( -ONE, *Q[i], *C[i], ONE, X );

      // Update MX, with the least number of applications of Op as possible
      // Update MX. If we have MQ, use it. Otherwise, just multiply by Op
      if (this->_hasOp) {
        if (MQ[i] == Teuchos::null) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
          *out << "Updating MX via M*X...\n";
#endif
          OPT::Apply( *(this->_Op), X, *MX);
          this->_OpCounter += MVT::GetNumberVecs(X);
        }
        else {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
          *out << "Updating MX via M*Q...\n";
#endif
          MVT::MvTimesMatAddMv( -ONE, *MQ[i], *C[i], ONE, *MX );
        }
      }
    }

    // Compute new Op-norms
    std::vector<ScalarType> newDot(xc);
    MVT::MvDot( X, *MX, newDot );
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    *out << "newDot = { ";
    std::copy(newDot.begin(), newDot.end(), std::ostream_iterator<ScalarType>(*out, " "));
    *out << "}\n";
#endif

    // determine (individually) whether to do another step of classical Gram-Schmidt
    for (int j = 0; j < xc; ++j) {
      
      if ( SCT::magnitude(kappa_*newDot[j]) < SCT::magnitude(oldDot[j]) ) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
        *out << "kappa_*newDot[" <<j<< "] == " << kappa_*newDot[j] << "... another step of Gram-Schmidt.\n";
#endif
#ifdef ANASAZI_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor lcltimer( *timerReortho_ );
#endif
        for (int i=0; i<nq; i++) {
          Teuchos::SerialDenseMatrix<int,ScalarType> C2(*C[i]);
          
          // Apply another step of classical Gram-Schmidt
          MatOrthoManager<ScalarType,MV,OP>::innerProdMat(*Q[i],X,C2,MQ[i],MX);
          *C[i] += C2;
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
          *out << "Applying projector P_Q[" << i << "]...\n";
#endif
          MVT::MvTimesMatAddMv( -ONE, *Q[i], C2, ONE, X );
          
          // Update MX as above
          if (this->_hasOp) {
            if (MQ[i] == Teuchos::null) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
              *out << "Updating MX via M*X...\n";
#endif
              OPT::Apply( *(this->_Op), X, *MX);
              this->_OpCounter += MVT::GetNumberVecs(X);
            }
            else {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
              *out << "Updating MX via M*Q...\n";
#endif
              MVT::MvTimesMatAddMv( -ONE, *MQ[i], C2, ONE, *MX );
            }
          }
        }
        break;
      } // if (kappa_*newDot[j] < oldDot[j])
    } // for (int j = 0; j < xc; ++j)

#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    *out << "Leaving Anasazi::BasicOrthoManager::projectMat(...)\n";
#endif
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X), with rank numvectors(X)
  template<class ScalarType, class MV, class OP>
  int BasicOrthoManager<ScalarType, MV, OP>::normalizeMat(
        MV &X, 
        Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
        Teuchos::RCP<MV> MX) const {
    // call findBasis(), with the instruction to try to generate a basis of rank numvecs(X)
    // findBasis() requires MX

    int xc = MVT::GetNumberVecs(X);
    ptrdiff_t xr = MVText::GetGlobalLength(X);

    // if Op==null, MX == X (via pointer)
    // Otherwise, either the user passed in MX or we will allocated and compute it
    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
        MX = MVT::Clone(X,xc);
        OPT::Apply(*(this->_Op),X,*MX);
        this->_OpCounter += MVT::GetNumberVecs(X);
      }
    }

    // if the user doesn't want to store the coefficients, 
    // allocate some local memory for them 
    if ( B == Teuchos::null ) {
      B = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xc,xc) );
    }

    int mxc = (this->_hasOp) ? MVT::GetNumberVecs( *MX ) : xc;
    ptrdiff_t mxr = (this->_hasOp) ? MVText::GetGlobalLength( *MX )  : xr;

    // check size of C, B
    TEUCHOS_TEST_FOR_EXCEPTION( xc == 0 || xr == 0, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::normalizeMat(): X must be non-empty" );
    TEUCHOS_TEST_FOR_EXCEPTION( B->numRows() != xc || B->numCols() != xc, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::normalizeMat(): Size of X not consistent with size of B" );
    TEUCHOS_TEST_FOR_EXCEPTION( xc != mxc || xr != mxr, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::normalizeMat(): Size of X not consistent with size of MX" );
    TEUCHOS_TEST_FOR_EXCEPTION( static_cast<ptrdiff_t>(xc) > xr, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::normalizeMat(): Size of X not feasible for normalization" );

    return findBasis(X, MX, *B, true );
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X) - span(W)
  template<class ScalarType, class MV, class OP>
  int BasicOrthoManager<ScalarType, MV, OP>::projectAndNormalizeMat(
          MV &X, 
          Teuchos::Array<Teuchos::RCP<const MV> >  Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
          Teuchos::RCP<MV> MX,
          Teuchos::Array<Teuchos::RCP<const MV> > MQ
      ) const {

#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    // Get a FancyOStream from out_arg or create a new one ...
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    *out << "Entering Anasazi::BasicOrthoManager::projectAndNormalizeMat(...)\n";
#endif

    int nq = Q.length();
    int xc = MVT::GetNumberVecs( X );
    ptrdiff_t xr = MVText::GetGlobalLength( X );
    int rank;

    /* if the user doesn't want to store the coefficients, 
     * allocate some local memory for them 
     */
    if ( B == Teuchos::null ) {
      B = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xc,xc) );
    }

    /******   DO NO MODIFY *MX IF _hasOp == false   ******/
    if (this->_hasOp) {
      if (MX == Teuchos::null) {
        // we need to allocate space for MX
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
        *out << "Allocating MX...\n";
#endif
        MX = MVT::Clone(X,MVT::GetNumberVecs(X));
        OPT::Apply(*(this->_Op),X,*MX);
        this->_OpCounter += MVT::GetNumberVecs(X);
      }
    }
    else {
      // Op == I  -->  MX = X (ignore it if the user passed it in)
      MX = Teuchos::rcpFromRef(X);
    }

    int mxc = MVT::GetNumberVecs( *MX );
    ptrdiff_t mxr = MVText::GetGlobalLength( *MX );

    TEUCHOS_TEST_FOR_EXCEPTION( xc == 0 || xr == 0, std::invalid_argument, "Anasazi::BasicOrthoManager::projectAndNormalizeMat(): X must be non-empty" );

    ptrdiff_t numbas = 0;
    for (int i=0; i<nq; i++) {
      numbas += MVT::GetNumberVecs( *Q[i] );
    }

    // check size of B
    TEUCHOS_TEST_FOR_EXCEPTION( B->numRows() != xc || B->numCols() != xc, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::projectAndNormalizeMat(): Size of X must be consistent with size of B" );
    // check size of X and MX
    TEUCHOS_TEST_FOR_EXCEPTION( xc<0 || xr<0 || mxc<0 || mxr<0, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::projectAndNormalizeMat(): MVT returned negative dimensions for X,MX" );
    // check size of X w.r.t. MX 
    TEUCHOS_TEST_FOR_EXCEPTION( xc!=mxc || xr!=mxr, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::projectAndNormalizeMat(): Size of X must be consistent with size of MX" );
    // check feasibility
    TEUCHOS_TEST_FOR_EXCEPTION( numbas+xc > xr, std::invalid_argument, 
                        "Anasazi::BasicOrthoManager::projectAndNormalizeMat(): Orthogonality constraints not feasible" );

    // orthogonalize all of X against Q
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    *out << "Orthogonalizing X against Q...\n";
#endif
    projectMat(X,Q,C,MX,MQ);

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
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
      *out << "Attempting to find orthonormal basis for X...\n";
#endif
      rank = findBasis(X,MX,*B,false,curxsize);

      if (oldrank != -1 && rank != oldrank) {
        // we had previously stopped before, after operating on vector oldrank
        // we saved its coefficients, augmented it with a random vector, and
        // then called findBasis() again, which proceeded to add vector oldrank
        // to the basis. 
        // now, restore the saved coefficients into B
        for (int i=0; i<xc; i++) {
          (*B)(i,oldrank) = oldCoeff(i,0);
        }
      }

      if (rank < xc) {
        if (rank != oldrank) {
          // we quit on this vector and will augment it with random below
          // this is the first time that we have quit on this vector
          // therefor, (*B)(:,rank) contains the actual coefficients of the 
          // input vectors with respect to the previous vectors in the basis
          // save these values, as (*B)(:,rank) will be overwritten by our next
          // call to findBasis()
          // we will restore it after we are done working on this vector
          for (int i=0; i<xc; i++) {
            oldCoeff(i,0) = (*B)(i,rank);
          }
        }
      }

      if (rank == xc) {
        // we are done
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
        *out << "Finished computing basis.\n";
#endif
        break;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION( rank < oldrank, OrthoError,   
                            "Anasazi::BasicOrthoManager::projectAndNormalizeMat(): basis lost rank; this shouldn't happen");

        if (rank != oldrank) {
          // we added a vector to the basis; reset the chance counter
          numTries = 10;
          // store old rank
          oldrank = rank;
        }
        else {
          // has this vector run out of chances to escape degeneracy?
          if (numTries <= 0) {
            break;
          }
        }
        // use one of this vector's chances
        numTries--;

        // randomize troubled direction
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
        *out << "Randomizing X[" << rank << "]...\n";
#endif
        Teuchos::RCP<MV> curX, curMX;
        std::vector<int> ind(1);
        ind[0] = rank;
        curX  = MVT::CloneViewNonConst(X,ind);
        MVT::MvRandom(*curX);
        if (this->_hasOp) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
          *out << "Applying operator to random vector.\n";
#endif
          curMX = MVT::CloneViewNonConst(*MX,ind);
          OPT::Apply( *(this->_Op), *curX, *curMX );
          this->_OpCounter += MVT::GetNumberVecs(*curX);
        }

        // orthogonalize against Q
        // if !this->_hasOp, the curMX will be ignored.
        // we don't care about these coefficients
        // on the contrary, we need to preserve the previous coeffs
        {
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummyC(0);
          projectMat(*curX,Q,dummyC,curMX,MQ);
        }
      }
    } while (1);

    // this should never raise an exception; but our post-conditions oblige us to check
    TEUCHOS_TEST_FOR_EXCEPTION( rank > xc || rank < 0, std::logic_error, 
                        "Anasazi::BasicOrthoManager::projectAndNormalizeMat(): Debug error in rank variable." );

#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    *out << "Leaving Anasazi::BasicOrthoManager::projectAndNormalizeMat(...)\n";
#endif

    return rank;
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X), with the option of extending the subspace so that 
  // the rank is numvectors(X)
  template<class ScalarType, class MV, class OP>
  int BasicOrthoManager<ScalarType, MV, OP>::findBasis(
                MV &X, Teuchos::RCP<MV> MX, 
                Teuchos::SerialDenseMatrix<int,ScalarType> &B,
                bool completeBasis, int howMany ) const {

    // For the inner product defined by the operator Op or the identity (Op == 0)
    //   -> Orthonormalize X 
    // Modify MX accordingly
    //
    // Note that when Op is 0, MX is not referenced
    //
    // Parameter variables
    //
    // X  : Vectors to be orthonormalized
    //
    // MX : Image of the multivector X under the operator Op
    //
    // Op  : Pointer to the operator for the inner product
    //

#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    // Get a FancyOStream from out_arg or create a new one ...
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    *out << "Entering Anasazi::BasicOrthoManager::findBasis(...)\n";
#endif

    const ScalarType ONE  = SCT::one();
    const MagnitudeType ZERO = SCT::magnitude(SCT::zero());

    int xc = MVT::GetNumberVecs( X );

    if (howMany == -1) {
      howMany = xc;
    }

    /*******************************************************
     *  If _hasOp == false, we will not reference MX below *
     *******************************************************/
    TEUCHOS_TEST_FOR_EXCEPTION(this->_hasOp == true && MX == Teuchos::null, std::logic_error,
        "Anasazi::BasicOrthoManager::findBasis(): calling routine did not specify MS.");
    TEUCHOS_TEST_FOR_EXCEPTION( howMany < 0 || howMany > xc, std::logic_error, 
                        "Anasazi::BasicOrthoManager::findBasis(): Invalid howMany parameter" );

    /* xstart is which column we are starting the process with, based on howMany
     * columns before xstart are assumed to be Op-orthonormal already
     */
    int xstart = xc - howMany;

    for (int j = xstart; j < xc; j++) {

      // numX represents the number of currently orthonormal columns of X
      int numX = j;
      // j represents the index of the current column of X
      // these are different interpretations of the same value

      // 
      // set the lower triangular part of B to zero
      for (int i=j+1; i<xc; ++i) {
        B(i,j) = ZERO;
      }

      // Get a view of the vector currently being worked on.
      std::vector<int> index(1);
      index[0] = j;
      Teuchos::RCP<MV> Xj = MVT::CloneViewNonConst( X, index );
      Teuchos::RCP<MV> MXj;
      if ((this->_hasOp)) {
        // MXj is a view of the current vector in MX
        MXj = MVT::CloneViewNonConst( *MX, index );
      }
      else {
        // MXj is a pointer to Xj, and MUST NOT be modified
        MXj = Xj;
      }

      // Get a view of the previous vectors.
      std::vector<int> prev_idx( numX );
      Teuchos::RCP<const MV> prevX, prevMX;

      if (numX > 0) {
        for (int i=0; i<numX; ++i) prev_idx[i] = i;
        prevX = MVT::CloneViewNonConst( X, prev_idx );
        if (this->_hasOp) {
          prevMX = MVT::CloneViewNonConst( *MX, prev_idx );
        }
      } 

      bool rankDef = true;
      /* numTrials>0 will denote that the current vector was randomized for the purpose
       * of finding a basis vector, and that the coefficients of that vector should
       * not be stored in B
       */
      for (int numTrials = 0; numTrials < 10; numTrials++) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
        *out << "Trial " << numTrials << " for vector " << j << "\n";
#endif

        // Make storage for these Gram-Schmidt iterations.
        Teuchos::SerialDenseMatrix<int,ScalarType> product(numX, 1);
        std::vector<MagnitudeType> origNorm(1), newNorm(1), newNorm2(1);

        //
        // Save old MXj vector and compute Op-norm
        //
        Teuchos::RCP<MV> oldMXj = MVT::CloneCopy( *MXj ); 
        MatOrthoManager<ScalarType,MV,OP>::normMat(*Xj,origNorm,MXj);
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
        *out << "origNorm = " << origNorm[0] << "\n";
#endif

        if (numX > 0) {
          // Apply the first step of Gram-Schmidt

          // product <- prevX^T MXj
          MatOrthoManager<ScalarType,MV,OP>::innerProdMat(*prevX,*Xj,product,Teuchos::null,MXj);

          // Xj <- Xj - prevX prevX^T MXj   
          //     = Xj - prevX product
#ifdef ANASAZI_BASIC_ORTHO_DEBUG 
          *out << "Orthogonalizing X[" << j << "]...\n";
#endif
          MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *Xj );

          // Update MXj
          if (this->_hasOp) {
            // MXj <- Op*Xj_new
            //      = Op*(Xj_old - prevX prevX^T MXj)
            //      = MXj - prevMX product
#ifdef ANASAZI_BASIC_ORTHO_DEBUG 
            *out << "Updating MX[" << j << "]...\n";
#endif
            MVT::MvTimesMatAddMv( -ONE, *prevMX, product, ONE, *MXj );
          }

          // Compute new Op-norm
          MatOrthoManager<ScalarType,MV,OP>::normMat(*Xj,newNorm,MXj);
          MagnitudeType product_norm = product.normOne();
          
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
          *out << "newNorm = " << newNorm[0] << "\n";
          *out << "prodoct_norm = " << product_norm << "\n";
#endif

          // Check if a correction is needed.
          if ( product_norm/newNorm[0] >= tol_ || newNorm[0] < eps_*origNorm[0]) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
            if (product_norm/newNorm[0] >= tol_) {
              *out << "product_norm/newNorm == " << product_norm/newNorm[0] << "... another step of Gram-Schmidt.\n";
            }
            else {
              *out << "eps*origNorm == " << eps_*origNorm[0] << "... another step of Gram-Schmidt.\n";
            }
#endif
            // Apply the second step of Gram-Schmidt
            // This is the same as above
            Teuchos::SerialDenseMatrix<int,ScalarType> P2(numX,1);
            MatOrthoManager<ScalarType,MV,OP>::innerProdMat(*prevX,*Xj,P2,Teuchos::null,MXj);
            product += P2;
#ifdef ANASAZI_BASIC_ORTHO_DEBUG 
            *out << "Orthogonalizing X[" << j << "]...\n";
#endif
            MVT::MvTimesMatAddMv( -ONE, *prevX, P2, ONE, *Xj );
            if ((this->_hasOp)) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG 
              *out << "Updating MX[" << j << "]...\n";
#endif
              MVT::MvTimesMatAddMv( -ONE, *prevMX, P2, ONE, *MXj );
            }
            // Compute new Op-norms
            MatOrthoManager<ScalarType,MV,OP>::normMat(*Xj,newNorm2,MXj);
            product_norm = P2.normOne();
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
            *out << "newNorm2 = " << newNorm2[0] << "\n";
            *out << "product_norm = " << product_norm << "\n";
#endif
            if ( product_norm/newNorm2[0] >= tol_ || newNorm2[0] < eps_*origNorm[0] ) {
              // we don't do another GS, we just set it to zero.
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
              if (product_norm/newNorm2[0] >= tol_) {
                *out << "product_norm/newNorm2 == " << product_norm/newNorm2[0] << "... setting vector to zero.\n";
              }
              else if (newNorm[0] < newNorm2[0]) {
                *out << "newNorm2 > newNorm... setting vector to zero.\n";
              }
              else {
                *out << "eps*origNorm == " << eps_*origNorm[0] << "... setting vector to zero.\n";
              }
#endif
              MVT::MvInit(*Xj,ZERO);
              if ((this->_hasOp)) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
                *out << "Setting MX[" << j << "] to zero as well.\n";
#endif
                MVT::MvInit(*MXj,ZERO);
              }
            }
          } 
        } // if (numX > 0) do GS

        // save the coefficients, if we are working on the original vector and not a randomly generated one
        if (numTrials == 0) {
          for (int i=0; i<numX; i++) {
            B(i,j) = product(i,0);
          }
        }

        // Check if Xj has any directional information left after the orthogonalization.
        MatOrthoManager<ScalarType,MV,OP>::normMat(*Xj,newNorm,MXj);
        if ( newNorm[0] != ZERO && newNorm[0] > SCT::sfmin() ) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
          *out << "Normalizing X[" << j << "], norm(X[" << j << "]) = " << newNorm[0] << "\n";
#endif
          // Normalize Xj.
          // Xj <- Xj / norm
          MVT::MvScale( *Xj, ONE/newNorm[0]);
          if (this->_hasOp) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
            *out << "Normalizing M*X[" << j << "]...\n";
#endif
            // Update MXj.
            MVT::MvScale( *MXj, ONE/newNorm[0]);
          }

          // save it, if it corresponds to the original vector and not a randomly generated one
          if (numTrials == 0) {
            B(j,j) = newNorm[0];
          }

          // We are not rank deficient in this vector. Move on to the next vector in X.
          rankDef = false;
          break;
        }
        else {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
          *out << "Not normalizing M*X[" << j << "]...\n";
#endif
          // There was nothing left in Xj after orthogonalizing against previous columns in X.
          // X is rank deficient.
          // reflect this in the coefficients
          B(j,j) = ZERO;

          if (completeBasis) {
            // Fill it with random information and keep going.
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
            *out << "Inserting random vector in X[" << j << "]...\n";
#endif
            MVT::MvRandom( *Xj );
            if (this->_hasOp) {
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
              *out << "Updating M*X[" << j << "]...\n";
#endif
              OPT::Apply( *(this->_Op), *Xj, *MXj );
              this->_OpCounter += MVT::GetNumberVecs(*Xj);
            }
          }
          else {
            rankDef = true;
            break;
          }
        }
      }  // for (numTrials = 0; numTrials < 10; ++numTrials)

      // if rankDef == true, then quit and notify user of rank obtained
      if (rankDef == true) {
        TEUCHOS_TEST_FOR_EXCEPTION( rankDef && completeBasis, OrthoError, 
                            "Anasazi::BasicOrthoManager::findBasis(): Unable to complete basis" );
#ifdef ANASAZI_BASIC_ORTHO_DEBUG
        *out << "Returning early, rank " << j << " from Anasazi::BasicOrthoManager::findBasis(...)\n";
#endif
        return j;
      }

    } // for (j = 0; j < xc; ++j)

#ifdef ANASAZI_BASIC_ORTHO_DEBUG
    *out << "Returning " << xc << " from Anasazi::BasicOrthoManager::findBasis(...)\n";
#endif
    return xc;
  }

} // namespace Anasazi

#endif // ANASAZI_BASIC_ORTHOMANAGER_HPP

