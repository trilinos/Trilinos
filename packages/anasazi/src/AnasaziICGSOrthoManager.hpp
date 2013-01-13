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


/*! \file AnasaziICGSOrthoManager.hpp
  \brief Basic implementation of the Anasazi::OrthoManager class
*/

#ifndef ANASAZI_ICSG_ORTHOMANAGER_HPP
#define ANASAZI_ICSG_ORTHOMANAGER_HPP

/*!   \class Anasazi::ICGSOrthoManager
      \brief An implementation of the Anasazi::GenOrthoManager that performs orthogonalization
      using iterated classical Gram-Schmidt.
      
      \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziGenOrthoManager.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#ifdef ANASAZI_ICGS_DEBUG
#  include <Teuchos_FancyOStream.hpp>
#endif

namespace Anasazi {

  template<class ScalarType, class MV, class OP>
  class ICGSOrthoManager : public GenOrthoManager<ScalarType,MV,OP> {

  private:
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<ScalarType>  SCT;
    typedef MultiVecTraits<ScalarType,MV>      MVT;
    typedef MultiVecTraitsExt<ScalarType,MV>   MVText;
    typedef OperatorTraits<ScalarType,MV,OP>   OPT;

  public:

    //! @name Constructor/Destructor
    //@{ 
    //! Constructor specifying the operator defining the inner product as well as the number of orthogonalization iterations.
    ICGSOrthoManager( Teuchos::RCP<const OP> Op = Teuchos::null, int numIters = 2,
                      typename Teuchos::ScalarTraits<ScalarType>::magnitudeType eps = 0.0,
                      typename Teuchos::ScalarTraits<ScalarType>::magnitudeType tol = 0.20 );


    //! Destructor
    ~ICGSOrthoManager() {}
    //@}


    //! @name Methods implementing Anasazi::GenOrthoManager
    //@{ 

    /*! \brief Applies a series of generic projectors.
     *
     * Given a list of bases <tt>X[i]</tt> and <tt>Y[i]</tt> (a projection pair), this method
     * takes a multivector \c S and applies the projectors
     * \f[
     *    P_{X[i],Y[i]} S = S - X[i] \langle Y[i], X[i] \rangle^{-1} \langle Y[i], S \rangle\ .
     * \f]
     * This operation projects \c S onto the space orthogonal to the <tt>Y[i]</tt>, 
     * along the range of the <tt>X[i]</tt>. The inner product specified by \f$\langle \cdot, 
     * \cdot \rangle\f$ is given by innerProd(). 
     *
     * \note The call 
     * \code 
     * projectGen(S, tuple(X1,X2), tuple(Y1,Y2)) 
     * \endcode 
     * is equivalent to the call
     * \code 
     * projectGen(S, tuple(X2,X1), tuple(Y2,Y1))
     * \endcode
     *
     * The method also returns the coefficients <tt>C[i]</tt> associated with each projection pair, so that 
     * \f[
     *    S_{in} = S_{out} + \sum_i X[i] C[i]
     * \f]
     * and therefore
     * \f[
     *    C[i] = \langle Y[i], X[i] \rangle^{-1} \langle Y[i], S \rangle\ .
     * \f]
     *
     * Lastly, for reasons of efficiency, the user must specify whether the projection pairs are bi-orthonormal with
     * respect to innerProd(), i.e., whether \f$\langle Y[i], X[i] \rangle = I\f$. In the case that the bases are specified
     * to be biorthogonal, the inverse \f$\langle Y, X \rangle^{-1}\f$ will not be computed. Furthermore, the user may optionally 
     * specifiy the image of \c S and the projection pairs under the inner product operator getOp().
     *
     * projectGen() is implemented to apply the projectors via an iterated Classical Gram-Schmidt, where the iteration is performed
     * getNumIters() number of times.
     *
     @param S [in/out] The multivector to be modified.<br>
       On output, the columns of \c S will be orthogonal to each <tt>Y[i]</tt>, satisfying
       \f[
          \langle Y[i], S_{out} \rangle = 0
       \f]
       Also, 
       \f[
          S_{in} = S_{out} + \sum_i X[i] C[i]
       \f]

     @param X [in] Multivectors for bases under which \f$S_{in}\f$ is modified.
        
     @param Y [in] Multivectors for bases to which \f$S_{out}\f$ should be orthogonal.

     @param isBiortho [in] A flag specifying whether the bases <tt>X[i]</tt>
            and <tt>Y[i]</tt> are biorthonormal, i.e,. whether \f$\langle Y[i],
            X[i]\rangle == I\f$.

     @param C [out] Coefficients for reconstructing \f$S_{in}\f$ via the bases <tt>X[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c S and <tt>X[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>.<br>
       If <tt>C[i]</tt> points to a Teuchos::SerialDenseMatrix with size
       inconsistent with \c S and \c <tt>X[i]</tt>, then a std::invalid_argument
       exception will be thrown.<br>
       Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null pointer,
       the caller will not have access to the computed coefficients <tt>C[i]</tt>.

     @param MS [in/out] If specified by the user, on input \c MS is required to be the image of \c S under the operator getOp(). 
     On output, \c MS will be updated to reflect the changes in \c S.
     
     @param MX [in] If specified by the user, on <tt>MX[i]</tt> is required to be the image of <tt>X[i]</tt> under the operator getOp().
     @param MY [in] If specified by the user, on <tt>MY[i]</tt> is required to be the image of <tt>Y[i]</tt> under the operator getOp().

     \pre 
     <ul>
     <li>If <tt>X[i] != Teuchos::null</tt> or <tt>Y[i] != Teuchos::null</tt>, then <tt>X[i]</tt> and <tt>Y[i]</tt> are required to 
         have the same number of columns, and each should have the same number of rows as \c S.
     <li>For any <tt>i != j</tt>, \f$\langle Y[i], X[j] \rangle == 0\f$.
     <li>If <tt>biOrtho == true</tt>, \f$\langle Y[i], X[i]\rangle == I\f$
     <li>Otherwise, if <tt>biOrtho == false</tt>, then \f$\langle Y[i], X[i]\rangle\f$ should be Hermitian positive-definite.
     <li>If <tt>X[i]</tt> and <tt>Y[i]</tt> have \f$xc_i\f$ columns and \c S has \f$sc\f$ columns, then <tt>C[i]</tt> if specified must
         be \f$xc_i \times sc\f$.
     </ul>
     */
    void projectGen( 
          MV &S,
          Teuchos::Array<Teuchos::RCP<const MV> > X,
          Teuchos::Array<Teuchos::RCP<const MV> > Y,
          bool isBiOrtho,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
          Teuchos::RCP<MV> MS                        = Teuchos::null,
          Teuchos::Array<Teuchos::RCP<const MV> > MX = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null)),
          Teuchos::Array<Teuchos::RCP<const MV> > MY = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null))
        ) const;


    /*! \brief Applies a series of generic projectors and returns an orthonormal basis for the residual data.
     *
     * Given a list of bases <tt>X[i]</tt> and <tt>Y[i]</tt> (a projection pair), this method
     * takes a multivector \c S and applies the projectors
     * \f[
     *    P_{X[i],Y[i]} S = S - X[i] \langle Y[i], X[i] \rangle^{-1} \langle Y[i], S \rangle\ .
     * \f]
     * These operation project \c S onto the space orthogonal to the range of the <tt>Y[i]</tt>, 
     * along the range of \c X[i]. The inner product specified by \f$\langle \cdot, \cdot \rangle\f$
     * is given by innerProd(). 
     *
     * The method returns in \c S an orthonormal basis for the residual
     * \f[
     *    \left( \prod_{i} P_{X[i],Y[i]} \right) S_{in} = S_{out} B\ ,
     * \f]
     * where \c B contains the (not necessarily triangular) coefficients of the residual with respect to the 
     * new basis.
     *
     * The method also returns the coefficients <tt>C[i]</tt> and \c B associated with each projection pair, so that 
     * \f[
     *    S_{in} = S_{out} B + \sum_i X[i] C[i]
     * \f]
     * and 
     * \f[
     *    C[i] = \langle Y[i], X[i] \rangle^{-1} \langle Y[i], S \rangle\ .
     * \f]
     *
     * Lastly, for reasons of efficiency, the user must specify whether the projection pairs are bi-orthonormal with
     * respect to innerProd(), i.e., whether \f$\langle Y[i], X[i] \rangle = I\f$. Furthermore, the user may optionally 
     * specifiy the image of \c S and the projection pairs under the inner product operator getOp().
     @param S [in/out] The multivector to be modified.<br>
       On output, the columns of \c S will be orthogonal to each <tt>Y[i]</tt>, satisfying
       \f[
          \langle Y[i], S_{out} \rangle = 0
       \f]
       Also, 
       \f[
          S_{in}(1:m,1:n) = S_{out}(1:m,1:rank) B(1:rank,1:n) + \sum_i X[i] C[i]\ ,
       \f]
       where \c m is the number of rows in \c S, \c n is the number of
       columns in \c S, and \c rank is the value returned from the method.

     @param X [in] Multivectors for bases under which \f$S_{in}\f$ is modified.
        
     @param Y [in] Multivectors for bases to which \f$S_{out}\f$ should be orthogonal.

     @param isBiortho [in] A flag specifying whether the bases <tt>X[i]</tt>
            and <tt>Y[i]</tt> are biorthonormal, i.e,. whether \f$\langle Y[i],
            X[i]\rangle == I\f$.

     @param C [out] Coefficients for reconstructing \f$S_{in}\f$ via the bases <tt>X[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c X and <tt>Q[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>.<br>
       If <tt>C[i]</tt> points to a Teuchos::SerialDenseMatrix with size
       inconsistent with \c S and \c <tt>X[i]</tt>, then a std::invalid_argument
       exception will be thrown.<br>
       Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null pointer,
       the caller will not have access to the computed coefficients <tt>C[i]</tt>.

     @param B [out] The coefficients of the original \c S with respect to the computed basis. If \c B is a non-null pointer and 
     \c B matches the dimensions of \c B, then the
     coefficients computed during the orthogonalization routine will be stored in \c B, similar to calling 
       \code
          innerProd( Sout, Sin, B );
       \endcode
     If \c B points to a Teuchos::SerialDenseMatrix with size inconsistent with
     \c S, then a std::invalid_argument exception will be thrown.<br>
     Otherwise, if \c B is null, the caller will not have access to the computed
     coefficients.<br>
     The normalization uses classical Gram-Schmidt iteration, so that \c B is an upper triangular matrix with positive diagonal elements.

     @param MS [in/out] If specified by the user, on input \c MS is required to be the image of \c S under the operator getOp(). 
     On output, \c MS will be updated to reflect the changes in \c S.
     
     @param MX [in] If specified by the user, on <tt>MX[i]</tt> is required to be the image of <tt>X[i]</tt> under the operator getOp().
     @param MY [in] If specified by the user, on <tt>MY[i]</tt> is required to be the image of <tt>Y[i]</tt> under the operator getOp().

     \pre 
     <ul>
     <li>If <tt>X[i] != Teuchos::null</tt> or <tt>Y[i] != Teuchos::null</tt>, then <tt>X[i]</tt> and <tt>Y[i]</tt> are required to 
         have the same number of columns, and each should have the same number of rows as \c S.
     <li>For any <tt>i != j</tt>, \f$\langle Y[i], X[j] \rangle == 0\f$.
     <li>If <tt>biOrtho == true</tt>, \f$\langle Y[i], X[i]\rangle == I\f$
     <li>Otherwise, if <tt>biOrtho == false</tt>, then \f$\langle Y[i], X[i]\rangle\f$ should be Hermitian positive-definite.
     <li>If <tt>X[i]</tt> and <tt>Y[i]</tt> have \f$xc_i\f$ columns and \c S has \f$sc\f$ columns, then <tt>C[i]</tt> if specified must
         be \f$xc_i \times sc\f$.
     <li>If <tt>S</tt> has \f$sc\f$ columns, then \c B if specified must be \f$sc \times sc \f$.
     </ul>

     @return Rank of the basis computed by this method.
     */
    int projectAndNormalizeGen (
          MV &S,
          Teuchos::Array<Teuchos::RCP<const MV> > X,
          Teuchos::Array<Teuchos::RCP<const MV> > Y,
          bool isBiOrtho,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null,
          Teuchos::RCP<MV> MS                        = Teuchos::null,
          Teuchos::Array<Teuchos::RCP<const MV> > MX = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null)),
          Teuchos::Array<Teuchos::RCP<const MV> > MY = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null))
        ) const;


    //@}


    //! @name Methods implementing Anasazi::MatOrthoManager
    //@{ 


    /*! \brief Given a list of mutually orthogonal and internally orthonormal bases \c Q, this method
     * projects a multivector \c X onto the space orthogonal to the individual <tt>Q[i]</tt>, 
     * optionally returning the coefficients of \c X for the individual <tt>Q[i]</tt>. All of this is done with respect
     * to the inner product innerProd().
     * 
     * This method calls projectGen() as follows:
     * \code
     *    projectGen(X,Q,Q,true,C,MX,MQ,MQ);
     * \endcode
     * See projectGen() for argument requirements.
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
     * This method calls projectAndNormalizeGen() as follows:
     * \code
     *    projectAndNormalizeGen(X,empty,empty,true,empty,B,MX);
     * \endcode
     * See projectAndNormalizeGen() for argument requirements.
    */
    int normalizeMat ( 
          MV &X, 
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null,
          Teuchos::RCP<MV> MX                                         = Teuchos::null
        ) const;


    /*! \brief Given a set of bases <tt>Q[i]</tt> and a multivector \c X, this method computes an orthonormal basis for \f$colspan(X) - \sum_i colspan(Q[i])\f$.
     *
     * This method calls projectAndNormalizeGen() as follows:
     * \code
     *    projectAndNormalizeGen(X,Q,Q,true,C,B,MX,MQ,MQ);
     * \endcode
     * See projectAndNormalizeGen() for argument requirements.
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
    void setNumIters( int numIters ) {
      numIters_ = numIters;
      TEUCHOS_TEST_FOR_EXCEPTION(numIters_ < 1,std::invalid_argument,
          "Anasazi::ICGSOrthoManager::setNumIters(): input must be >= 1.");
    }

    //! Return parameter for re-orthogonalization threshold.
    int getNumIters() const { return numIters_; } 

    //@} 

  private:
    MagnitudeType eps_;
    MagnitudeType tol_;
    
    //! Parameter for re-orthogonalization.
    int numIters_;
  
    // ! Routine to find an orthonormal basis
    int findBasis(MV &X, Teuchos::RCP<MV> MX, 
                         Teuchos::SerialDenseMatrix<int,ScalarType> &B, 
                         bool completeBasis, int howMany = -1) const;
  };



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template<class ScalarType, class MV, class OP>
  ICGSOrthoManager<ScalarType,MV,OP>::ICGSOrthoManager( Teuchos::RCP<const OP> Op,
                                                        int numIters,
                                                        typename Teuchos::ScalarTraits<ScalarType>::magnitudeType eps,
                                                        typename Teuchos::ScalarTraits<ScalarType>::magnitudeType tol) : 
    GenOrthoManager<ScalarType,MV,OP>(Op), eps_(eps), tol_(tol)
  {
    setNumIters(numIters);
    TEUCHOS_TEST_FOR_EXCEPTION(eps_ < SCT::magnitude(SCT::zero()),std::invalid_argument,
        "Anasazi::ICGSOrthoManager::ICGSOrthoManager(): argument \"eps\" must be non-negative.");
    if (eps_ == 0) {
      Teuchos::LAPACK<int,MagnitudeType> lapack;
      eps_ = lapack.LAMCH('E');
      eps_ = Teuchos::ScalarTraits<MagnitudeType>::pow(eps_,.50);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
        tol_ < SCT::magnitude(SCT::zero()) || tol_ > SCT::magnitude(SCT::one()),
        std::invalid_argument,
        "Anasazi::ICGSOrthoManager::ICGSOrthoManager(): argument \"tol\" must be in [0,1].");
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute the distance from orthonormality
  template<class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
  ICGSOrthoManager<ScalarType,MV,OP>::orthonormErrorMat(const MV &X, Teuchos::RCP<const MV> MX) const {
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
  ICGSOrthoManager<ScalarType,MV,OP>::orthogErrorMat(const MV &X1, const MV &X2, Teuchos::RCP<const MV> MX1, Teuchos::RCP<const MV> MX2) const {
    int r1 = MVT::GetNumberVecs(X1);
    int r2  = MVT::GetNumberVecs(X2);
    Teuchos::SerialDenseMatrix<int,ScalarType> xTx(r1,r2);
    MatOrthoManager<ScalarType,MV,OP>::innerProdMat(X1,X2,xTx,MX1,MX2);
    return xTx.normFrobenius();
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  template<class ScalarType, class MV, class OP>
  void ICGSOrthoManager<ScalarType, MV, OP>::projectMat(
          MV &X, 
          Teuchos::Array<Teuchos::RCP<const MV> > Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
          Teuchos::RCP<MV> MX,
          Teuchos::Array<Teuchos::RCP<const MV> > MQ
      ) const 
  {
    projectGen(X,Q,Q,true,C,MX,MQ,MQ);
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X), with rank numvectors(X)
  template<class ScalarType, class MV, class OP>
  int ICGSOrthoManager<ScalarType, MV, OP>::normalizeMat(
        MV &X, 
        Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
        Teuchos::RCP<MV> MX) const 
  {
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
                        "Anasazi::ICGSOrthoManager::normalizeMat(): X must be non-empty" );
    TEUCHOS_TEST_FOR_EXCEPTION( B->numRows() != xc || B->numCols() != xc, std::invalid_argument, 
                        "Anasazi::ICGSOrthoManager::normalizeMat(): Size of X not consistent with size of B" );
    TEUCHOS_TEST_FOR_EXCEPTION( xc != mxc || xr != mxr, std::invalid_argument, 
                        "Anasazi::ICGSOrthoManager::normalizeMat(): Size of X not consistent with size of MX" );
    TEUCHOS_TEST_FOR_EXCEPTION( static_cast<ptrdiff_t>(xc) > xr, std::invalid_argument, 
                        "Anasazi::ICGSOrthoManager::normalizeMat(): Size of X not feasible for normalization" );

    return findBasis(X, MX, *B, true );
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X) - span(W)
  template<class ScalarType, class MV, class OP>
  int ICGSOrthoManager<ScalarType, MV, OP>::projectAndNormalizeMat(
          MV &X, 
          Teuchos::Array<Teuchos::RCP<const MV> >  Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
          Teuchos::RCP<MV> MX,
          Teuchos::Array<Teuchos::RCP<const MV> > MQ
      ) const 
  {
    return projectAndNormalizeGen(X,Q,Q,true,C,B,MX,MQ,MQ);
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  template<class ScalarType, class MV, class OP>
  void ICGSOrthoManager<ScalarType, MV, OP>::projectGen(
          MV &S,
          Teuchos::Array<Teuchos::RCP<const MV> > X,
          Teuchos::Array<Teuchos::RCP<const MV> > Y,
          bool isBiortho,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
          Teuchos::RCP<MV> MS,
          Teuchos::Array<Teuchos::RCP<const MV> > MX,
          Teuchos::Array<Teuchos::RCP<const MV> > MY
      ) const 
  {
    // For the inner product defined by the operator Op or the identity (Op == 0)
    //   -> Orthogonalize S against each Y[i], modifying it in the range of X[i]
    // Modify MS accordingly
    //
    // Note that when Op is 0, MS is not referenced
    //
    // Parameter variables
    //
    // S  : Multivector to be transformed
    //
    // MS : Image of the block vector S by the mass matrix
    //
    // X,Y: Bases to orthogonalize against/via.
    //

#ifdef ANASAZI_ICGS_DEBUG
    // Get a FancyOStream from out_arg or create a new one ...
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    *out << "Entering Anasazi::ICGSOrthoManager::projectGen(...)\n";
#endif

    const ScalarType     ONE = SCT::one();
    const MagnitudeType ZERO = SCT::magnitude(SCT::zero());
    Teuchos::LAPACK<int,ScalarType> lapack;
    Teuchos::BLAS<int,ScalarType> blas;

    int sc = MVT::GetNumberVecs( S );
    ptrdiff_t sr = MVText::GetGlobalLength( S );
    int numxy = X.length();
    TEUCHOS_TEST_FOR_EXCEPTION(X.length() != Y.length(),std::invalid_argument, 
        "Anasazi::ICGSOrthoManager::projectGen(): X and Y must contain the same number of multivectors.");
    std::vector<int> xyc(numxy);
    // short-circuit
    if (numxy == 0 || sc == 0 || sr == 0) {
#ifdef ANASAZI_ICGS_DEBUG
      *out << "Leaving Anasazi::ICGSOrthoManager::projectGen(...)\n";
#endif
      return;
    }
    // if we don't have enough C, expand it with null references
    // if we have too many, resize to throw away the latter ones
    // if we have exactly as many as we have X,Y this call has no effect
    //
    // same for MX, MY
    C.resize(numxy);
    MX.resize(numxy);
    MY.resize(numxy);

    // check size of S w.r.t. common sense
    TEUCHOS_TEST_FOR_EXCEPTION( sc<0 || sr<0, std::invalid_argument, 
                        "Anasazi::ICGSOrthoManager::projectGen(): MVT returned negative dimensions for S." );

    // check size of MS
    if (this->_hasOp == true) {
      if (MS != Teuchos::null) {
        TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*MS) != sr, std::invalid_argument, 
            "Anasazi::ICGSOrthoManager::projectGen(): MS length not consistent with S." );
        TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*MS) != sc, std::invalid_argument,
            "Anasazi::ICGSOrthoManager::projectGen(): MS width not consistent with S." );
      }
    }

    // tally up size of all X,Y and check/allocate C
    ptrdiff_t sumxyc = 0;
    int MYmissing = 0;
    int MXmissing = 0;
    for (int i=0; i<numxy; i++) {
      if (X[i] != Teuchos::null && Y[i] != Teuchos::null) {
        TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*X[i]) != sr, std::invalid_argument, 
            "Anasazi::ICGSOrthoManager::projectGen(): X[" << i << "] length not consistent with S." );
        TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*Y[i]) != sr, std::invalid_argument, 
            "Anasazi::ICGSOrthoManager::projectGen(): Y[" << i << "] length not consistent with S." );
        TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*X[i]) != MVT::GetNumberVecs(*Y[i]), std::invalid_argument,
            "Anasazi::ICGSOrthoManager::projectGen(): X[" << i << "] and Y[" << i << "] widths not consistent." );

        xyc[i] = MVT::GetNumberVecs( *X[i] );
        TEUCHOS_TEST_FOR_EXCEPTION( sr < static_cast<ptrdiff_t>(xyc[i]), std::invalid_argument, 
            "Anasazi::ICGSOrthoManager::projectGen(): X[" << i << "],Y[" << i << "] have less rows than columns, and therefore cannot be full rank." );
        sumxyc += xyc[i];

        // check size of C[i]
        if ( C[i] == Teuchos::null ) {
          C[i] = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xyc[i],sc) );
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION( C[i]->numRows() != xyc[i] || C[i]->numCols() != sc , std::invalid_argument, 
              "Anasazi::ICGSOrthoManager::projectGen(): Size of Q not consistent with size of C." );
        }
        // check sizes of MX[i], MY[i] if present
        // if not present, count their absence
        if (MX[i] != Teuchos::null) {
          TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*MX[i]) != sr || MVT::GetNumberVecs(*MX[i]) != xyc[i], std::invalid_argument,
              "Anasazi::ICGSOrthoManager::projectGen(): Size of MX[" << i << "] not consistent with size of X[" << i << "]." );
        }
        else {
          MXmissing += xyc[i];
        }
        if (MY[i] != Teuchos::null) {
          TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*MY[i]) != sr || MVT::GetNumberVecs(*MY[i]) != xyc[i], std::invalid_argument,
              "Anasazi::ICGSOrthoManager::projectGen(): Size of MY[" << i << "] not consistent with size of Y[" << i << "]." );
        }
        else {
          MYmissing += xyc[i];
        }
      }
      else {
        // if one is null and the other is not... the user may have made a mistake
        TEUCHOS_TEST_FOR_EXCEPTION(X[i] != Teuchos::null || Y[i] != Teuchos::null, std::invalid_argument,
            "Anasazi::ICGSOrthoManager::projectGen(): " 
            << (X[i] == Teuchos::null ? "Y[" : "X[") << i << "] was provided but " 
            << (X[i] == Teuchos::null ? "X[" : "Y[") << i << "] was not.");
      }
    }

    // is this operation even feasible?
    TEUCHOS_TEST_FOR_EXCEPTION(sumxyc > sr, std::invalid_argument, 
        "Anasazi::ICGSOrthoManager::projectGen(): dimension of all X[i],Y[i] is "
        << sumxyc << ", but length of vectors is only " << sr << ". This is infeasible.");


    /******   DO NO MODIFY *MS IF _hasOp == false   
     * if _hasOp == false, we don't need MS, MX or MY
     *
     * if _hasOp == true, we need MS (for S M-norms) and
     * therefore, we must also update MS, regardless of whether
     * it gets returned to the user (i.e., whether the user provided it)
     * we may need to allocate and compute MX or MY
     *
     * let BXY denote:
     *    "X" - we have all M*X[i]
     *    "Y" - we have all M*Y[i]
     *    "B" - we have biorthogonality (for all X[i],Y[i])
     * Encode these as values 
     *   X = 1
     *   Y = 2
     *   B = 4
     * We have 8 possibilities, 0-7
     *
     * We must allocate storage and computer the following, lest 
     * innerProdMat do it for us:
     *  none (0) - allocate MX, for inv(<X,Y>) and M*S
     *
     * for the following, we will compute M*S manually before returning
     *   B(4)  BY(6)  Y(2)                    -->  updateMS = 1
     * for the following, we will update M*S as we go, using M*X
     *   XY(3)  X(1)  none(0)  BXY(7)  BX(5)  -->  updateMS = 2
     * these choices favor applications of M over allocation of memory
     *
     */
    int updateMS = -1;
    if (this->_hasOp) {
      int whichAlloc = 0;
      if (MXmissing == 0) {
        whichAlloc += 1;
      }
      if (MYmissing == 0) {
        whichAlloc += 2;
      }
      if (isBiortho) {
        whichAlloc += 4;
      }

      switch (whichAlloc) {
        case 2:
        case 4:
        case 6:
          updateMS = 1;
          break;
        case 0:
        case 1:
        case 3:
        case 5:
        case 7:
          updateMS = 2;
          break;
      }

      // produce MS
      if (MS == Teuchos::null) {
#ifdef ANASAZI_ICGS_DEBUG
    *out << "Allocating MS...\n";
#endif
        MS = MVT::Clone(S,MVT::GetNumberVecs(S));
        OPT::Apply(*(this->_Op),S,*MS);
        this->_OpCounter += MVT::GetNumberVecs(S);
      }

      // allocate the rest
      if (whichAlloc == 0) {
        // allocate and compute missing MX
        for (int i=0; i<numxy; i++) {
          if (MX[i] == Teuchos::null) {
#ifdef ANASAZI_ICGS_DEBUG
            *out << "Allocating MX[" << i << "]...\n";
#endif
            Teuchos::RCP<MV> tmpMX = MVT::Clone(*X[i],xyc[i]);
            OPT::Apply(*(this->_Op),*X[i],*tmpMX);
            MX[i] = tmpMX;
            this->_OpCounter += xyc[i];
          }
        }
      }
    }
    else {
      // Op == I  -->  MS == S 
      MS = Teuchos::rcpFromRef(S);
      updateMS = 0;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(updateMS == -1,std::logic_error,
        "Anasazi::ICGSOrthoManager::projectGen(): Error in updateMS logic.");


    ////////////////////////////////////////////////////////////////////
    // Perform the Gram-Schmidt transformation for a block of vectors
    ////////////////////////////////////////////////////////////////////

    // Compute Cholesky factorizations for the Y'*M*X
    // YMX stores the YMX (initially) and their Cholesky factorizations (utlimately)
    Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > YMX(numxy);
    if (isBiortho == false) {
      for (int i=0; i<numxy; i++) {
#ifdef ANASAZI_ICGS_DEBUG
        *out << "Computing YMX[" << i << "] and its Cholesky factorization...\n";
#endif
        YMX[i] = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xyc[i],xyc[i]) );
        MatOrthoManager<ScalarType,MV,OP>::innerProdMat(*Y[i],*X[i],*YMX[i],MY[i],MX[i]);
#ifdef ANASAZI_ICGS_DEBUG
        // YMX should be symmetric positive definite
        // the cholesky will check the positive definiteness, but it looks only at the upper half
        // we will check the symmetry by removing the upper half from the lower half, which should 
        // result in zeros
        // also, diagonal of YMX should be real; consider the complex part as error
        {
          MagnitudeType err = ZERO;
          for (int jj=0; jj<YMX[i]->numCols(); jj++) {
            err =+ SCT::magnitude(SCT::imag((*YMX[i])(jj,jj)));
            for (int ii=jj; ii<YMX[i]->numRows(); ii++) {
              err += SCT::magnitude( (*YMX[i])(ii,jj) - SCT::conjugate((*YMX[i])(jj,ii)) );
            }
          }
          *out << "Symmetry error in YMX[" << i << "] == " << err << "\n";
        }
#endif
        // take the LU factorization
        int info;
        lapack.POTRF('U',YMX[i]->numRows(),YMX[i]->values(),YMX[i]->stride(),&info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error,
            "Anasazi::ICGSOrthoManager::projectGen(): Error computing Cholesky factorization of Y[i]^T * M * X[i] using POTRF: returned info " << info);
      }
    }

    // Compute the initial Op-norms
#ifdef ANASAZI_ICGS_DEBUG
    std::vector<MagnitudeType> oldNorms(sc);
    MatOrthoManager<ScalarType,MV,OP>::normMat(S,oldNorms,MS);
    *out << "oldNorms = { ";
    std::copy(oldNorms.begin(), oldNorms.end(), std::ostream_iterator<MagnitudeType>(*out, " "));
    *out << "}\n";
#endif


    // clear the C[i] and allocate Ccur
    Teuchos::Array<Teuchos::SerialDenseMatrix<int,ScalarType> > Ccur(numxy);
    for (int i=0; i<numxy; i++) {
      C[i]->putScalar(ZERO);
      Ccur[i].reshape(C[i]->numRows(),C[i]->numCols());
    }

    for (int iter=0; iter < numIters_; iter++) {
#ifdef ANASAZI_ICGS_DEBUG
      *out << "beginning iteration " << iter+1 << "\n";
#endif

      // Define the product Y[i]'*Op*S
      for (int i=0; i<numxy; i++) {
        // Compute Y[i]'*M*S
        MatOrthoManager<ScalarType,MV,OP>::innerProdMat(*Y[i],S,Ccur[i],MY[i],MS);
        if (isBiortho == false) {
          // C[i] = inv(YMX[i])*(Y[i]'*M*S)
          int info;
          lapack.POTRS('U',YMX[i]->numCols(),Ccur[i].numCols(),
              YMX[i]->values(),YMX[i]->stride(),
              Ccur[i].values(),Ccur[i].stride(), &info);
          TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::logic_error,
              "Anasazi::ICGSOrthoManager::projectGen(): Error code " << info << " from lapack::POTRS." );
        }

        // Multiply by X[i] and subtract the result in S
#ifdef ANASAZI_ICGS_DEBUG
        *out << "Applying projector P_{X[" << i << "],Y[" << i << "]}...\n";
#endif
        MVT::MvTimesMatAddMv( -ONE, *X[i], Ccur[i], ONE, S );

        // Accumulate coeffs into previous step
        *C[i] += Ccur[i];

        // Update MS as required
        if (updateMS == 1) {
#ifdef ANASAZI_ICGS_DEBUG
          *out << "Updating MS...\n";
#endif
          OPT::Apply( *(this->_Op), S, *MS);
          this->_OpCounter += sc;
        }
        else if (updateMS == 2) {
#ifdef ANASAZI_ICGS_DEBUG
          *out << "Updating MS...\n";
#endif
          MVT::MvTimesMatAddMv( -ONE, *MX[i], Ccur[i], ONE, *MS );
        }
      }

      // Compute new Op-norms
#ifdef ANASAZI_ICGS_DEBUG
      std::vector<MagnitudeType> newNorms(sc);
      MatOrthoManager<ScalarType,MV,OP>::normMat(S,newNorms,MS);
      *out << "newNorms = { ";
      std::copy(newNorms.begin(), newNorms.end(), std::ostream_iterator<MagnitudeType>(*out, " "));
      *out << "}\n";
#endif
    }

#ifdef ANASAZI_ICGS_DEBUG
    *out << "Leaving Anasazi::ICGSOrthoManager::projectGen(...)\n";
#endif
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(S) - span(Y)
  template<class ScalarType, class MV, class OP>
  int ICGSOrthoManager<ScalarType, MV, OP>::projectAndNormalizeGen(
          MV &S,
          Teuchos::Array<Teuchos::RCP<const MV> > X,
          Teuchos::Array<Teuchos::RCP<const MV> > Y,
          bool isBiortho,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B,
          Teuchos::RCP<MV> MS,
          Teuchos::Array<Teuchos::RCP<const MV> > MX,
          Teuchos::Array<Teuchos::RCP<const MV> > MY
      ) const {
    // For the inner product defined by the operator Op or the identity (Op == 0)
    //   -> Orthogonalize S against each Y[i], modifying it in the range of X[i]
    // Modify MS accordingly
    // Then construct a M-orthonormal basis for the remaining part
    //
    // Note that when Op is 0, MS is not referenced
    //
    // Parameter variables
    //
    // S  : Multivector to be transformed
    //
    // MS : Image of the block vector S by the mass matrix
    //
    // X,Y: Bases to orthogonalize against/via.
    //

#ifdef ANASAZI_ICGS_DEBUG
    // Get a FancyOStream from out_arg or create a new one ...
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    *out << "Entering Anasazi::ICGSOrthoManager::projectAndNormalizeGen(...)\n";
#endif

    int rank;
    int sc = MVT::GetNumberVecs( S );
    ptrdiff_t sr = MVText::GetGlobalLength( S );
    int numxy = X.length();
    TEUCHOS_TEST_FOR_EXCEPTION(X.length() != Y.length(),std::invalid_argument, 
        "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): X and Y must contain the same number of multivectors.");
    std::vector<int> xyc(numxy);
    // short-circuit
    if (sc == 0 || sr == 0) {
#ifdef ANASAZI_ICGS_DEBUG
      *out << "Leaving Anasazi::ICGSOrthoManager::projectGen(...)\n";
#endif
      return 0;
    }
    // if we don't have enough C, expand it with null references
    // if we have too many, resize to throw away the latter ones
    // if we have exactly as many as we have X,Y this call has no effect
    //
    // same for MX, MY
    C.resize(numxy);
    MX.resize(numxy);
    MY.resize(numxy);

    // check size of S w.r.t. common sense
    TEUCHOS_TEST_FOR_EXCEPTION( sc<0 || sr<0, std::invalid_argument, 
                        "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): MVT returned negative dimensions for S." );

    // check size of MS
    if (this->_hasOp == true) {
      if (MS != Teuchos::null) {
        TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*MS) != sr, std::invalid_argument, 
            "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): MS length not consistent with S." );
        TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*MS) != sc, std::invalid_argument,
            "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): MS width not consistent with S." );
      }
    }

    // tally up size of all X,Y and check/allocate C
    ptrdiff_t sumxyc = 0;
    int MYmissing = 0;
    int MXmissing = 0;
    for (int i=0; i<numxy; i++) {
      if (X[i] != Teuchos::null && Y[i] != Teuchos::null) {
        TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*X[i]) != sr, std::invalid_argument, 
            "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): X[" << i << "] length not consistent with S." );
        TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*Y[i]) != sr, std::invalid_argument, 
            "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): Y[" << i << "] length not consistent with S." );
        TEUCHOS_TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*X[i]) != MVT::GetNumberVecs(*Y[i]), std::invalid_argument,
            "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): X[" << i << "] and Y[" << i << "] widths not consistent." );

        xyc[i] = MVT::GetNumberVecs( *X[i] );
        TEUCHOS_TEST_FOR_EXCEPTION( sr < static_cast<ptrdiff_t>(xyc[i]), std::invalid_argument, 
            "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): X[" << i << "],Y[" << i << "] have less rows than columns, and therefore cannot be full rank." );
        sumxyc += xyc[i];

        // check size of C[i]
        if ( C[i] == Teuchos::null ) {
          C[i] = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(xyc[i],sc) );
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION( C[i]->numRows() != xyc[i] || C[i]->numCols() != sc , std::invalid_argument, 
              "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): Size of Q not consistent with size of C." );
        }
        // check sizes of MX[i], MY[i] if present
        // if not present, count their absence
        if (MX[i] != Teuchos::null) {
          TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*MX[i]) != sr || MVT::GetNumberVecs(*MX[i]) != xyc[i], std::invalid_argument,
              "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): Size of MX[" << i << "] not consistent with size of X[" << i << "]." );
        }
        else {
          MXmissing += xyc[i];
        }
        if (MY[i] != Teuchos::null) {
          TEUCHOS_TEST_FOR_EXCEPTION( MVText::GetGlobalLength(*MY[i]) != sr || MVT::GetNumberVecs(*MY[i]) != xyc[i], std::invalid_argument,
              "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): Size of MY[" << i << "] not consistent with size of Y[" << i << "]." );
        }
        else {
          MYmissing += xyc[i];
        }
      }
      else {
        // if one is null and the other is not... the user may have made a mistake
        TEUCHOS_TEST_FOR_EXCEPTION(X[i] != Teuchos::null || Y[i] != Teuchos::null, std::invalid_argument,
            "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): " 
            << (X[i] == Teuchos::null ? "Y[" : "X[") << i << "] was provided but " 
            << (X[i] == Teuchos::null ? "X[" : "Y[") << i << "] was not.");
      }
    }

    // is this operation even feasible?
    TEUCHOS_TEST_FOR_EXCEPTION(sumxyc + sc > sr, std::invalid_argument, 
        "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): dimension of all X[i],Y[i] is "
        << sumxyc << " and requested " << sc << "-dimensional basis, but length of vectors is only " 
        << sr << ". This is infeasible.");


    /******   DO NO MODIFY *MS IF _hasOp == false   
     * if _hasOp == false, we don't need MS, MX or MY
     *
     * if _hasOp == true, we need MS (for S M-norms and normalization) and
     * therefore, we must also update MS, regardless of whether
     * it gets returned to the user (i.e., whether the user provided it)
     * we may need to allocate and compute MX or MY
     *
     * let BXY denote:
     *    "X" - we have all M*X[i]
     *    "Y" - we have all M*Y[i]
     *    "B" - we have biorthogonality (for all X[i],Y[i])
     * Encode these as values 
     *   X = 1
     *   Y = 2
     *   B = 4
     * We have 8 possibilities, 0-7
     *
     * We must allocate storage and computer the following, lest 
     * innerProdMat do it for us:
     *  none (0) - allocate MX, for inv(<X,Y>) and M*S
     *
     * for the following, we will compute M*S manually before returning
     *   B(4)  BY(6)  Y(2)                    -->  updateMS = 1
     * for the following, we will update M*S as we go, using M*X
     *   XY(3)  X(1)  none(0)  BXY(7)  BX(5)  -->  updateMS = 2
     * these choices favor applications of M over allocation of memory
     *
     */
    int updateMS = -1;
    if (this->_hasOp) {
      int whichAlloc = 0;
      if (MXmissing == 0) {
        whichAlloc += 1;
      }
      if (MYmissing == 0) {
        whichAlloc += 2;
      }
      if (isBiortho) {
        whichAlloc += 4;
      }

      switch (whichAlloc) {
        case 2:
        case 4:
        case 6:
          updateMS = 1;
          break;
        case 0:
        case 1:
        case 3:
        case 5:
        case 7:
          updateMS = 2;
          break;
      }

      // produce MS
      if (MS == Teuchos::null) {
#ifdef ANASAZI_ICGS_DEBUG
    *out << "Allocating MS...\n";
#endif
        MS = MVT::Clone(S,MVT::GetNumberVecs(S));
        OPT::Apply(*(this->_Op),S,*MS);
        this->_OpCounter += MVT::GetNumberVecs(S);
      }

      // allocate the rest
      if (whichAlloc == 0) {
        // allocate and compute missing MX
        for (int i=0; i<numxy; i++) {
          if (MX[i] == Teuchos::null) {
#ifdef ANASAZI_ICGS_DEBUG
            *out << "Allocating MX[" << i << "]...\n";
#endif
            Teuchos::RCP<MV> tmpMX = MVT::Clone(*X[i],xyc[i]);
            OPT::Apply(*(this->_Op),*X[i],*tmpMX);
            MX[i] = tmpMX;
            this->_OpCounter += xyc[i];
          }
        }
      }
    }
    else {
      // Op == I  -->  MS == S 
      MS = Teuchos::rcpFromRef(S);
      updateMS = 0;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(updateMS == -1,std::logic_error,
        "Anasazi::ICGSOrthoManager::projectGen(): Error in updateMS logic.");


    // if the user doesn't want to store the coefficients, 
    // allocate some local memory for them 
    if ( B == Teuchos::null ) {
      B = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(sc,sc) );
    }
    else {
      // check size of B
      TEUCHOS_TEST_FOR_EXCEPTION( B->numRows() != sc || B->numCols() != sc, std::invalid_argument, 
          "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): Size of S must be consistent with size of B" );
    }


    // orthogonalize all of S against X,Y
    projectGen(S,X,Y,isBiortho,C,MS,MX,MY);

    Teuchos::SerialDenseMatrix<int,ScalarType> oldCoeff(sc,1);
    // start working
    rank = 0;
    int numTries = 10;   // each vector in X gets 10 random chances to escape degeneracy
    int oldrank = -1;
    do {
      int curssize = sc - rank;

      // orthonormalize S, but quit if it is rank deficient
      // we can't let findBasis generated random vectors to complete the basis,
      // because it doesn't know about Q; we will do this ourselves below
#ifdef ANASAZI_ICGS_DEBUG
      *out << "Attempting to find orthonormal basis for X...\n";
#endif
      rank = findBasis(S,MS,*B,false,curssize);

      if (oldrank != -1 && rank != oldrank) {
        // we had previously stopped before, after operating on vector oldrank
        // we saved its coefficients, augmented it with a random vector, and
        // then called findBasis() again, which proceeded to add vector oldrank
        // to the basis. 
        // now, restore the saved coefficients into B
        for (int i=0; i<sc; i++) {
          (*B)(i,oldrank) = oldCoeff(i,0);
        }
      }

      if (rank < sc) {
        if (rank != oldrank) {
          // we quit on this vector and will augment it with random below
          // this is the first time that we have quit on this vector
          // therefor, (*B)(:,rank) contains the actual coefficients of the 
          // input vectors with respect to the previous vectors in the basis
          // save these values, as (*B)(:,rank) will be overwritten by our next
          // call to findBasis()
          // we will restore it after we are done working on this vector
          for (int i=0; i<sc; i++) {
            oldCoeff(i,0) = (*B)(i,rank);
          }
        }
      }

      if (rank == sc) {
        // we are done
#ifdef ANASAZI_ICGS_DEBUG
        *out << "Finished computing basis.\n";
#endif
        break;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION( rank < oldrank, OrthoError,   
                            "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): basis lost rank; this shouldn't happen");

        if (rank != oldrank) {
          // we added a basis vector from random info; reset the chance counter
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
#ifdef ANASAZI_ICGS_DEBUG
        *out << "Inserting random vector in X[" << rank << "]. Attempt " << 10-numTries << ".\n";
#endif
        Teuchos::RCP<MV> curS, curMS;
        std::vector<int> ind(1);
        ind[0] = rank;
        curS  = MVT::CloneViewNonConst(S,ind);
        MVT::MvRandom(*curS);
        if (this->_hasOp) {
#ifdef ANASAZI_ICGS_DEBUG
          *out << "Applying operator to random vector.\n";
#endif
          curMS = MVT::CloneViewNonConst(*MS,ind);
          OPT::Apply( *(this->_Op), *curS, *curMS );
          this->_OpCounter += MVT::GetNumberVecs(*curS);
        }

        // orthogonalize against X,Y
        // if !this->_hasOp, the curMS will be ignored.
        // we don't care about these coefficients
        // on the contrary, we need to preserve the previous coeffs
        {
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummyC(0);
          projectGen(*curS,X,Y,isBiortho,dummyC,curMS,MX,MY);
        }
      }
    } while (1);

    // this should never raise an exception; but our post-conditions oblige us to check
    TEUCHOS_TEST_FOR_EXCEPTION( rank > sc || rank < 0, std::logic_error, 
                        "Anasazi::ICGSOrthoManager::projectAndNormalizeGen(): Debug error in rank variable." );

#ifdef ANASAZI_ICGS_DEBUG
    *out << "Returning " << rank << " from Anasazi::ICGSOrthoManager::projectAndNormalizeGen(...)\n";
#endif

    return rank;
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Find an Op-orthonormal basis for span(X), with the option of extending the subspace so that 
  // the rank is numvectors(X)
  template<class ScalarType, class MV, class OP>
  int ICGSOrthoManager<ScalarType, MV, OP>::findBasis(
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

#ifdef ANASAZI_ICGS_DEBUG
    // Get a FancyOStream from out_arg or create a new one ...
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    *out << "Entering Anasazi::ICGSOrthoManager::findBasis(...)\n";
#endif

    const ScalarType     ONE = SCT::one();
    const MagnitudeType ZERO = SCT::magnitude(SCT::zero());

    int xc = MVT::GetNumberVecs( X );

    if (howMany == -1) {
      howMany = xc;
    }

    /*******************************************************
     *  If _hasOp == false, we will not reference MX below *
     *******************************************************/
    TEUCHOS_TEST_FOR_EXCEPTION(this->_hasOp == true && MX == Teuchos::null, std::logic_error,
        "Anasazi::ICGSOrthoManager::findBasis(): calling routine did not specify MS.");
    TEUCHOS_TEST_FOR_EXCEPTION( howMany < 0 || howMany > xc, std::logic_error, 
                        "Anasazi::ICGSOrthoManager::findBasis(): Invalid howMany parameter" );

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
        prevX = MVT::CloneView( X, prev_idx );
        if (this->_hasOp) {
          prevMX = MVT::CloneView( *MX, prev_idx );
        }
      } 

      bool rankDef = true;
      /* numTrials>0 will denote that the current vector was randomized for the purpose
       * of finding a basis vector, and that the coefficients of that vector should
       * not be stored in B
       */
      for (int numTrials = 0; numTrials < 10; numTrials++) {
#ifdef ANASAZI_ICGS_DEBUG
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
#ifdef ANASAZI_ICGS_DEBUG
        *out << "origNorm = " << origNorm[0] << "\n";
#endif

        if (numX > 0) {
          // Apply the first step of Gram-Schmidt

          // product <- prevX^T MXj
          MatOrthoManager<ScalarType,MV,OP>::innerProdMat(*prevX,*Xj,product,Teuchos::null,MXj);

          // Xj <- Xj - prevX prevX^T MXj   
          //     = Xj - prevX product
#ifdef ANASAZI_ICGS_DEBUG
          *out << "Orthogonalizing X[" << j << "]...\n";
#endif
          MVT::MvTimesMatAddMv( -ONE, *prevX, product, ONE, *Xj );

          // Update MXj
          if (this->_hasOp) {
            // MXj <- Op*Xj_new
            //      = Op*(Xj_old - prevX prevX^T MXj)
            //      = MXj - prevMX product
#ifdef ANASAZI_ICGS_DEBUG
            *out << "Updating MX[" << j << "]...\n";
#endif
            MVT::MvTimesMatAddMv( -ONE, *prevMX, product, ONE, *MXj );
          }

          // Compute new Op-norm
          MatOrthoManager<ScalarType,MV,OP>::normMat(*Xj,newNorm,MXj);
          MagnitudeType product_norm = product.normOne();
          
#ifdef ANASAZI_ICGS_DEBUG
          *out << "newNorm = " << newNorm[0] << "\n";
          *out << "prodoct_norm = " << product_norm << "\n";
#endif

          // Check if a correction is needed.
          if ( product_norm/newNorm[0] >= tol_ || newNorm[0] < eps_*origNorm[0]) {
#ifdef ANASAZI_ICGS_DEBUG
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
#ifdef ANASAZI_ICGS_DEBUG
            *out << "Orthogonalizing X[" << j << "]...\n";
#endif
            MVT::MvTimesMatAddMv( -ONE, *prevX, P2, ONE, *Xj );
            if ((this->_hasOp)) {
#ifdef ANASAZI_ICGS_DEBUG
              *out << "Updating MX[" << j << "]...\n";
#endif
              MVT::MvTimesMatAddMv( -ONE, *prevMX, P2, ONE, *MXj );
            }
            // Compute new Op-norms
            MatOrthoManager<ScalarType,MV,OP>::normMat(*Xj,newNorm2,MXj);
            product_norm = P2.normOne();
#ifdef ANASAZI_ICGS_DEBUG
            *out << "newNorm2 = " << newNorm2[0] << "\n";
            *out << "product_norm = " << product_norm << "\n";
#endif
            if ( product_norm/newNorm2[0] >= tol_ || newNorm2[0] < eps_*origNorm[0] ) {
              // we don't do another GS, we just set it to zero.
#ifdef ANASAZI_ICGS_DEBUG
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
#ifdef ANASAZI_ICGS_DEBUG
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
#ifdef ANASAZI_ICGS_DEBUG
          *out << "Normalizing X[" << j << "], norm(X[" << j << "]) = " << newNorm[0] << "\n";
#endif
          // Normalize Xj.
          // Xj <- Xj / norm
          MVT::MvScale( *Xj, ONE/newNorm[0]);
          if (this->_hasOp) {
#ifdef ANASAZI_ICGS_DEBUG
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
#ifdef ANASAZI_ICGS_DEBUG
          *out << "Not normalizing M*X[" << j << "]...\n";
#endif
          // There was nothing left in Xj after orthogonalizing against previous columns in X.
          // X is rank deficient.
          // reflect this in the coefficients
          B(j,j) = ZERO;

          if (completeBasis) {
            // Fill it with random information and keep going.
#ifdef ANASAZI_ICGS_DEBUG
            *out << "Inserting random vector in X[" << j << "]...\n";
#endif
            MVT::MvRandom( *Xj );
            if (this->_hasOp) {
#ifdef ANASAZI_ICGS_DEBUG
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
                            "Anasazi::ICGSOrthoManager::findBasis(): Unable to complete basis" );
#ifdef ANASAZI_ICGS_DEBUG
        *out << "Returning early, rank " << j << " from Anasazi::ICGSOrthoManager::findBasis(...)\n";
#endif
        return j;
      }

    } // for (j = 0; j < xc; ++j)

#ifdef ANASAZI_ICGS_DEBUG
    *out << "Returning " << xc << " from Anasazi::ICGSOrthoManager::findBasis(...)\n";
#endif
    return xc;
  }

} // namespace Anasazi

#endif // ANASAZI_ICSG_ORTHOMANAGER_HPP

