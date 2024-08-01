// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziMatOrthoManager.hpp
  \brief  Templated virtual class for providing orthogonalization/orthonormalization methods with matrix-based
          inner products.
*/

#ifndef ANASAZI_MATORTHOMANAGER_HPP
#define ANASAZI_MATORTHOMANAGER_HPP

/*! \class Anasazi::MatOrthoManager

  \brief Anasazi's templated virtual class for providing routines for
  orthogonalization and orthonormalization of multivectors using matrix-based
  inner products.

  This class extends Anasazi::OrthoManager by providing extra calling arguments
  to orthogonalization routines, to reduce the cost of applying the inner
  product in cases where the user already has the image of target multivectors
  under the inner product matrix.

  A concrete implementation of this class is necessary. The user can create
  their own implementation if those supplied are not suitable for their needs.

  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziOrthoManager.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"

namespace Anasazi {

  template <class ScalarType, class MV, class OP>
  class MatOrthoManager : public OrthoManager<ScalarType,MV> {
  public:
    //! @name Constructor/Destructor
    //@{
    //! Default constructor.
    MatOrthoManager(Teuchos::RCP<const OP> Op = Teuchos::null);

    //! Destructor.
    virtual ~MatOrthoManager() {};
    //@}

    //! @name Accessor routines
    //@{

    //! Set operator used for inner product.
    virtual void setOp( Teuchos::RCP<const OP> Op );

    //! Get operator used for inner product.
    virtual Teuchos::RCP<const OP> getOp() const;

    //! Retrieve operator counter.
    /*! This counter returns the number of applications of the operator specifying the inner
     * product. When the operator is applied to a multivector, the counter is incremented by the
     * number of vectors in the multivector. If the operator is not specified, the counter is never
     * incremented.
     */
    int getOpCounter() const;

    //! Reset the operator counter to zero.
    /*! See getOpCounter() for more details.
     */
    void resetOpCounter();

    //@}

    //! @name Matrix-based Orthogonality Methods
    //@{

    /*! \brief Provides a matrix-based inner product.
     *
     * Provides the inner product
     * \f[
     *    \langle x, y \rangle = x^H M y
     * \f]
     * Optionally allows the provision of \f$ M y\f$ and/or \f$ M x\f$. See OrthoManager::innerProd() for more details.
     *
     */
    void innerProdMat(
          const MV& X, const MV& Y,
          Teuchos::SerialDenseMatrix<int,ScalarType>& Z,
          Teuchos::RCP<const MV> MX = Teuchos::null,
          Teuchos::RCP<const MV> MY = Teuchos::null
        ) const;

    /*! \brief Provides the norm induced by the matrix-based inner product.
     *
     *  Provides the norm:
     *  \f[
     *     \|x\|_M = \sqrt{x^H M y}
     *  \f]
     *  Optionally allows the provision of \f$ M x\f$. See OrthoManager::norm() for more details.
     */
    void normMat(
          const MV& X,
          std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > &normvec,
          Teuchos::RCP<const MV> MX = Teuchos::null
        ) const;

    /*! \brief Provides matrix-based projection method.
     *
     * This method optionally allows the provision of \f$ M X\f$ and/or the \f$ M Q[i]\f$. See OrthoManager::project() for more details.
     @param X, Q, C [in/out] As in OrthoManager::project()

     @param MX [in/out] If specified by the user, on input \c MX is required to be the image of \c X under the operator getOp().
     On output, \c MX will be updated to reflect the changes in \c X.

     @param MQ [in] If specified by the user, on \c MQ[i] is required to be the image of <tt>Q[i]</tt> under the operator getOp().
     */
    virtual void projectMat (
          MV &X,
          Teuchos::Array<Teuchos::RCP<const MV> >  Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
          Teuchos::RCP<MV> MX                                                          = Teuchos::null,
          Teuchos::Array<Teuchos::RCP<const MV> > MQ
              = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null))
        ) const = 0;

    /*! \brief Provides matrix-based orthonormalization method.
     *
     * This method optionally allows the provision of \f$ M X\f$. See orthoManager::normalize() for more details.
     @param X, B [in/out] As in OrthoManager::normalize()

     @param MX [in/out] If specified by the user, on input \c MX is required to be the image of \c X under the operator getOp().
     On output, \c MX will be updated to reflect the changes in \c X.

     @return Rank of the basis computed by this method, less than or equal to
       the number of columns in \c X. This specifies how many columns in the
       returned \c X and \c MX and rows in the returned \c B are valid.
    */
    virtual int normalizeMat (
          MV &X,
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null,
          Teuchos::RCP<MV> MX                                         = Teuchos::null
        ) const = 0;


    /*! \brief Provides matrix-based projection/orthonormalization method.
     *
     * This method optionally allows the provision of \f$ M X\f$ and/or the \f$ M Q[i]\f$. See orthoManager::projectAndNormalize() for more details.
     @param X, Q, C, B [in/out] As in OrthoManager::projectAndNormalize()

     @param MX [in/out] If specified by the user, on input \c MX is required to be the image of \c X under the operator getOp().
     On output, \c MX will be updated to reflect the changes in \c X.

     @param MQ [in] If specified by the user, on \c MQ[i] is required to be the image of <tt>Q[i]</tt> under the operator getOp().

     @return Rank of the basis computed by this method, less than or equal to
       the number of columns in \c X. This specifies how many columns in the
       returned \c X and \c MX and rows in the returned \c B are valid.
    */
    virtual int projectAndNormalizeMat (
          MV &X,
          Teuchos::Array<Teuchos::RCP<const MV> >  Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B                  = Teuchos::null,
          Teuchos::RCP<MV> MX                                                          = Teuchos::null,
          Teuchos::Array<Teuchos::RCP<const MV> > MQ
              = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null))
        ) const = 0;

    /*! \brief This method computes the error in orthonormality of a multivector.
     *
     *  This method optionally allows optionally exploits a caller-provided \c MX.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthonormErrorMat(const MV &X, Teuchos::RCP<const MV> MX = Teuchos::null) const = 0;

    /*! \brief This method computes the error in orthogonality of two multivectors.
     *
     *  This method optionally allows optionally exploits a caller-provided \c MX and/or \c MY.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthogErrorMat(
          const MV &X,
          const MV &Y,
          Teuchos::RCP<const MV> MX = Teuchos::null,
          Teuchos::RCP<const MV> MY = Teuchos::null
        ) const = 0;

    //@}

    //! @name Methods implementing Anasazi::OrthoManager
    //@{

    /*! \brief Implements the interface OrthoManager::innerProd().
     *
     * This method calls
     * \code
     * innerProdMat(X,Y,Z);
     * \endcode
     */
    void innerProd( const MV& X, const MV& Y, Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const;

    /*! \brief Implements the interface OrthoManager::norm().
     *
     * This method calls
     * \code
     * normMat(X,normvec);
     * \endcode
     */
    void norm( const MV& X, std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > &normvec ) const;

    /*! \brief Implements the interface OrthoManager::project().
     *
     * This method calls
     * \code
     * projectMat(X,Q,C);
     * \endcode
     */
    void project (
          MV &X,
          Teuchos::Array<Teuchos::RCP<const MV> > Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null))
        ) const;

    /*! \brief Implements the interface OrthoManager::normalize().
     *
     * This method calls
     * \code
     * normalizeMat(X,B);
     * \endcode
     */
    int normalize ( MV &X, Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null) const;

    /*! \brief Implements the interface OrthoManager::projectAndNormalize().
     *
     * This method calls
     * \code
     * projectAndNormalizeMat(X,Q,C,B);
     * \endcode
     */
    int projectAndNormalize (
          MV &X,
          Teuchos::Array<Teuchos::RCP<const MV> > Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null
        ) const;

    /*! \brief Implements the interface OrthoManager::orthonormError().
     *
     * This method calls
     * \code
     * orthonormErrorMat(X);
     * \endcode
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthonormError(const MV &X) const;

    /*! \brief Implements the interface OrthoManager::orthogError().
     *
     * This method calls
     * \code
     * orthogErrorMat(X1,X2);
     * \endcode
     */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
    orthogError(const MV &X1, const MV &X2) const;

    //@}

  protected:
    Teuchos::RCP<const OP> _Op;
    bool _hasOp;
    mutable int _OpCounter;

  };

  template <class ScalarType, class MV, class OP>
  MatOrthoManager<ScalarType,MV,OP>::MatOrthoManager(Teuchos::RCP<const OP> Op)
      : _Op(Op), _hasOp(Op!=Teuchos::null), _OpCounter(0) {}

  template <class ScalarType, class MV, class OP>
  void MatOrthoManager<ScalarType,MV,OP>::setOp( Teuchos::RCP<const OP> Op )
  {
    _Op = Op;
    _hasOp = (_Op != Teuchos::null);
  }

  template <class ScalarType, class MV, class OP>
  Teuchos::RCP<const OP> MatOrthoManager<ScalarType,MV,OP>::getOp() const
  {
    return _Op;
  }

  template <class ScalarType, class MV, class OP>
  int MatOrthoManager<ScalarType,MV,OP>::getOpCounter() const
  {
    return _OpCounter;
  }

  template <class ScalarType, class MV, class OP>
  void MatOrthoManager<ScalarType,MV,OP>::resetOpCounter()
  {
    _OpCounter = 0;
  }

  template <class ScalarType, class MV, class OP>
  void MatOrthoManager<ScalarType,MV,OP>::innerProd(
      const MV& X, const MV& Y, Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const
  {
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef MultiVecTraits<ScalarType,MV>     MVT;
    typedef OperatorTraits<ScalarType,MV,OP>  OPT;

    Teuchos::RCP<const MV> P,Q;
    Teuchos::RCP<MV> R;

    if (_hasOp) {
      // attempt to minimize the amount of work in applying
      if ( MVT::GetNumberVecs(X) < MVT::GetNumberVecs(Y) ) {
        R = MVT::Clone(X,MVT::GetNumberVecs(X));
        OPT::Apply(*_Op,X,*R);
        _OpCounter += MVT::GetNumberVecs(X);
        P = R;
        Q = Teuchos::rcpFromRef(Y);
      }
      else {
        P = Teuchos::rcpFromRef(X);
        R = MVT::Clone(Y,MVT::GetNumberVecs(Y));
        OPT::Apply(*_Op,Y,*R);
        _OpCounter += MVT::GetNumberVecs(Y);
        Q = R;
      }
    }
    else {
      P = Teuchos::rcpFromRef(X);
      Q = Teuchos::rcpFromRef(Y);
    }

    MVT::MvTransMv(SCT::one(),*P,*Q,Z);
  }

  template <class ScalarType, class MV, class OP>
  void MatOrthoManager<ScalarType,MV,OP>::innerProdMat(
      const MV& X, const MV& Y, Teuchos::SerialDenseMatrix<int,ScalarType>& Z, Teuchos::RCP<const MV> MX, Teuchos::RCP<const MV> MY) const
  {
    (void) MX;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef MultiVecTraits<ScalarType,MV>     MVT;
    // typedef OperatorTraits<ScalarType,MV,OP>  OPT; // unused

    Teuchos::RCP<MV> P,Q;

    if ( MY == Teuchos::null ) {
      innerProd(X,Y,Z);
    }
    else if ( _hasOp ) {
      // the user has done the matrix vector for us
      MVT::MvTransMv(SCT::one(),X,*MY,Z);
    }
    else {
      // there is no matrix vector
      MVT::MvTransMv(SCT::one(),X,Y,Z);
    }
#ifdef TEUCHOS_DEBUG
    for (int j=0; j<Z.numCols(); j++) {
      for (int i=0; i<Z.numRows(); i++) {
        TEUCHOS_TEST_FOR_EXCEPTION(SCT::isnaninf(Z(i,j)), std::logic_error,
            "Anasazi::MatOrthoManager::innerProdMat(): detected NaN/inf.");
      }
    }
#endif
  }

  template <class ScalarType, class MV, class OP>
  void MatOrthoManager<ScalarType,MV,OP>::norm(
      const MV& X, std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > &normvec ) const
  {
    this->normMat(X,normvec);
  }

  template <class ScalarType, class MV, class OP>
  void MatOrthoManager<ScalarType,MV,OP>::normMat(
      const MV& X,
      std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > &normvec,
      Teuchos::RCP<const MV> MX) const
  {
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef Teuchos::ScalarTraits<typename SCT::magnitudeType> MT;
    typedef MultiVecTraits<ScalarType,MV>     MVT;
    typedef OperatorTraits<ScalarType,MV,OP>  OPT;
    
    int nvecs = MVT::GetNumberVecs(X);
    
    // Make sure that normvec has enough entries to hold the norms
    // of all the columns of X.  std::vector<T>::size_type is
    // unsigned, so do the appropriate cast to avoid signed/unsigned
    // comparisons that trigger compiler warnings.
    if (normvec.size() < static_cast<size_t>(nvecs))
      normvec.resize (nvecs);
    
    if (!_hasOp) {
      // X == MX, since the operator M is the identity.
      MX = Teuchos::rcp(&X, false);
      MVT::MvNorm(X, normvec);
    }
    else {
      // The caller didn't give us a previously computed MX, so
      // apply the operator.  We assign to MX only after applying
      // the operator, so that if the application fails, MX won't be
      // modified.
      if(MX == Teuchos::null) {
        Teuchos::RCP<MV> tempVec = MVT::Clone(X,nvecs);
        OPT::Apply(*_Op,X,*tempVec);
        _OpCounter += nvecs;
        MX = tempVec;
      }
      else {
        // The caller gave us a previously computed MX.  Make sure
        // that it has at least as many columns as X.
        const int numColsMX = MVT::GetNumberVecs(*MX);
        TEUCHOS_TEST_FOR_EXCEPTION(numColsMX < nvecs, std::invalid_argument,
                           "MatOrthoManager::norm(X, MX, normvec): "
                           "MX has fewer columns than X: "
                           "MX has " << numColsMX << " columns, "
                           "and X has " << nvecs << " columns.");
      }
      
      std::vector<ScalarType> dotvec(nvecs);
      MVT::MvDot(X,*MX,dotvec);
      for (int i=0; i<nvecs; i++) {
        normvec[i] = MT::squareroot( SCT::magnitude(dotvec[i]) );
      }
    }
  }

  template <class ScalarType, class MV, class OP>
  void MatOrthoManager<ScalarType,MV,OP>::project (
        MV &X,
        Teuchos::Array<Teuchos::RCP<const MV> > Q,
        Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
      ) const
  {
    this->projectMat(X,Q,C);
  }

  template <class ScalarType, class MV, class OP>
  int MatOrthoManager<ScalarType,MV,OP>::normalize (
      MV &X, Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B ) const
  {
    return this->normalizeMat(X,B);
  }

  template <class ScalarType, class MV, class OP>
  int MatOrthoManager<ScalarType,MV,OP>::projectAndNormalize (
        MV &X,
        Teuchos::Array<Teuchos::RCP<const MV> > Q,
        Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C,
        Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B
      ) const
  {
    return this->projectAndNormalizeMat(X,Q,C,B);
  }

  template <class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
  MatOrthoManager<ScalarType,MV,OP>::orthonormError(const MV &X) const
  {
    return this->orthonormErrorMat(X,Teuchos::null);
  }

  template <class ScalarType, class MV, class OP>
  typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
  MatOrthoManager<ScalarType,MV,OP>::orthogError(const MV &X1, const MV &X2) const
  {
    return this->orthogErrorMat(X1,X2);
  }

} // end of Anasazi namespace


#endif

// end of file AnasaziMatOrthoManager.hpp
