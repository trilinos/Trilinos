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

/*! \file BelosMatOrthoManager.hpp
  \brief  Templated virtual class for providing orthogonalization/orthonormalization methods with matrix-based 
          inner products.
*/

#ifndef BELOS_MATORTHOMANAGER_HPP
#define BELOS_MATORTHOMANAGER_HPP

/*! \class Belos::MatOrthoManager
  
  \brief Belos's templated virtual class for providing routines for orthogonalization and 
  orthonormzalition of multivectors using matrix-based inner products.

  This class extends Belos::OrthoManager by providing extra calling arguments to orthogonalization
  routines, to reduce the cost of applying the inner product in cases where the user already
  has the image of the source multivector under the inner product matrix.

  A concrete implementation of this class is necessary. The user can create
  their own implementation if those supplied are not suitable for their needs.
  
  \author Chris Baker, Teri Barth, and Heidi Thornquist
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosOrthoManager.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosOperatorTraits.hpp"

namespace Belos {

  template <class ScalarType, class MV, class OP>
  class MatOrthoManager : public OrthoManager<ScalarType,MV> {
  protected:
    Teuchos::RCP<const OP> _Op;
    bool _hasOp;

  public:
    //! @name Constructor/Destructor
    //@{ 
    //! Default constructor.
    MatOrthoManager(Teuchos::RCP<const OP> Op = Teuchos::null) : _Op(Op), _hasOp(Op!=Teuchos::null) {};

    //! Destructor.
    virtual ~MatOrthoManager() {};
    //@}

    //! @name Accessor routines
    //@{ 

    //! Set operator.
    void setOp( Teuchos::RCP<const OP> Op ) { 
      _Op = Op; 
      _hasOp = (_Op != Teuchos::null);
    };

    //! Get operator.
    Teuchos::RCP<const OP> getOp() const { return _Op; } 

    //@}


    //! @name Orthogonalization methods
    //@{ 

    /*! \brief Provides the inner product defining the orthogonality concepts, using the provided operator.

    All concepts of orthogonality discussed in this class are with respect to this inner product.
     */
    void innerProd( const MV& X, const MV& Y, 
                                  Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const {
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
          P = R;
          Q = Teuchos::rcp( &Y, false );
        }
        else {
          P = Teuchos::rcp( &X, false );
          R = MVT::Clone(Y,MVT::GetNumberVecs(Y));
          OPT::Apply(*_Op,Y,*R);
          Q = R;
        }
      }
      else {
        P = Teuchos::rcp( &X, false );
        Q = Teuchos::rcp( &Y, false );
      }
      
      MVT::MvTransMv(SCT::one(),*P,*Q,Z);
    }

    /*! \brief Provides the inner product defining the orthogonality concepts, using the provided operator.
     *  The method has the options of exploiting a caller-provided \c MX.
     *
     *  If pointer \c MY is null, then this routine calls innerProd(X,Y,Z). Otherwise, it forgoes the 
     *  operator application and uses \c *MY in the computation of the inner product.
     */
    void innerProd( const MV& X, const MV& Y, Teuchos::RCP<const MV> MY, 
                            Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const {
      typedef Teuchos::ScalarTraits<ScalarType> SCT;
      typedef MultiVecTraits<ScalarType,MV>     MVT;

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
    }

    /*! \brief Provides the norm induced by innerProd().
     */
    void norm( const MV& X, std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType >& normvec ) const {
      norm(X,Teuchos::null,normvec);
    }

    /// \brief Compute norm of each column of X
    ///
    /// Compute the norm of each column of X; where the norm is that
    /// induced by innerProd().  If you already have MX = M*X (where M
    /// is the operator defining the inner product), you may pass it in
    /// to avoid applying the operator again.
    ///
    /// \param X [in] Multivector for which to compute the norm of
    ///   each column
    /// \param MX [in] If not null and the inner product operator M is
    ///   not the identity, MX is assumed to be M*X (the result of
    ///   applying the operator M to X).  MX may be null, in which
    ///   case if M is not the identity, it is applied to X.
    /// \param normvec [out] On output, normvec[j] is the norm of
    ///   column j of X.  normvec is resized if it has fewer entries
    ///   than the number of columns in X.
    ///
    void 
    norm (const MV& X, 
	  Teuchos::RCP<const MV> MX, 
	  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& normvec) const 
    {
      typedef Teuchos::ScalarTraits<ScalarType> SCT;
      typedef MultiVecTraits<ScalarType,MV>     MVT;
      typedef OperatorTraits<ScalarType,MV,OP>  OPT;

      const int numColsX = MVT::GetNumberVecs(X);
      if (!_hasOp) {
	// X == MX, since the operator M is the identity.
        MX = Teuchos::rcp(&X, false);
      }
      else if (MX.is_null()) {
	// The caller didn't give us a previously computed MX, so
	// apply the operator.  We assign to MX only after applying
	// the operator, so that if the application fails, MX won't be
	// modified.
        Teuchos::RCP<MV> R = MVT::Clone(X, numColsX);
        OPT::Apply(*_Op,X,*R);
        MX = R;
      } 
      else {
	// The caller gave us a previously computed MX.  Make sure
	// that it has at least as many columns as X.
	const int numColsMX = MVT::GetNumberVecs(*MX);
	TEUCHOS_TEST_FOR_EXCEPTION(numColsMX < numColsX, std::invalid_argument,
			   "MatOrthoManager::norm(X, MX, normvec): "
			   "MX has fewer columns than X: "
			   "MX has " << numColsMX << " columns, "
			   "and X has " << numColsX << " columns.");
      }

      // Make sure that normvec has enough entries to hold the norms
      // of all the columns of X.  std::vector<T>::size_type is
      // unsigned, so do the appropriate cast to avoid signed/unsigned
      // comparisons that trigger compiler warnings.
      if (normvec.size() < static_cast<size_t>(numColsX))
	normvec.resize (numColsX);

      Teuchos::SerialDenseMatrix<int,ScalarType> z(1,1);
      Teuchos::RCP<const MV> Xi, MXi;
      std::vector<int> ind(1);
      for (int i = 0; i < numColsX; ++i) {
        ind[0] = i;
        Xi = MVT::CloneView(X,ind);
        MXi = MVT::CloneView(*MX,ind);
        MVT::MvTransMv(SCT::one(),*Xi,*MXi,z);
        normvec[i] = SCT::magnitude( SCT::squareroot( z(0,0) ) );
      }
    }


    /*! \brief Given a list of (mutually and internally) orthonormal bases \c Q, this method
     * takes a multivector \c X and projects it onto the space orthogonal to the individual <tt>Q[i]</tt>, 
     * optionally returning the coefficients of \c X for the individual <tt>Q[i]</tt>. All of this is done with respect
     * to the inner product innerProd().
     *  
     * After calling this routine, \c X will be orthogonal to each of the <tt>Q[i]</tt>.
     *
     @param X [in/out] The multivector to be modified.
       On output, \c X will be orthogonal to <tt>Q[i]</tt> with respect to innerProd().

     @param MX [in] The image of the multivector under the specified operator. If \c MX is null, it is not used.

     @param C [out] The coefficients of \c X in the \c *Q[i], with respect to innerProd(). If <tt>C[i]</tt> is a non-null pointer 
       and \c *C[i] matches the dimensions of \c X and \c *Q[i], then the coefficients computed during the orthogonalization
       routine will be stored in the matrix \c *C[i]. If <tt>C[i]</tt> is a non-null pointer whose size does not match the dimensions of 
       \c X and \c *Q[i], then a std::invalid_argument std::exception will be thrown. Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null
       pointer, then the orthogonalization manager will declare storage for the coefficients and the user will not have access to them.

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against. Each <tt>Q[i]</tt> is assumed to have
     orthonormal columns, and the <tt>Q[i]</tt> are assumed to be mutually orthogonal.
    */
    virtual void project ( MV &X, Teuchos::RCP<MV> MX, 
                           Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
                           Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const = 0;


    
    /*! \brief This method calls project(X,Teuchos::null,C,Q); see documentation for that function.
    */
    virtual void project ( MV &X, 
                           Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
                           Teuchos::ArrayView<Teuchos::RCP<const MV> > Q) const {
      project(X,Teuchos::null,C,Q);
    }

    /*! \brief This method takes a multivector \c X and attempts to compute an orthonormal basis for \f$colspan(X)\f$, with respect to innerProd().
     *
     * This routine returns an integer \c rank stating the rank of the computed basis. If \c X does not have full rank and the normalize() routine does 
     * not attempt to augment the subspace, then \c rank may be smaller than the number of columns in \c X. In this case, only the first \c rank columns of 
     * output \c X and first \c rank rows of \c B will be valid.
     *  
     @param X [in/out] The multivector to the modified. 
       On output, \c X will have some number of orthonormal columns (with respect to innerProd()).

     @param MX [in/out] The image of the multivector under the specified operator. If \c MX is null, it is not used.
                        On output, it returns the image of the valid basis vectors under the specified operator.

     @param B [out] The coefficients of \c X in the computed basis. If \c B is a non-null pointer 
       and \c *B has appropriate dimensions, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix \c *B. If \c B is a non-null pointer whose size does not match the dimensions of 
       \c X, then a std::invalid_argument std::exception will be thrown. Otherwise, 
       the orthogonalization manager will declare storage for the coefficients and the user will not have access to them. <b>This matrix may or may not be triangular; see 
       documentation for individual orthogonalization managers.</b>

     @return Rank of the basis computed by this method.
    */
    virtual int normalize ( MV &X, Teuchos::RCP<MV> MX, 
                            Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B ) const = 0;


    /*! \brief This method calls normalize(X,Teuchos::null,B); see documentation for that function.
    */
    virtual int normalize ( MV &X, Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B ) const {
      return normalize(X,Teuchos::null,B);
    }


    /*! \brief Given a set of bases <tt>Q[i]</tt> and a multivector \c X, this method computes an orthonormal basis for \f$colspan(X) - \sum_i colspan(Q[i])\f$.
     *
     *  This routine returns an integer \c rank stating the rank of the computed basis. If the subspace \f$colspan(X) - \sum_i colspan(Q[i])\f$ does not 
     *  have dimension as large as the number of columns of \c X and the orthogonalization manager doe not attempt to augment the subspace, then \c rank 
     *  may be smaller than the number of columns of \c X. In this case, only the first \c rank columns of output \c X and first \c rank rows of \c B will 
     *  be valid.
     *
     * \note This routine guarantees both the orthgonality constraints against the <tt>Q[i]</tt> as well as the orthonormality constraints. Therefore, this method 
     * is not necessarily equivalent to calling project() followed by a call to normalize(); see the documentation for specific orthogonalization managers.
     *
     @param X [in/out] The multivector to the modified. 
       On output, the relevant rows of \c X will be orthogonal to the <tt>Q[i]</tt> and will have orthonormal columns (with respect to innerProd()).

     @param MX [in/out] The image of the multivector under the specified operator. If \c MX is null, it is not used.
                        On output, it returns the image of the valid basis vectors under the specified operator.

     @param C [out] The coefficients of the original \c X in the \c *Q[i], with respect to innerProd(). If <tt>C[i]</tt> is a non-null pointer 
       and \c *C[i] matches the dimensions of \c X and \c *Q[i], then the coefficients computed during the orthogonalization
       routine will be stored in the matrix \c *C[i]. If <tt>C[i]</tt> is a non-null pointer whose size does not match the dimensions of 
       \c X and \c *Q[i], then a std::invalid_argument std::exception will be thrown. Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null
       pointer, then the orthogonalization manager will declare storage for the coefficients and the user will not have access to them.

     @param B [out] The coefficients of \c X in the computed basis. If \c B is a non-null pointer 
       and \c *B has appropriate dimensions, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix \c *B. If \c B is a non-null pointer whose size does not match the dimensions of 
       \c X, then a std::invalid_argument std::exception will be thrown. Otherwise, 
       the orthogonalization manager will declare storage for the coefficients and the user will not have access to them. <b>This matrix may or may not be triangular; see 
       documentation for individual orthogonalization managers.</b>

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against. Each <tt>Q[i]</tt> is assumed to have
     orthonormal columns, and the <tt>Q[i]</tt> are assumed to be mutually orthogonal.

     @return Rank of the basis computed by this method.
    */
    virtual int projectAndNormalize ( MV &X, Teuchos::RCP<MV> MX, 
                                      Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
                                      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B, 
                                      Teuchos::ArrayView<Teuchos::RCP<const MV> > Q ) const = 0;

    /*! \brief This method calls projectAndNormalize(X,Teuchos::null,C,B,Q); see documentation for that function.
    */
    virtual int projectAndNormalize ( MV &X, 
                                      Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
                                      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B, 
                                      Teuchos::ArrayView<Teuchos::RCP<const MV> > Q ) const {
      return projectAndNormalize(X,Teuchos::null,C,B,Q);
    }

    //@}

    //! @name Error methods
    //@{ 

    /*! \brief This method computes the error in orthonormality of a multivector.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthonormError(const MV &X) const {
      return orthonormError(X,Teuchos::null);
    }

    /*! \brief This method computes the error in orthonormality of a multivector.
     *  The method has the option of exploiting a caller-provided \c MX.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthonormError(const MV &X, Teuchos::RCP<const MV> MX) const = 0;

    /*! \brief This method computes the error in orthogonality of two multivectors. This method 
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthogError(const MV &X1, const MV &X2) const {
      return orthogError(X1,Teuchos::null,X2);
    }

    /*! \brief This method computes the error in orthogonality of two multivectors.
     *  The method has the option of
     *  exploiting a caller-provided \c MX.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType 
    orthogError(const MV &X1, Teuchos::RCP<const MV> MX1, const MV &X2) const = 0;

    //@}

  };

} // end of Belos namespace


#endif

// end of file BelosMatOrthoManager.hpp
