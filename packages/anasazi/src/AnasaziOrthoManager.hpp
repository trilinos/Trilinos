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

/*! \file AnasaziOrthoManager.hpp
  \brief  Templated virtual class for providing orthogonalization/orthonormalization methods.
*/

#ifndef ANASAZI_ORTHOMANAGER_HPP
#define ANASAZI_ORTHOMANAGER_HPP

/*! \class Anasazi::OrthoManager
  
  \brief Anasazi's templated virtual class for providing routines for orthogonalization and 
  orthonormzalition of multivectors. 

  This class defines concepts of orthogonality through the definition of an
  inner product. It also provides computational routines for orthogonalization.

  A concrete implementation of this class is necessary. The user can create
  their own implementation if those supplied are not suitable for their needs.
  
  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_Array.hpp"




namespace Anasazi {


  //! @name OrthoManager Exceptions
  //@{ 

  /** \brief Exception thrown to signal error in an orthogonalization manager method.
   */
  class OrthoError : public AnasaziError
  {public: OrthoError(const std::string& what_arg) : AnasaziError(what_arg) {}};

  //@}

  template <class ScalarType, class MV>
  class OrthoManager {
  public:
    //! @name Constructor/Destructor
    //@{ 
    //! Default constructor.
    OrthoManager() {};

    //! Destructor.
    virtual ~OrthoManager() {};
    //@}

    //! @name Orthogonalization methods
    //@{ 

    /*! \brief Provides the inner product defining the orthogonality concepts.

    All concepts of orthogonality discussed in this class are defined with respect to this inner product.

    \note This can be different than the \c MvTransMv method from the multivector class. For example,
    if there is a mass matrix \c M, then this might be the \c M inner product: \f$x^HMx\f$. 

    @param Z [out] <tt>Z(i,j)</tt> contains the inner product of <tt>X[i]</tt> and <tt>Y[i]</tt>:
    \f[
      Z(i,j) = \langle X[i], Y[i] \rangle
    \f]
    
     */
    virtual void innerProd( const MV &X, const MV &Y, Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const = 0;


    /*! \brief Provides the norm induced by innerProd().
     *
     * This computes the norm for each column of a multivector. The norm is that induced by innerProd(), where
     * \f[ \|x\| = \sqrt{\langle x, x \rangle} \f]
     *
     * @param normvec [out] Pointer to a vector of norms, whose \c i-th entry corresponds to the \c i-th column of \c X
     */
    virtual void norm( const MV& X, std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > * normvec ) const = 0;

    /*! \brief Given a list of (mutually and internally) orthonormal bases \c Q, this method
     * takes a multivector \c X and projects it onto the space orthogonal to the individual <tt>Q[i]</tt>, 
     * optionally returning the coefficients of \c X for the individual <tt>Q[i]</tt>. All of this is done with respect
     * to the inner product innerProd().
     *  
     * After calling this routine, \c X will be orthogonal to each of the <tt>Q[i]</tt>.
     *
     @param X [in/out] The multivector to be modified.
       On output, the columns of \c X will be orthogonal to each <tt>Q[i]</tt> 

     @param C [out] The coefficients of \c X in the <tt>Q[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c X and <tt>Q[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>, similar to calling
       \code
          innerProd( Q[i], X, C[i] );
       \endcode
       If <tt>C[i]</tt> is a non-null pointer whose size does not match the dimensions of 
       \c X and \c <tt>Q[i]</tt>, then a std::invalid_argument exception will be thrown. Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null
       pointer, the caller will not have access to the computed coefficients.

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against. Each <tt>Q[i]</tt> is assumed to have
     orthonormal columns, and the <tt>Q[i]</tt> are assumed to be mutually orthogonal.
    */
    virtual void project ( MV &X, 
                           Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
                           Teuchos::Array<Teuchos::RCP<const MV> > Q) const = 0;

    /*! \brief This method takes a multivector \c X and attempts to compute an orthonormal basis for \f$colspan(X)\f$, with respect to innerProd().
     *
     * This routine returns an integer \c rank stating the rank of the computed basis. If \c X does not have full rank and the normalize() routine does 
     * not attempt to augment the subspace, then \c rank may be smaller than the number of columns in \c X. In this case, only the first \c rank columns of 
     * output \c X and first \c rank rows of \c B will be valid.
     *  
     @param X [in/out] 
       On output, \c X will have some number of orthonormal columns equal to the return value of normalize():
       \f[
          \langle X[i], X[j] \rangle = \delta_{ij}
       \f]

     @param B [out] The coefficients of the original \c X with respect to the computed basis. This matrix is not necessarily triangular; see the documentation
       for specific orthogonalization managers, similar to calling
       \code
          innerProd( Xout, Xin, B );
       \endcode

     @return Rank of the basis computed by this method, less than or equal to the number of columns in \c X.
    */
    virtual int normalize ( MV &X, Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B ) const = 0;


    /*! \brief Given a set of bases <tt>Q[i]</tt> and a multivector \c X, this method computes an orthonormal basis for \f$colspan(X) - \sum_i colspan(Q[i])\f$.
     *
     *  This routine returns an integer \c rank stating the rank of the computed basis. If the subspace \f$colspan(X) - \sum_i colspan(Q[i])\f$ does not 
     *  have dimension as large as the number of columns of \c X and the orthogonalization manager does not attempt to augment the subspace, then \c rank 
     *  may be smaller than the number of columns of \c X. In this case, only the first \c rank columns of output \c X and first \c rank rows of \c B will 
     *  be valid.
     *
     * \note This routine guarantees both the orthgonality of the returned basis against the <tt>Q[i]</tt> as well as the orthonormality of the returned basis. Therefore, 
     * this method is not necessarily equivalent to calling project() followed by a call to normalize(); see the documentation for specific orthogonalization managers.
     *
     @param X [in/out] 
       On output, the first \c rank columns of \c X satisfy
       \f[
            \langle X[i], X[j] \rangle = \delta_{ij} \quad \textrm{and} \quad \langle X, Q[i] \rangle = 0\ .
       \f]
       and
       \f[
          X_{in} = X_{out} B + \sum_i Q[i] C[i] 
       \f]

     @param C [out] The coefficients of \c X in the <tt>Q[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c X and <tt>Q[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>, similar to calling
       \code
          innerProd( Q[i], X, C[i] );
       \endcode
       If <tt>C[i]</tt> is a non-null pointer whose size does not match the dimensions of 
       \c X and \c <tt>Q[i]</tt>, then a std::invalid_argument exception will be thrown. Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null
       pointer, the caller will not have access to the computed coefficients.

     @param B [out] The coefficients of the original \c X with respect to the computed basis. This matrix is not necessarily triangular; see the documentation
       for specific orthogonalization managers, similar to calling
       \code
          innerProd( Xout, Xin, B );
       \endcode

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against. Each <tt>Q[i]</tt> is assumed to have
     orthonormal columns, and the <tt>Q[i]</tt> are assumed to be mutually orthogonal.

     @return Rank of the basis computed by this method.
    */
    virtual int projectAndNormalize ( MV &X, 
                                      Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
                                      Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B, 
                                      Teuchos::Array<Teuchos::RCP<const MV> > Q ) const = 0;

    //@}

    //! @name Error methods
    //@{ 

    /*! \brief This method computes the error in orthonormality of a multivector.
     *
     * This method return some measure \f$\| \langle X, X \rangle - I \| \f$. See the documentation for the specific orthogonalization manager.
     */
    virtual typename Teuchos::ScalarTraits< ScalarType >::magnitudeType orthonormError(const MV &X) const = 0;

    /*! \brief This method computes the error in orthogonality of two multivectors.
     *
     * This method return some measure \f$\| \langle X1, X2 \rangle - 0 \| \f$. See the documentation for the specific orthogonalization manager.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType orthogError(const MV &X1, const MV &X2) const = 0;

    //@}

  };

} // end of Anasazi namespace


#endif

// end of file AnasaziOrthoManager.hpp
