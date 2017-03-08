// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ***********************************************************************
// @HEADER

/*! \file AnasaziOrthoManager.hpp
  \brief  Templated virtual class for providing orthogonalization/orthonormalization methods.
*/

#ifndef ANASAZI_ORTHOMANAGER_HPP
#define ANASAZI_ORTHOMANAGER_HPP

/*! \class Anasazi::OrthoManager
  
  \brief Anasazi's templated virtual class for providing routines for orthogonalization and 
  orthonormalization of multivectors. 

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

    \note This is potentially different from MultiVecTraits::MvTransMv(). For example, it is customary in many
    eigensolvers to exploit a mass matrix \c M for the inner product: \f$x^HMx\f$. 

    @param Z [out] <tt>Z(i,j)</tt> contains the inner product of <tt>X[i]</tt> and <tt>Y[i]</tt>:
    \f[
      Z(i,j) = \langle X[i], Y[i] \rangle
    \f]
    
     */
    virtual void innerProd( const MV &X, const MV &Y, Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const = 0;


    /*! \brief Provides the norm induced by innerProd().
     *
     * This computes the norm for each column of a multivector. This is the norm induced by innerProd():
     * \f[ \|x\| = \sqrt{\langle x, x \rangle} \f]
     *
     @param normvec [out] Vector of norms, whose \c i-th entry corresponds to the \c i-th column of \c X
      
     \pre 
     <ul>
     <li><tt>normvec.size() == GetNumberVecs(X)</tt>
     </ul>
     */
    virtual void norm( const MV& X, std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > &normvec ) const = 0;



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
          \langle Q[i], X_{out} \rangle = 0
       \f]
       Also, 
       \f[
          X_{out} = X_{in} - \sum_i Q[i] \langle Q[i], X_{in} \rangle
       \f]

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against, satisfying 
     \f[
        \langle Q[i], Q[j] \rangle = I \quad\textrm{if}\quad i=j
     \f]
     and
     \f[
        \langle Q[i], Q[j] \rangle = 0 \quad\textrm{if}\quad i \neq j\ .
     \f]

     @param C [out] The coefficients of \c X in the bases <tt>Q[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c X and <tt>Q[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>, similar to calling
       \code
          innerProd( Q[i], X, C[i] );
       \endcode
       If <tt>C[i]</tt> points to a Teuchos::SerialDenseMatrix with size
       inconsistent with \c X and \c <tt>Q[i]</tt>, then a std::invalid_argument
       exception will be thrown.<br>
       Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null pointer,
       the caller will not have access to the computed coefficients.
    */
    virtual void project ( 
          MV &X, 
          Teuchos::Array<Teuchos::RCP<const MV> > Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C 
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null))
        ) const = 0;



    /*! \brief This method takes a multivector \c X and attempts to compute a basis 
     * for \f$colspan(X)\f$. This basis is orthonormal with respect to innerProd().
     *
     * This routine returns an integer \c rank stating the rank of the computed
     * basis. If \c X does not have full rank and the normalize() routine does
     * not attempt to augment the subspace, then \c rank may be smaller than
     * the number of columns in \c X. In this case, only the first \c rank
     * columns of output \c X and first \c rank rows of \c B will be valid.
     *
     @param X [in/out] The multivector to be modified.<br>
       On output, the first \c rank columns of \c X satisfy
       \f[
          \langle X[i], X[j] \rangle = \delta_{ij}\ .
       \f]
       Also,
       \f[
          X_{in}(1:m,1:n) = X_{out}(1:m,1:rank) B(1:rank,1:n)\ ,
       \f]
       where \c m is the number of rows in \c X and \c n is the number of
       columns in \c X.

     @param B [out] The coefficients of the original \c X with respect to the
       computed basis. If \c B is a non-null pointer and \c B matches the
       dimensions of \c B, then the coefficients computed during the
       orthogonalization routine will be stored in \c B, similar to calling 
       \code
          innerProd( X_{out}, X_{in}, B );
       \endcode
       If \c B points to a Teuchos::SerialDenseMatrix with
       size inconsistent with \c X, then a std::invalid_argument exception will
       be thrown.<br>
       Otherwise, if \c B is null, the caller will not have access
       to the computed coefficients.<br>
       \note This matrix is not necessarily triangular (as in a QR factorization);
       see the documentation of specific orthogonalization managers.

     @return Rank of the basis computed by this method, less than or equal to
       the number of columns in \c X. This specifies how many columns in the
       returned \c X and rows in the returned \c B are valid.
    */
    virtual int normalize ( 
        MV &X, 
        Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null) const = 0;


    /*! \brief Given a set of bases <tt>Q[i]</tt> and a multivector \c X, this method computes an orthonormal basis for \f$colspan(X) - \sum_i colspan(Q[i])\f$.
     *    
     *  This routine returns an integer \c rank stating the rank of the
     *  computed basis. If the subspace \f$colspan(X) - \sum_i colspan(Q[i])\f$
     *  does not have dimension as large as the number of columns of \c X and
     *  the orthogonalization manager does not attempt to augment the subspace,
     *  then \c rank may be smaller than the number of columns of \c X. In this
     *  case, only the first \c rank columns of output \c X and first \c rank
     *  rows of \c B will be valid.
     *
     * \note This routine guarantees both the orthogonality of the returned
     * basis against the <tt>Q[i]</tt> as well as the orthonormality of the
     * returned basis. Therefore, this method is not necessarily equivalent to
     * calling project() followed by a call to normalize(); see the
     * documentation for specific orthogonalization managers.
     * 
     @param X [in/out] 
       On output, the first \c rank columns of \c X satisfy
       \f[
            \langle X[i], X[j] \rangle = \delta_{ij} \quad \textrm{and} \quad \langle X, Q[i] \rangle = 0\ .
       \f]
       Also, 
       \f[
          X_{in}(1:m,1:n) = X_{out}(1:m,1:rank) B(1:rank,1:n) + \sum_i Q[i] C[i]
       \f]
       where \c m is the number of rows in \c X and \c n is the number of columns in \c X.

     @param Q [in] A list of multivector bases specifying the subspaces to be orthogonalized against, satisfying 
     \f[
        \langle Q[i], Q[j] \rangle = I \quad\textrm{if}\quad i=j
     \f]
     and
     \f[
        \langle Q[i], Q[j] \rangle = 0 \quad\textrm{if}\quad i \neq j\ .
     \f]

     @param C [out] The coefficients of \c X in the <tt>Q[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c X and <tt>Q[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>, similar to calling
       \code
          innerProd( Q[i], X, C[i] );
       \endcode
       If <tt>C[i]</tt> points to a Teuchos::SerialDenseMatrix with size
       inconsistent with \c X and \c <tt>Q[i]</tt>, then a std::invalid_argument
       exception will be thrown.<br>
       Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null pointer,
       the caller will not have access to the computed coefficients.

     @param B [out] The coefficients of the original \c X with respect to the computed basis. If \c B is a non-null pointer and \c B matches the dimensions of \c B, then the
     coefficients computed during the orthogonalization routine will be stored in \c B, similar to calling 
       \code
          innerProd( Sout, Sin, B );
       \endcode
     If \c B points to a Teuchos::SerialDenseMatrix with size inconsistent with
     \c X, then a std::invalid_argument exception will be thrown.<br>
     Otherwise, if \c B is null, the caller will not have access to the computed
     coefficients.<br>
     \note This matrix is not necessarily triangular (as in a QR
     factorization); see the documentation of specific orthogonalization
     managers.

     @return Rank of the basis computed by this method, less than or equal to
       the number of columns in \c X. This specifies how many columns in the
       returned \c X and rows in the returned \c B are valid.
    */
    virtual int projectAndNormalize ( 
          MV &X, 
          Teuchos::Array<Teuchos::RCP<const MV> > Q,
          Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C 
              = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
          Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null
        ) const = 0;

    //@}

    //! @name Error methods
    //@{ 

    /*! \brief This method computes the error in orthonormality of a multivector.
     *
     * This method return some measure of \f$\| \langle X, X \rangle - I \| \f$. <br>
     * See the documentation of specific orthogonalization managers.
     */
    virtual typename Teuchos::ScalarTraits< ScalarType >::magnitudeType orthonormError(const MV &X) const = 0;

    /*! \brief This method computes the error in orthogonality of two multivectors.
     *
     * This method return some measure of \f$\| \langle X1, X2 \rangle \| \f$. <br>
     * See the documentation of specific orthogonalization managers.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType orthogError(const MV &X1, const MV &X2) const = 0;

    //@}

  };

} // end of Anasazi namespace


#endif

// end of file AnasaziOrthoManager.hpp
