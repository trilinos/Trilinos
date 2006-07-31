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
  
  A concrete implementation of this class is necessary. The user can create
  their own implementation if those supplied are not suitable for their needs.
  
  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RefCountPtr.hpp"
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

    \note This can be different than the MvTransMv method from the multivector class. For example,
    if there is a mass matrix \c M, then this might be the \c M inner product (\f$x^HMx\f$)
     */
    virtual void innerProd( const MV &X, const MV &Y, Teuchos::SerialDenseMatrix<int,ScalarType>& Z ) const = 0;


    /*! \brief Provides the norm induced by innerProd().
     */
    virtual void norm( const MV& X, std::vector< typename Teuchos::ScalarTraits<ScalarType>::magnitudeType > * normvec ) const = 0;

    /*! \brief Given a list of (mutually and independently) orthonormal bases \c Q, this method
     * takes a multivector \c X and projects it onto the space orthogonal to the individual \c Q[i], 
     * returning optionally the coefficients of \c X in the individual \c Q[i]. All of this is done with respect
     * to innerProd().
     *  
     @param X [in/out] The multivector to be modified.
       On output, \c X will be orthogonal to \c Q[i] with respect to innerProd().

     @param C [out] The coefficients of \c X in the \c Q[i], with respect to innerProd().

     @param Q [in] A list of multivectors specifying the bases to be orthogonalized against. Each \c Q[i] is assumed to have
     orthonormal columns, and the \c Q[i] are assumed to be mutually orthogonal.

     @return Code specifying failure of the routine, as defined by the implementation.
    */
    virtual void project ( MV &X, 
                           Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
                           Teuchos::Array<Teuchos::RefCountPtr<const MV> > Q) const = 0;

    /*! \brief This method takes a multivector and orthonormalizes the columns, with respect to innerProd().
     *  
     @param X [in/out] The multivector to the modified. 
       On output, \c X will have orthonormal columns with respect to innerProd().

     @param R [out] The coefficients of the original X with respect to the produced basis.

     @return Rank of the basis computed by this method.
    */
    virtual int normalize ( MV &X, Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > R ) const = 0;


    /*! \brief This method takes a multivector and projects it onto the space orthogonal to 
     *  another given multivector.  It also orthonormalizes the 
     *  columns of the resulting multivector. Both of these operations are conducted 
     *  with respect to innerProd().
     *  
     @param X [in/out] The multivector to the modified. 
       On output, \c X will be orthogonal to \c Q and will have orthonormal columns, with respect to innerProd().

     @param C [out] The coefficients of \c X in the \c Q[i], with respect to innerProd().

     @param R [out] The coefficients of the original X with respect to the produced basis.

     @param Q [in] A list of multivectors specifying the bases to be orthogonalized against. Each \c Q[i] is assumed to have
     orthonormal columns, and the \c Q[i] are assumed to be mutually orthogonal.

     @return Rank of the basis computed by this method.
    */
    virtual int projectAndNormalize ( MV &X, 
                                      Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > C, 
                                      Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > R, 
                                      Teuchos::Array<Teuchos::RefCountPtr<const MV> > Q ) const = 0;

    //@}

    //! @name Error methods
    //@{ 

    /*! \brief This method computes the error in orthonormality of a multivector.
     */
    virtual typename Teuchos::ScalarTraits< ScalarType >::magnitudeType orthonormError(const MV &X) const = 0;

    /*! \brief This method computes the error in orthogonality of two multivectors.
     */
    virtual typename Teuchos::ScalarTraits<ScalarType>::magnitudeType orthogError(const MV &X1, const MV &X2) const = 0;

    //@}

  };

} // end of Anasazi namespace


#endif

// end of file AnasaziOrthoManager.hpp
