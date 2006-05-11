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
  
  A concrete implementation of this class is necessary.  The user can create
  their own implementation if those supplied are not suitable for their needs.
  
  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziTypes.hpp"

namespace Anasazi {

  template <class ScalarType, class MV, class OP>
  class OrthoManager {
  public:
    //@{ \name Constructor/Destructor.
    //! Default constructor.
    OrthoManager() {};

    //! Destructor.
    virtual ~OrthoManager() {};
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
    virtual ReturnType project ( MV &X, MV &MX, const OP *M, const MV &Q ) const = 0;



    /*! \brief This method takes a multivector and orthonormalizes the columns using a specified inner product.
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
    virtual ReturnType normalize ( MV &X, MV &MX, const OP *M, int &rank ) const = 0;



    /*! \brief This method takes a multivector and projects it onto the space orthogonal to 
     *  another given multivector, in a specified inner product. It also orthonormalizes the 
     *  columns of the resulting multivector.
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
    virtual ReturnType projectAndNormalize ( MV &X, MV &MX, const OP *M, const MV &Q, int &rank ) const = 0;

    //@}
  };

} // end of Anasazi namespace


#endif

// end of file AnasaziOrthoManager.hpp
