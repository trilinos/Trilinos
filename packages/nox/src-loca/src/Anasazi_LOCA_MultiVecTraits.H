// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef ANASAZI_LOCA_MULTIVECTRAITS_H
#define ANASAZI_LOCA_MULTIVECTRAITS_H

// Anasazi includes
#include "AnasaziMultiVecTraits.hpp"

// NOX/LOCA includes
#include "NOX_Abstract_MultiVector.H"

namespace Anasazi {
  
  /*!
   * \brief Implementation of the Anasazi::MultiVecTraits for 
   * NOX::Abstract::MultiVectors.
   */
  /*!
   * This class provides an implementation of the Anasazi::MultiVecTraits class
   * for NOX::Abstract::MultiVector's, to be used in conjunction with the
   * Anasazi::LOCA::Matrix class.
   */
  template <>
  class MultiVecTraits<double,NOX::Abstract::MultiVector> {
  public:

    //! Synony for double
    typedef double ScalarType;

    //! Synonym for NOX::Abstract::MultiVector
    typedef NOX::Abstract::MultiVector MV;

    //@{ \name Creation methods

    /*! 
     * \brief Creates a new empty \c MV containing \c numvecs columns.
     *
     * \return Reference-counted pointer to the new multivector of type \c MV.
     */
    static Teuchos::RefCountPtr<MV> Clone( const MV& mv, const int numvecs ) 
    {
      return mv.clone(numvecs);
    }
      
    /*! 
     * \brief Creates a new \c MV and copies contents of \c mv into the new 
     * vector (deep copy).
     *
     * \return Reference-counted pointer to the new multivector of type \c MV.
     */
    static Teuchos::RefCountPtr<MV> CloneCopy( const MV& mv ) 
    {
      return mv.clone(NOX::DeepCopy);
    }

    /*! 
     * \brief Creates a new \c MV and copies the selected contents of \c mv 
     * into the new vector (deep copy).  
     *
     * The copied vectors from \c mv are indicated by the \c index.size() 
     * indices in \c index.      
     * \return Reference-counted pointer to the new multivector of type \c MV.
     */
    static Teuchos::RefCountPtr<MV> CloneCopy( const MV& mv, 
					       const std::vector<int>& index )
    {
      return mv.subCopy(index);
    }

    /*! 
     * \brief Creates a new \c MV that shares the selected contents of 
     * \c mv (shallow copy).
     * 
     * The index of the \c numvecs vectors shallow copied from \c mv are 
     * indicated by the indices given in \c index.
     * \return Reference-counted pointer to the new multivector of type \c MV.
     */      
    static Teuchos::RefCountPtr<MV> CloneView( MV& mv, 
					       const std::vector<int>& index )
    {
      return mv.subView(index);
    }

    /*! 
     * \brief Creates a new const \c MV that shares the selected contents of 
     * \c mv (shallow copy).
     *
     * The index of the \c numvecs vectors shallow copied from \c mv are 
     * indicated by the indices given in \c index.
     * \return Reference-counted pointer to the new const multivector of 
     * type \c MV.
     */      
    static Teuchos::RefCountPtr<const MV> CloneView( const MV& mv, 
						     const std::vector<int>& index ) 
    {
      return mv.subView(index);
    }

    //@}

    //@{ \name Attribute methods

    //! Obtain the vector length of \c mv.
    static int GetVecLength( const MV& mv ) 
    {
      return mv.length();
    }

    //! Obtain the number of vectors in \c mv
    static int GetNumberVecs( const MV& mv ) 
    {
      return mv.numVectors();
    }

    //@}

    //@{ \name Update methods

    /*! 
     * \brief Update \c mv with \f$ \alpha AB + \beta mv \f$.
     */
    static void MvTimesMatAddMv( const ScalarType alpha, const MV& A, 
				 const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
				 const ScalarType beta, MV& mv ) 
    {
      mv.update(Teuchos::NO_TRANS, alpha, A, B, beta);
    }

    /*! 
     * \brief Replace \c mv with \f$\alpha A + \beta B\f$.
     */
    static void MvAddMv( const ScalarType alpha, const MV& A, 
			 const ScalarType beta, const MV& B, MV& mv ) 
    {
      mv.update(alpha, A, beta, B);
    }

    /*! 
     * \brief Compute a dense matrix \c B through the matrix-matrix multiply 
     * \f$ \alpha A^Hmv \f$.
     */
    static void MvTransMv( const ScalarType alpha, const MV& A, const MV& mv, 
			   Teuchos::SerialDenseMatrix<int,ScalarType>& B) 
    {
      mv.multiply(alpha, A, B);
    }
    
    /*! 
     * \brief Compute a vector \c b where the components are the individual 
     * dot-products of the \c i-th columns of \c A and \c mv, i.e.
     * \f$b[i] = A[i]^Hmv[i]\f$.
     */
    static void MvDot ( const MV& mv, const MV& A, 
			std::vector<ScalarType>* b) 
    {
      for (unsigned int i=0; i<b->size(); i++)
	(*b)[i] = A[i].innerProduct(mv[i]);
    }

    //@}
    //@{ \name Norm method

    /*! 
     * \brief Compute the 2-norm of each individual vector of \c mv.  
     * Upon return, \c normvec[i] holds the value of \f$||mv_i||_2\f$, 
     * the \c i-th column of \c mv.
     */
    static void MvNorm( const MV& mv, 
			std::vector<double>* normvec ) 
    {
      mv.norm(*normvec);
    }

    //@}

    //@{ \name Initialization methods

    /*! 
     * \brief Copy the vectors in \c A to a set of vectors in \c mv indicated 
     * by the indices given in \c index.
     *
     * The \c numvecs vectors in \c A are copied to a subset of vectors in 
     * \c mv indicated by the indices given in \c index,
     * i.e.<tt> mv[index[i]] = A[i]</tt>.
     */
    static void SetBlock( const MV& A, const std::vector<int>& index, MV& mv )
    {
      mv.setBlock(A, index);
    }

    /*! 
     * \brief Replace the vectors in \c mv with random vectors.
     */
    static void MvRandom( MV& mv ) 
    {
      mv.random();
    }

    /*! 
     * \brief Replace each element of the vectors in \c mv with \c alpha.
     */
    static void MvInit( MV& mv, 
			const ScalarType alpha = 0.0 ) 
    {
      mv.init(alpha);
    }

    //@}

    //@{ \name Print method

    /*! 
     * \brief Print the \c mv multi-vector to the \c os output stream.
     */
    static void MvPrint( const MV& mv, ostream& os ) 
    {
      mv.print(os);
    }

    //@}
  }; // class MultiVecTraits
} // namespace Anasazi

#endif 
