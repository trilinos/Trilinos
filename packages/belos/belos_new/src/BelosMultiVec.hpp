// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef BELOS_MULTI_VEC_HPP
#define BELOS_MULTI_VEC_HPP

#include "BelosMultiVecTraits.hpp"
#include "BelosTypes.hpp"
#include "BelosConfigDefs.hpp"

/*! 	\class Belos::MultiVec

	\brief Belos's templated pure virtual class for constructing multivectors that 
	are used by the eigensolver.

	A concrete implementation of this class is necessary.  The user can create
	their own implementation if those supplied are not suitable for their needs.

	\author Michael Heroux, Rich Lehoucq, Heidi Thornquist
*/

namespace Belos {

template <class ScalarType>
class MultiVec {
public:
  //! @name Constructor/Destructor
	//@{ 
	//! %Belos::MultiVec constructor.
	MultiVec() {};

	//! %Belos::MultiVec destructor.
	virtual ~MultiVec () {};

	//@}
  //! @name Creation methods for new multivectors
	//@{ 

	/*! \brief Creates a new empty %Belos::MultiVec containing \c numvecs columns.

	    \return Pointer to the new multivector with uninitialized values	
	*/

	virtual MultiVec<ScalarType> * Clone ( const int numvecs ) const = 0;

	/*! \brief Creates a new %Belos::MultiVec and copies contents of \c *this into
	    the new vector (deep copy).
	
	    \return Pointer to the new multivector	
	*/
	
	virtual MultiVec<ScalarType> * CloneCopy () const = 0;
	
	/*! \brief Creates a new %Belos::MultiVec and copies the selected contents of \c *this 
	    into the new vector (deep copy).  The copied 
	    vectors from \c *this are indicated by the \c index.size() indices in \c index.

	    \return Pointer to the new multivector	
	*/

	virtual MultiVec<ScalarType> * CloneCopy ( const std::vector<int>& index ) const = 0;
	
	/*! \brief Creates a new %Belos::MultiVec that shares the selected contents of \c *this.
	    The index of the \c numvecs vectors copied from \c *this are indicated by the
	    indices given in \c index.

	    \return Pointer to the new multivector	
	*/

	virtual MultiVec<ScalarType> * CloneView ( const std::vector<int>& index ) = 0;
	//@}

  //! @name Dimension information methods	
	//@{ 
	//! Obtain the vector length of *this multivector block.

	virtual int GetVecLength () const = 0;

	//! Obtain the number of vectors in *this multivector block.

	virtual int GetNumberVecs () const = 0;

	//@}
  //! @name Update methods
	//@{ 
	/*! \brief Update \c *this with \c alpha * \c A * \c B + \c beta * (\c *this).
	*/

	virtual void MvTimesMatAddMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, 
		const Teuchos::SerialDenseMatrix<int,ScalarType>& B, const ScalarType beta ) = 0;

	/*! \brief Replace \c *this with \c alpha * \c A + \c beta * \c B.
	*/

	virtual void MvAddMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, const ScalarType beta, const MultiVec<ScalarType>& B ) = 0;

	/*! \brief Compute a dense matrix \c B through the matrix-matrix multiply 
	   \c alpha * \c A^T * (\c *this).
	*/

	virtual void MvTransMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, Teuchos::SerialDenseMatrix<int,ScalarType>& B) const = 0;

	/*! \brief Compute a vector \c b where the components are the individual dot-products, i.e.\c b[i] = \c A[i]^T*\c this[i] where \c A[i] is the i-th column of A.
	*/

	virtual void MvDot ( const MultiVec<ScalarType>& A, std::vector<ScalarType>* b ) const = 0;

	//@}
  //! @name Norm method
	//@{ 

	/*! \brief Compute the 2-norm of each individual vector of \c *this.  
	   Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
	*/

        virtual void MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>* normvec, NormType type = TwoNorm ) const = 0;

	//@}
  //! @name Initialization methods
	//@{ 
	/*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  The \c 
  	    numvecs vectors in \c A are copied to a subset of vectors in \c *this
	    indicated by the indices given in \c index.
	*/

	virtual void SetBlock ( const MultiVec<ScalarType>& A, const std::vector<int>& index ) = 0;
	
	/*! \brief Replace the vectors in \c *this with random vectors.
	*/

	virtual void MvRandom () = 0;

	/*! \brief Replace each element of the vectors in \c *this with \c alpha.
	*/

	virtual void MvInit ( const ScalarType alpha ) = 0;

	//@}
  //! @name Print method
	//@{ 
	/*! \brief Print the \c *this multivector.
	*/
	virtual void MvPrint ( ostream& os ) const = 0;
	//@}
};


  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for Belos::MultiVec.
  //
  ////////////////////////////////////////////////////////////////////


  template<class ScalarType>
  class MultiVecTraits<ScalarType,MultiVec<ScalarType> >
  {
  public:
    ///
    static Teuchos::RefCountPtr<MultiVec<ScalarType> > Clone( const MultiVec<ScalarType>& mv, const int numvecs )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).Clone(numvecs) ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<ScalarType> > CloneCopy( const MultiVec<ScalarType>& mv )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneCopy() ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<ScalarType> > CloneCopy( const MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneCopy(index) ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<ScalarType> > CloneView( MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( mv.CloneView(index) ); }
    ///
    static Teuchos::RefCountPtr<const MultiVec<ScalarType> > CloneView( const MultiVec<ScalarType>& mv, const std::vector<int>& index )
    { return Teuchos::rcp( const_cast<MultiVec<ScalarType>&>(mv).CloneView(index) ); }
    ///
    static int GetVecLength( const MultiVec<ScalarType>& mv )
    { return mv.GetVecLength(); }
    ///
    static int GetNumberVecs( const MultiVec<ScalarType>& mv )
    { return mv.GetNumberVecs(); }
    ///
    static void MvTimesMatAddMv( ScalarType alpha, const MultiVec<ScalarType>& A, 
				 const Teuchos::SerialDenseMatrix<int,ScalarType>& B, 
				 ScalarType beta, MultiVec<ScalarType>& mv )
    { mv.MvTimesMatAddMv(alpha, A, B, beta); }
    ///
    static void MvAddMv( ScalarType alpha, const MultiVec<ScalarType>& A, ScalarType beta, const MultiVec<ScalarType>& B, MultiVec<ScalarType>& mv )
    { mv.MvAddMv(alpha, A, beta, B); }
    ///
    static void MvTransMv( ScalarType alpha, const MultiVec<ScalarType>& A, const MultiVec<ScalarType>& mv, Teuchos::SerialDenseMatrix<int,ScalarType>& B )
    { mv.MvTransMv(alpha, A, B); }
    ///
    static void MvDot( const MultiVec<ScalarType>& mv, const MultiVec<ScalarType>& A, std::vector<ScalarType>* b )
    { mv.MvDot( A, b ); }
    ///
    static void MvNorm( const MultiVec<ScalarType>& mv, std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>* normvec, NormType type = TwoNorm )
    { mv.MvNorm(normvec,type); }
    ///
    static void SetBlock( const MultiVec<ScalarType>& A, const std::vector<int>& index, MultiVec<ScalarType>& mv )
    { mv.SetBlock(A, index); }
    ///
    static void MvRandom( MultiVec<ScalarType>& mv )
    { mv.MvRandom(); }
    ///
    static void MvInit( MultiVec<ScalarType>& mv, ScalarType alpha = Teuchos::ScalarTraits<ScalarType>::zero() )
    { mv.MvInit(alpha); }
    ///
    static void MvPrint( const MultiVec<ScalarType>& mv, ostream& os )
    { mv.MvPrint(os); }
    
  };


} // namespace Belos

#endif

// end of file BelosMultiVec.hpp
