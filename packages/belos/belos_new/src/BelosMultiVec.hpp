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

/*! \file BelosMultiVec.hpp
    \brief Virtual base class which defines the multivector interface required 
	by the iterative linear solver.
*/

#include "Teuchos_SerialDenseMatrix.hpp"
#include "BelosReturnType.hpp"
#include "BelosConfigDefs.hpp"

/*! 	\class Belos::MultiVec

	\brief Belos's templated pure virtual class for constructing multivectors that 
	are used by the eigensolver.

	A concrete implementation of this class is necessary.  The user can create
	their own implementation if those supplied are not suitable for their needs.

	\author Michael Heroux, Rich Lehoucq, and Heidi Thornquist
*/

namespace Belos {

  //! \enum Enumerated list for describing the multivector norm type.
  enum NormType {   OneNorm,       /*!< Compute the one-norm \f$\sum_{i=1}^{n}(|x_i w_i|)\f$ for each vector. */
		    TwoNorm,       /*!< Compute the two-norm *\f$\sqrt(\sum_{i=1}^{n}((x_i w_i)^2)\f$ for each vector. */
		    InfNorm,      /*!< Compute the infinity-norm \f$(\max_{i=1}^{n}\{|x_i w_i|\})\f$ for each vector. */
  };

template <class TYPE>
class MultiVec {
public:
	//@{ \name Constructor/Destructor

	//! Default Constructor.
	MultiVec() {};

	//! Destructor.
	virtual ~MultiVec () {};

	//@}

	//@{ \name Creation methods

	/*! \brief Creates a new empty %Belos::MultiVec containing \c numvecs columns.

	    \return Pointer to the new MultiVec	
	*/
	virtual MultiVec<TYPE> * Clone ( const int numvecs ) = 0;

	/*! \brief Creates a new %Belos::MultiVec and copies contents of \c *this into
	    the new vector (deep copy).
	
	    \return Pointer to the new MultiVec	
	*/
	virtual MultiVec<TYPE> * CloneCopy () = 0;
	
	/*! \brief Creates a new %Belos::MultiVec and copies the selected contents of \c *this 
	    into the new vector (deep copy).  The number (\c numvecs) of copied 
	    vectors from \c *this are indicated by the indices in \c index.

	    \return Pointer to the new MultiVec	
	*/
	virtual MultiVec<TYPE> * CloneCopy ( int index[], int numvecs ) = 0;
	
	/*! \brief Creates a new %Belos::MultiVec that shares the selected contents of \c *this (shallow copy).
	    The index of the \c numvecs vectors copied from \c *this are indicated by the
	    indices given in \c index.

	    \return Pointer to the new MultiVec	
	*/
	virtual MultiVec<TYPE> * CloneView ( int index[], int numvecs ) = 0;

	//@}

	//@{ \name  Accessor methods	

	//! Obtain the vector length of *this multivector block.

	virtual int GetVecLength () const = 0;

	//! Obtain the number of vectors in *this multivector block.

	virtual int GetNumberVecs () const = 0;

	//@}
	//@{ \name Update methods

	/*! \brief Update \c *this with \c alpha * \c A * \c B + \c beta * (\c *this).
	 */
	virtual void MvTimesMatAddMv ( TYPE alpha, MultiVec<TYPE>& A, 
		Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta ) = 0;

        /*! \brief Replace \c *this with \c alpha * \c A + \c beta * \c B.
	 */
	virtual void MvAddMv ( TYPE alpha, MultiVec<TYPE>& A, TYPE beta, MultiVec<TYPE>& B ) = 0;

	/*! \brief Compute a dense matrix \c B through the matrix-matrix multiply 
	   \c alpha * \c A^T * (\c *this).
	*/

	virtual void MvTransMv ( TYPE alpha, MultiVec<TYPE>& A, Teuchos::SerialDenseMatrix<int,TYPE>& B) = 0;

	//@}
	//@{ \name Norm method

	/*! \brief Compute the norm of each individual vector of \c *this.  

	   Upon return, \c normvec[i] holds the given norm of the \c i-th vector of \c *this
	   
	   \note <b> The two-norm must be supported by any implementation of this virtual base class.</b>
	   The other norm implementations should be implemented to allow more flexibility in 
	   designing status tests for the iterative solvers.
	*/
	virtual ReturnType MvNorm ( TYPE *normvec, NormType norm_type = TwoNorm ) = 0;

	//@}
	//@{ \name Initialization methods

	/*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  The \c 
  	    numvecs vectors in \c A are copied to a subset of vectors in \c *this
	    indicated by the indices given in \c index.
	*/
	virtual void SetBlock ( MultiVec<TYPE>& A, int index[], int numvecs ) = 0;
	
	/*! \brief Replace the vectors in \c *this with random vectors.
	*/
	virtual void MvRandom () = 0;

	/*! \brief Replace each element of the vectors in \c *this with \c alpha.
	*/
	virtual void MvInit ( TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() ) = 0;

	//@}
	//@{ \name Print method.

	/*! \brief Print the \c *this multivector.  

	    \note This does not have to be implemented, a general statement will be sent to the output stream if it is not.
	*/
        virtual void MvPrint (ostream& os) { os << "Belos::MultiVec::MvPrint() is not supported." <<endl; };
	//@}
};

}
#endif
// end of file BelosMultiVec.hpp
